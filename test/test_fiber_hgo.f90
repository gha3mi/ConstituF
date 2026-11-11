! SPDX-License-Identifier: Apache-2.0

program test_fiber_hgo
    !! Test convention: stress and stress-like material parameters in Pa.
    !! Strain and deformation gradients are dimensionless. Explicit rescaling
    !! cases state their factors.
    !!
    !! Two-fiber tension-only energy, first Piola consistency,
    !! frame indifference, exact tangent and parameter derivatives.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
    use constituf, only: wp, material_t, verify_tangent, cf_ok, cf_invalid_argument, cf_evaluation_failed
    use constituf_fiber_hgo, only: fiber_hgo_t
    use constituf_neo_hookean, only: neo_hookean_t
    implicit none (type, external)

    integer, parameter :: ndim = 3, nfull = 9, noutput = 10, npoints = 2, nparameters = 4
    real(wp), parameter :: mu = 6.0_wp, lam = 8.0_wp, k1 = 4.0_wp, k2 = 3.0_wp
    real(wp), parameter :: base(nparameters) = [mu,lam,k1,k2]
    integer, parameter :: parameter_ids(nparameters) = [4,1,3,2]
    character(len=14), parameter :: names(nparameters) = [character(len=14) :: &
        "fiber_k2", "shear_modulus", "fiber_k1", "lame_lambda"]
    type(material_t) :: material, reference, perturbed
    real(wp) :: x(nfull,npoints), y(noutput,npoints), tangent(noutput,nfull,npoints)
    real(wp) :: dsdp(noutput,nparameters,npoints), p(nparameters)
    real(wp) :: h, err
    real(wp) :: plus(noutput,npoints), minus(noutput,npoints), elastic(nfull,npoints)
    real(wp) :: coarse(noutput,npoints), derivative(noutput,npoints), tolerance(noutput,npoints)
    real(wp) :: rotated_f(ndim,ndim), rotated_p(ndim,ndim)
    real(wp) :: rotation(ndim,ndim), first_f(ndim,ndim), first_p(ndim,ndim)
    real(wp) :: image(nfull,1), rotated_y(noutput,1)
    real(wp) :: directions(ndim,2)
    integer :: status, k, refinement, refinements

    x = 0.0_wp
    x(1,1) = 1.18_wp
    x(5,1) = 0.93_wp
    x(9,1) = 1.04_wp
    x(1,2) = 0.96_wp
    x(5,2) = 1.14_wp
    x(9,2) = 1.01_wp
    call material%init(fiber_hgo_t(mu,lam,k1,k2),npoints,status,sensitivity_parameters=names)
    if (status /= cf_ok) error stop "HGO initialization"
    if (material%input_size() /= nfull .or. material%output_size() /= noutput) error stop "HGO field dimensions"
    if (material%output_offset("energy.helmholtz") /= noutput) error stop "HGO energy field metadata"
    call material%eval(x,y,tangent,status,dsdp=dsdp)
    if (status /= cf_ok) error stop "Fiber HGO response"
    if (any(y(noutput,:) <= 0.0_wp)) error stop "Fiber HGO stored energies not positive"
    if (maxval(abs(tangent(noutput,:,1)-y(1:nfull,1))) > 1.0e-12_wp) then
        error stop "Fiber HGO energy gradient not Piola stress"
    end if

    call reference%init(neo_hookean_t(lam,mu),npoints,status)
    if (status /= cf_ok) error stop "HGO matrix-only reference setup"
    call reference%eval(x,elastic,status=status)
    if (status /= cf_ok) error stop "HGO matrix-only stress"
    if (y(1,1) <= elastic(1,1)) error stop "First fiber did not stiffen tension"
    if (y(5,2) <= elastic(5,2)) error stop "Second fiber did not stiffen tension"

    call perturbed%init(fiber_hgo_t(mu,lam,0.0_wp,k2),npoints,status)
    if (status /= cf_ok) error stop "Zero-fiber stiffness setup"
    call perturbed%eval(x,plus,status=status)
    if (status /= cf_ok .or. maxval(abs(plus(1:nfull,:)-elastic)) > 1.0e-10_wp) then
        error stop "Zero-fiber HGO does not reduce to Neo-Hookean"
    end if

    ! These stretches stay on their activation branches under the precision-scaled steps.
    call verify_tangent(material,x,err,status)
    if (status /= cf_ok .or. .not. ieee_is_finite(err) .or. err > 1.0_wp) then
        error stop "HGO material and energy tangent"
    end if

    do k = 1,nparameters
        h = 0.25_wp*base(parameter_ids(k))
        refinements = 1
        if (parameter_ids(k) == 4) then
            h = epsilon(1.0_wp)**(1.0_wp/3.0_wp)*base(4)
            refinements = 2
        end if
        ! Only k2 is nonlinear. The other parameter secants have no truncation error.
        do refinement = 1, refinements
            p = base
            p(parameter_ids(k)) = p(parameter_ids(k))+h
            call perturbed%init(fiber_hgo_t(p(1),p(2),p(3),p(4)),npoints,status)
            if (status /= cf_ok) error stop "HGO plus parameter setup"
            call perturbed%eval(x,plus,status=status)
            if (status /= cf_ok) error stop "HGO positive parameter perturbation"
            p(parameter_ids(k)) = base(parameter_ids(k))-h
            call perturbed%init(fiber_hgo_t(p(1),p(2),p(3),p(4)),npoints,status)
            if (status /= cf_ok) error stop "HGO minus parameter setup"
            call perturbed%eval(x,minus,status=status)
            if (status /= cf_ok) error stop "HGO negative parameter perturbation"
            if (refinement == 1) then
                coarse = (plus-minus)/(2.0_wp*h)
                if (refinements == 2) h = 0.5_wp*h
            end if
        end do
        derivative = (plus-minus)/(2.0_wp*h)
        tolerance = abs(derivative-coarse)/3.0_wp + &
            64.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(plus),abs(minus))/h
        if (.not. all(ieee_is_finite(derivative)) .or. .not. all(ieee_is_finite(coarse)) .or. &
            any(abs(derivative-dsdp(:,k,:)) > tolerance)) then
            error stop "Fiber HGO analytical parameter derivative"
        end if
    end do

    rotation = reshape([2.0_wp,1.0_wp,-2.0_wp,1.0_wp,2.0_wp,2.0_wp, &
                        2.0_wp,-2.0_wp,1.0_wp],[ndim,ndim])/3.0_wp
    first_f = reshape(x(:,1),[ndim,ndim])
    first_p = reshape(y(1:nfull,1),[ndim,ndim])
    rotated_f = matmul(rotation,first_f)
    image(:,1) = reshape(rotated_f,[nfull])
    call perturbed%init(fiber_hgo_t(mu,lam,k1,k2),1,status)
    if (status /= cf_ok) error stop "Frame-rotation HGO setup"
    call perturbed%eval(image,rotated_y,status=status)
    if (status /= cf_ok) error stop "Rotated HGO evaluation"
    rotated_p = matmul(rotation,first_p)
    if (maxval(abs(rotated_y(1:nfull,1)-reshape(rotated_p,[nfull]))) > &
        512.0_wp*epsilon(1.0_wp)*max(1.0_wp,maxval(abs(rotated_p)))) then
        error stop "Fiber HGO violates frame indifference"
    end if
    if (abs(rotated_y(noutput,1)-y(noutput,1)) > &
        512.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(y(noutput,1)))) error stop "Rotated HGO changed energy"

    call test_physical_fibers()
    call test_fiber_energy_limits()

    x = 0.0_wp
    x(1,:) = 1.0_wp
    x(5,:) = 1.0_wp
    x(9,:) = 1.0_wp
    call material%eval(x,y,status=status,dsdp=dsdp)
    if (status /= cf_ok .or. maxval(abs(y)) > 1.0e-10_wp) error stop "HGO reference configuration not stress-free"

    directions = reshape([2.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp,0.0_wp],[ndim,2])
    call perturbed%init(fiber_hgo_t(mu,lam,k1,k2,directions),1,status)
    if (status /= cf_invalid_argument) error stop "Nonunit HGO fiber direction accepted"
    call perturbed%init(fiber_hgo_t(mu,lam,k1,0.0_wp),1,status)
    if (status /= cf_invalid_argument) error stop "Nonpositive HGO exponential exponent accepted"

    x(1,2) = -1.0_wp
    call material%eval(x,y,status=status,dsdp=dsdp)
    if (status /= cf_evaluation_failed) error stop "HGO accepted a negative-J deformation"
    x(1,2) = ieee_value(0.0_wp,ieee_quiet_nan)
    call material%eval(x,y,status=status,dsdp=dsdp)
    if (status /= cf_evaluation_failed) error stop "HGO accepted a nonfinite deformation"

contains

    subroutine test_physical_fibers()
        type(material_t) :: fiber, matrix
        real(wp) :: driving(nfull,2), actual(noutput,2), jacobian(noutput,nfull,2)
        real(wp) :: matrix_stress(nfull,2), matrix_tangent(nfull,nfull,2), sensitivities(noutput,2,2)
        real(wp) :: axes(ndim,2), expected(noutput)
        real(wp) :: log_j, matrix_energy
        integer :: result, point

        call fiber%init(fiber_hgo_t(mu,lam,k1,k2),2,result, &
            sensitivity_parameters=[character(len=8) :: "fiber_k1","fiber_k2"])
        if (result /= cf_ok) error stop "Fiber branch setup"
        call matrix%init(neo_hookean_t(lam,mu),2,result)
        if (result /= cf_ok) error stop "Fiber branch matrix setup"
        driving = 0.0_wp
        driving([1,5,9],1) = 1.0_wp
        driving([1,5,9],2) = [0.75_wp,0.875_wp,1.0_wp]
        call fiber%eval(driving,actual,jacobian,result,dsdp=sensitivities)
        if (result /= cf_ok) error stop "Inactive fiber branch evaluation"
        call matrix%eval(driving,matrix_stress,matrix_tangent,result)
        if (result /= cf_ok) error stop "Inactive fiber matrix reference"
        if (maxval(abs(actual(1:nfull,:)-matrix_stress)) > &
            512.0_wp*epsilon(1.0_wp)*max(1.0_wp,maxval(abs(matrix_stress)))) then
            error stop "Compressed fibers produced stress"
        end if
        ! At I4 = 1 compare the documented inactive tangent without crossing the switch.
        if (maxval(abs(jacobian(1:nfull,:,:)-matrix_tangent)) > &
            512.0_wp*epsilon(1.0_wp)*max(1.0_wp,maxval(abs(matrix_tangent)))) then
            error stop "Inactive fiber tangent differs from matrix"
        end if
        if (.not. all(abs(sensitivities) <= 0.0_wp)) error stop "Inactive fiber parameter derivative"
        do point = 1, 2
            log_j = log(product(driving([1,5,9],point)))
            matrix_energy = 0.5_wp*mu*(sum(driving(:,point)**2)-3.0_wp-2.0_wp*log_j) + &
                            0.5_wp*lam*log_j**2
            if (abs(actual(noutput,point)-matrix_energy) > &
                512.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(matrix_energy))) then
                error stop "Inactive fibers produced stored energy"
            end if
        end do

        ! One off-axis fiber stretches with I4-1 = 9/20. The e3 family is compressed.
        axes = reshape([0.6_wp,0.8_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp],[ndim,2])
        call fiber%init(fiber_hgo_t(mu,lam,k1,k2,axes),1,result)
        if (result /= cf_ok) error stop "Off-axis fiber setup"
        call matrix%init(neo_hookean_t(lam,mu),1,result)
        if (result /= cf_ok) error stop "Off-axis matrix setup"
        driving(:,1) = 0.0_wp
        driving([1,5,9],1) = [1.5_wp,1.0_wp,0.75_wp]
        call fiber%eval(driving(:,1:1),actual(:,1:1),status=result)
        if (result /= cf_ok) error stop "Off-axis fiber evaluation"
        call matrix%eval(driving(:,1:1),matrix_stress(:,1:1),status=result)
        if (result /= cf_ok) error stop "Off-axis matrix reference"
        expected(1:nfull) = matrix_stress(:,1) + exp(0.6075_wp)* &
            [1.944_wp,1.728_wp,0.0_wp,2.592_wp,2.304_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp]
        log_j = log(9.0_wp/8.0_wp)
        expected(noutput) = 0.5_wp*mu*(13.0_wp/16.0_wp-2.0_wp*log_j) + &
            0.5_wp*lam*log_j**2 + (2.0_wp/3.0_wp)*(exp(0.6075_wp)-1.0_wp)
        if (maxval(abs(actual(:,1)-expected)) > &
            1024.0_wp*epsilon(1.0_wp)*max(1.0_wp,maxval(abs(expected)))) then
            error stop "Independent off-axis fiber stress or energy"
        end if
    end subroutine test_physical_fibers

    subroutine test_fiber_energy_limits()
        type(material_t) :: fiber, matrix
        real(wp) :: driving(nfull,1), actual(noutput,1), jacobian(noutput,nfull,1)
        real(wp) :: sensitivities(noutput,2,1), matrix_stress(nfull,1), matrix_tangent(nfull,nfull,1)
        real(wp) :: curvature, log_j, matrix_energy
        integer :: result

        curvature = epsilon(1.0_wp)/16.0_wp
        call fiber%init(fiber_hgo_t(mu,lam,k1,curvature),1,result, &
            sensitivity_parameters=[character(len=8) :: "fiber_k1","fiber_k2"])
        if (result /= cf_ok) error stop "Quadratic fiber limit setup"
        driving = 0.0_wp
        driving([1,5,9],1) = [1.25_wp,1.0_wp,1.0_wp]
        call fiber%eval(driving,actual,status=result,dsdp=sensitivities)
        if (result /= cf_ok) error stop "Quadratic fiber limit evaluation"
        log_j = log(1.25_wp)
        matrix_energy = 0.5_wp*mu*(9.0_wp/16.0_wp-2.0_wp*log_j)+0.5_wp*lam*log_j**2
        ! Activation = 9/16. Corrections to these k2 -> 0 limits are below roundoff.
        if (abs(actual(noutput,1)-matrix_energy-81.0_wp/128.0_wp) > &
            1024.0_wp*epsilon(1.0_wp)) error stop "Small-k2 fiber energy cancellation"
        if (abs(sensitivities(noutput,1,1)-81.0_wp/512.0_wp) > &
            1024.0_wp*epsilon(1.0_wp)) error stop "Small-k2 fiber stiffness derivative"
        if (abs(sensitivities(noutput,2,1)-6561.0_wp/65536.0_wp) > &
            1024.0_wp*epsilon(1.0_wp)) error stop "Small-k2 fiber curvature derivative"

        ! Absent fibers need no exponential, even if their unused exponent overflows.
        curvature = 2.0_wp*log(huge(1.0_wp))/9.0_wp
        call fiber%init(fiber_hgo_t(mu,lam,0.0_wp,curvature),1,result, &
            sensitivity_parameters=[character(len=8) :: "fiber_k2"])
        if (result /= cf_ok) error stop "Zero-fiber limit setup"
        driving(1,1) = 2.0_wp
        call fiber%eval(driving,actual,jacobian,result,dsdp=sensitivities(:,1:1,:))
        if (result /= cf_ok) error stop "Absent fiber exponential affected matrix limit"
        call matrix%init(neo_hookean_t(lam,mu),1,result)
        if (result /= cf_ok) error stop "Zero-fiber matrix setup"
        call matrix%eval(driving,matrix_stress,matrix_tangent,result)
        if (result /= cf_ok) error stop "Zero-fiber matrix reference"
        if (maxval(abs(actual(1:nfull,:)-matrix_stress)) > &
            512.0_wp*epsilon(1.0_wp)*maxval(abs(matrix_stress))) error stop "Zero-fiber stress limit"
        if (maxval(abs(jacobian(1:nfull,:,:)-matrix_tangent)) > &
            512.0_wp*epsilon(1.0_wp)*maxval(abs(matrix_tangent))) error stop "Zero-fiber tangent limit"
        log_j = log(2.0_wp)
        matrix_energy = 0.5_wp*mu*(3.0_wp-2.0_wp*log_j)+0.5_wp*lam*log_j**2
        if (abs(actual(noutput,1)-matrix_energy) > &
            512.0_wp*epsilon(1.0_wp)*abs(matrix_energy)) error stop "Zero-fiber energy limit"
        if (.not. all(abs(sensitivities(:,1,1)) <= 0.0_wp)) error stop "Zero-fiber curvature derivative"

        ! The derivative with respect to k1 still requires that unrepresentable exponential.
        call fiber%init(fiber_hgo_t(mu,lam,0.0_wp,curvature),1,result, &
            sensitivity_parameters=[character(len=8) :: "fiber_k1"])
        if (result /= cf_ok) error stop "Zero-fiber stiffness derivative setup"
        call fiber%eval(driving,actual,status=result,dsdp=sensitivities(:,1:1,:))
        if (result /= cf_evaluation_failed) error stop "Unrepresentable fiber stiffness derivative accepted"
    end subroutine test_fiber_energy_limits
end program test_fiber_hgo

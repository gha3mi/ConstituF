! SPDX-License-Identifier: Apache-2.0

program test_axisymmetric_elastic
    !! Test convention: stress and stress-like material parameters in Pa.
    !! Strain and deformation gradients are dimensionless. Explicit rescaling
    !! cases state their factors.
    !!
    !! Compare axisymmetric 4-component stress and sensitivities with
    !! independent 3D Hooke elastic values and the discrete tangent.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
    use constituf, only: wp, material_t, verify_tangent, cf_ok, cf_invalid_argument, cf_evaluation_failed
    use constituf_axisymmetric_elastic, only: axisymmetric_elastic_t
    use constituf_isotropic_elastic, only: isotropic_elastic_t
    implicit none (type, external)

    integer, parameter :: naxis = 4, nmandel = 6, npoints = 2
    real(wp), parameter :: bulk = 150000.0_wp, shear = 75000.0_wp
    real(wp), parameter :: reference_parameters(2) = [bulk,shear]
    integer, parameter :: mapping(naxis) = [1,2,3,6]
    type(material_t) :: axis, reference, perturbed
    real(wp) :: x(naxis,npoints), y(naxis,npoints), c(naxis,naxis,npoints)
    real(wp) :: x3(nmandel,npoints), y3(nmandel,npoints), c3(nmandel,nmandel,npoints)
    real(wp) :: dsdp(naxis,2,npoints), plus(naxis,npoints), minus(naxis,npoints)
    real(wp) :: changed(2), expected(naxis,npoints)
    real(wp) :: h, diff, error, tolerance
    integer :: status, k, i, j
    integer :: point_codes(npoints)

    call axis%init(axisymmetric_elastic_t(bulk,shear),npoints,status, &
        sensitivity_parameters=[character(len=13) :: "shear_modulus", "bulk_modulus"])
    if (status /= cf_ok) error stop "Axisymmetric initialization"
    if (.not. axis%field_matches("input","strain.small","custom.axisymmetric_mandel","small_strain","spatial",[3,3])) then
        error stop "Axisymmetric custom physical field specification"
    end if
    if (.not. axis%field_matches("output","stress.cauchy","custom.axisymmetric_mandel","cauchy_stress","spatial",[3,3])) then
        error stop "Axisymmetric stress field specification"
    end if
    if (axis%input_size() /= naxis .or. axis%output_size() /= naxis) error stop "Axisymmetric packed dimensions"

    x(:,1) = [0.006_wp,-0.002_wp,0.001_wp,0.0005_wp]
    x(:,2) = [-0.001_wp,0.003_wp,-0.004_wp,0.0009_wp]
    call axis%eval(x,y,c,status,dsdp=dsdp)
    if (status /= cf_ok) error stop "Axisymmetric response"
    if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(c)) .or. &
        .not. all(ieee_is_finite(dsdp))) error stop "Nonfinite axisymmetric result"

    call reference%init(isotropic_elastic_t(bulk,shear),npoints,status)
    if (status /= cf_ok) error stop "3D elastic reference initialization"
    x3 = 0.0_wp
    x3(mapping,:) = x
    call reference%eval(x3,y3,c3,status)
    if (status /= cf_ok) error stop "3D elastic reference"
    tolerance = 64.0_wp*epsilon(1.0_wp)*max(1.0_wp,maxval(abs(y3)))
    if (maxval(abs(y-y3(mapping,:))) > tolerance) error stop "Axisymmetric response differs from 3D"
    tolerance = 64.0_wp*epsilon(1.0_wp)*maxval(abs(c3))
    do j = 1, naxis
        do i = 1, naxis
            if (maxval(abs(c(i,j,:)-c3(mapping(i),mapping(j),:))) > tolerance) then
                error stop "Axisymmetric consistent 3D tangent mismatch"
            end if
        end do
    end do

    do k = 1, 2
        changed = reference_parameters
        ! Stress is affine in each modulus, so this finite secant has no truncation error.
        h = 0.25_wp*changed(3-k)
        changed(3-k) = changed(3-k)+h
        call perturbed%init(axisymmetric_elastic_t(changed(1),changed(2)),npoints,status)
        if (status /= cf_ok) error stop "Axisymmetric perturbed material"
        call perturbed%eval(x,plus,status=status)
        if (status /= cf_ok) error stop "Axisymmetric perturbed stress"
        changed(3-k) = changed(3-k)-2.0_wp*h
        call perturbed%init(axisymmetric_elastic_t(changed(1),changed(2)),npoints,status)
        if (status /= cf_ok) error stop "Axisymmetric negative perturbation setup"
        call perturbed%eval(x,minus,status=status)
        if (status /= cf_ok) error stop "Axisymmetric negative perturbation"
        diff = maxval(abs((plus-minus)/(2.0_wp*h)-dsdp(:,k,:)))
        tolerance = 64.0_wp*epsilon(1.0_wp)*max(1.0_wp,maxval(abs(plus)),maxval(abs(minus)))/h
        if (.not. ieee_is_finite(diff) .or. diff > tolerance) then
            error stop "Axisymmetric exact analytical parameter sensitivity"
        end if
    end do

    call verify_tangent(axis,x,error,status,step=[(0.001_wp,i=1,naxis)], &
                        relative_tolerance=256.0_wp*epsilon(1.0_wp))
    if (status /= cf_ok .or. .not. ieee_is_finite(error) .or. error > 1.0_wp) then
        error stop "Axisymmetric tangent verification"
    end if

    ! Isolate hoop extension and physical tensor shear epsilon_rz = 0.001.
    x(:,1) = [0.0_wp,0.0_wp,0.001_wp,0.0_wp]
    x(:,2) = [0.0_wp,0.0_wp,0.0_wp,sqrt(2.0_wp)*0.001_wp]
    expected(:,1) = [100.0_wp,100.0_wp,250.0_wp,0.0_wp]
    expected(:,2) = [0.0_wp,0.0_wp,0.0_wp,sqrt(2.0_wp)*150.0_wp]
    call axis%eval(x,y,status=status,dsdp=dsdp)
    if (status /= cf_ok) error stop "Axisymmetric physical reference evaluation"
    if (maxval(abs(y-expected)) > 64.0_wp*epsilon(1.0_wp)*maxval(abs(expected))) then
        error stop "Axisymmetric hoop or physical shear reference"
    end if

    ! The diagonal and shear stiffnesses are representable although 4*G overflows.
    call perturbed%init(axisymmetric_elastic_t(0.2_wp*huge(1.0_wp),0.3_wp*huge(1.0_wp)),npoints,status)
    if (status /= cf_ok) error stop "Representable axisymmetric stiffness rejected"
    x = 0.0_wp
    call perturbed%eval(x,y,c,status)
    if (status /= cf_ok .or. .not. all(ieee_is_finite(c))) error stop "Finite large axisymmetric stiffness"
    if (abs(c(1,1,1)/huge(1.0_wp)-0.6_wp) > 64.0_wp*epsilon(1.0_wp)) then
        error stop "Large axisymmetric normal stiffness"
    end if

    ! A small positive bulk modulus must survive the hydrostatic limit.
    call perturbed%init(axisymmetric_elastic_t(epsilon(1.0_wp)**2,1.0_wp),npoints,status)
    if (status /= cf_ok) error stop "Small-bulk axisymmetric initialization"
    x(1:3,:) = 1.0_wp/1024.0_wp
    call perturbed%eval(x,y,status=status)
    if (status /= cf_ok) error stop "Small-bulk axisymmetric evaluation"
    if (maxval(abs(y(1:3,:)/epsilon(1.0_wp)**2-3.0_wp/1024.0_wp)) > 64.0_wp*epsilon(1.0_wp)) then
        error stop "Axisymmetric hydrostatic response lost bulk modulus"
    end if

    ! Finite stresses do not make overflowed analytical sensitivities valid.
    call perturbed%init(axisymmetric_elastic_t(tiny(1.0_wp),tiny(1.0_wp)),npoints,status, &
        sensitivity_parameters=[character(len=13) :: "shear_modulus", "bulk_modulus"])
    if (status /= cf_ok) error stop "Tiny-modulus axisymmetric initialization"
    x = 0.0_wp
    x(:,1) = [0.75_wp*huge(1.0_wp),-0.75_wp*huge(1.0_wp),0.0_wp,0.75_wp*huge(1.0_wp)]
    call perturbed%eval(x,y,status=status,dsdp=dsdp)
    if (status /= cf_evaluation_failed) error stop "Axisymmetric accepted nonfinite sensitivity"
    call perturbed%point_statuses(point_codes,status)
    if (status /= cf_ok .or. any(point_codes /= [cf_invalid_argument,cf_ok])) then
        error stop "Axisymmetric sensitivity failure diagnostics"
    end if

    call perturbed%init(axisymmetric_elastic_t(bulk,0.0_wp),npoints,status)
    if (status /= cf_invalid_argument) error stop "Axisymmetric accepted zero shear"
    x = 0.0_wp
    x(2,2) = ieee_value(0.0_wp,ieee_quiet_nan)
    call axis%eval(x,y,status=status,dsdp=dsdp)
    if (status == cf_ok) error stop "Axisymmetric accepted nonfinite input"
end program test_axisymmetric_elastic

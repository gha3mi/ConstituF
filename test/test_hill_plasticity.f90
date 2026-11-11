! SPDX-License-Identifier: Apache-2.0

program test_hill_plasticity
    !! Test convention: stress and stress-like material parameters in Pa.
    !! Strain and deformation gradients are dimensionless. Explicit rescaling
    !! cases state their factors.
    !!
    !! Hill yield surface, isotropic J2 limit, return mapping and tangent.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf, only: wp, material_t, verify_tangent, cf_ok, cf_invalid_argument, cf_nondifferentiable
    use constituf_hill_plasticity, only: hill_plasticity_t
    use constituf_j2_plasticity, only: j2_plasticity_t
    implicit none (type, external)

    integer, parameter :: nmandel = 6, nstate = 7
    real(wp), parameter :: bulk = 150000.0_wp, shear = 75000.0_wp
    real(wp), parameter :: sy = 350.0_wp, hardening = 2000.0_wp
    real(wp), parameter :: reference_tolerance = 512.0_wp*epsilon(1.0_wp)
    real(wp), parameter :: anisotropy(nmandel) = [0.6_wp,0.4_wp,0.5_wp,1.6_wp,1.3_wp,1.8_wp]
    real(wp), parameter :: isotropic(nmandel) = [0.5_wp,0.5_wp,0.5_wp,1.5_wp,1.5_wp,1.5_wp]
    type(material_t) :: hill, j2, anisotropic_material, restored
    real(wp) :: x(nmandel,1), stress(nmandel,1), s_j2(nmandel,1), t(nmandel,nmandel,1), t_j2(nmandel,nmandel,1)
    real(wp) :: state(nstate,1), before(nstate,1), sigma(nmandel)
    real(wp) :: q, plastic_work, error, step_size, solve_scale, difference_tolerance
    integer :: status, step

    call test_nearly_incompressible_pressure()
    call test_pressure_and_directional_shear()
    call test_exact_shear(1.0_wp)
    call test_exact_shear(epsilon(1.0_wp)**2)
    call test_hydrostatic_pressure()
    call test_path_sensitivities()
    call hill%init(hill_plasticity_t(bulk,shear,sy,hardening,isotropic),1,status)
    if (status /= cf_ok) error stop "Hill isotropic initialization"
    call j2%init(j2_plasticity_t(bulk,shear,sy,hardening),1,status)
    if (status /= cf_ok) error stop "J2 reference initialization"
    do step = 1, 4
        x(:,1) = real(step,wp)*[0.004_wp,-0.002_wp,-0.002_wp,0.001_wp,0.0_wp,0.0005_wp]
        call hill%eval(x,stress,t,status)
        if (status /= cf_ok) error stop "Hill isotropic evaluation"
        call j2%eval(x,s_j2,t_j2,status)
        if (status /= cf_ok) error stop "J2 reference evaluation"
        if (.not. all(ieee_is_finite(stress)) .or. .not. all(ieee_is_finite(s_j2)) .or. &
            .not. all(ieee_is_finite(t)) .or. .not. all(ieee_is_finite(t_j2))) then
            error stop "Nonfinite Hill isotropic comparison"
        end if
        if (maxval(abs(stress-s_j2)) > reference_tolerance*max(sy,maxval(abs(s_j2)))) then
            error stop "Hill isotropic stress must reduce to J2"
        end if
        if (maxval(abs(t-t_j2)) > reference_tolerance*maxval(abs(t_j2))) then
            error stop "Hill isotropic tangent must reduce to J2"
        end if
        call hill%commit(status)
        if (status /= cf_ok) error stop "Hill isotropic commit"
        call j2%commit(status)
        if (status /= cf_ok) error stop "J2 reference commit"
    end do

    call anisotropic_material%init(hill_plasticity_t(bulk,shear,sy,hardening,anisotropy),1,status)
    if (status /= cf_ok) error stop "Anisotropic Hill initialization"
    x(:,1) = [0.013_wp,-0.006_wp,-0.007_wp,0.002_wp,0.0004_wp,-0.0009_wp]
    solve_scale = max(sy,hill_equivalent(2.0_wp*shear*x(:,1),anisotropy))
    call anisotropic_material%eval(x,stress,t,status)
    if (status /= cf_ok) error stop "Anisotropic Hill plastic return"
    call anisotropic_material%commit(status)
    if (status /= cf_ok) error stop "Anisotropic Hill commit"
    call anisotropic_material%export_committed(state,status)
    if (status /= cf_ok) error stop "Hill committed state export"
    if (.not. all(ieee_is_finite(state)) .or. .not. all(ieee_is_finite(stress))) then
        error stop "Nonfinite Hill plastic reference"
    end if
    if (state(7,1) <= 0.0_wp) error stop "Hill plasticity did not activate"
    if (abs(sum(state(1:3,1))) > reference_tolerance*maxval(abs(state(1:6,1)))) then
        error stop "Anisotropic flow must be isochoric"
    end if
    sigma = stress(:,1)
    q = hill_equivalent(sigma,anisotropy)
    if (.not. ieee_is_finite(q)) error stop "Nonfinite Hill equivalent stress"
    if (abs(q-(sy+hardening*state(7,1))) > 128.0_wp*epsilon(1.0_wp)*solve_scale) then
        error stop "Hill returned stress violates the yield surface"
    end if
    plastic_work = dot_product(sigma,state(1:6,1))
    if (.not. ieee_is_finite(plastic_work) .or. plastic_work < 0.0_wp) error stop "Negative Hill plastic dissipation"

    x(:,1) = [0.019_wp,-0.009_wp,-0.010_wp,0.003_wp,0.0007_wp,-0.0013_wp]
    step_size = maxval(abs(x))*epsilon(1.0_wp)**(1.0_wp/3.0_wp)
    solve_scale = max(sy+hardening*state(7,1), &
        hill_equivalent(2.0_wp*shear*(x(:,1)-state(1:6,1)),anisotropy))
    ! The return solve admits a 64*eps*trial-scale residual, amplified by the difference quotient.
    difference_tolerance = 192.0_wp*epsilon(1.0_wp)*solve_scale/step_size
    call verify_tangent(anisotropic_material,x,error,status, &
        step=[(step_size,step=1,6)],absolute_tolerance=difference_tolerance)
    if (status /= cf_ok .or. .not. ieee_is_finite(error) .or. error > 1.0_wp) then
        error stop "Hill consistent tangent failed"
    end if
    call anisotropic_material%export_committed(before,status)
    if (status /= cf_ok) error stop "Hill post-verification history read"
    if (.not. all(abs(before-state) <= 0.0_wp)) error stop "Tangent verification modified accepted history"

    call anisotropic_material%eval(x,stress,status=status)
    if (status /= cf_ok) error stop "Hill rejected trial evaluation"
    call anisotropic_material%revert()
    call anisotropic_material%export_committed(before,status)
    if (status /= cf_ok) error stop "Hill post-revert history read"
    if (.not. all(abs(before-state) <= 0.0_wp)) error stop "Hill revert changed committed plastic history"

    call restored%init(hill_plasticity_t(bulk,shear,sy,hardening,[anisotropy(1:5),0.0_wp]),1,status)
    if (status /= cf_invalid_argument) error stop "Invalid Hill yield coefficients accepted"
    call restored%init(hill_plasticity_t(bulk,shear,0.0_wp,hardening,anisotropy),1,status)
    if (status /= cf_invalid_argument) error stop "Zero Hill yield stress accepted"

contains

    subroutine test_pressure_and_directional_shear()
        integer, parameter :: npoints = 7, nselected = 7
        real(wp), parameter :: shear_ref = 100.0_wp, yield_ref = 0.1_wp, shear_input = 0.002_wp
        real(wp), parameter :: coefficients(nmandel) = [0.7_wp,0.4_wp,0.5_wp,1.6_wp,1.3_wp,1.8_wp]
        character(len=12), parameter :: names(nselected) = [character(len=12) :: &
            "bulk_modulus","hill_f","hill_g","hill_h","hill_l","hill_m","hill_n"]
        type(material_t) :: material
        real(wp) :: driving(nmandel,npoints), response(nmandel,npoints), jacobian(nmandel,nmandel,npoints)
        real(wp) :: gradient(nmandel,nselected,npoints), history(nstate,npoints)
        real(wp) :: history_gradient(nstate,nselected,npoints)
        real(wp) :: expected_response(nmandel,npoints), expected_gradient(nmandel,nselected,npoints)
        real(wp) :: expected_history(nstate,npoints), expected_history_gradient(nstate,nselected,npoints)
        real(wp) :: coefficient, root_coefficient
        integer :: result, axis, pressure_case, ipoint, component, parameter

        call material%init(hill_plasticity_t(1.0e6_wp,shear_ref,yield_ref,0.0_wp,coefficients), &
            npoints,result,sensitivity_parameters=names)
        if (result /= cf_ok) error stop "Hill pressure and directional reference setup"
        driving = 0.0_wp
        driving(1:3,1) = 0.001_wp
        driving(1:3,5:7) = 0.001_wp
        expected_response = 0.0_wp
        expected_response(1:3,1) = 3000.0_wp
        expected_response(1:3,5:7) = 3000.0_wp
        expected_gradient = 0.0_wp
        expected_gradient(1:3,1,1) = 0.003_wp
        expected_gradient(1:3,1,5:7) = 0.003_wp
        expected_history = 0.0_wp
        expected_history_gradient = 0.0_wp
        do axis = 4,6
            coefficient = coefficients(axis)
            root_coefficient = sqrt(coefficient)
            do pressure_case = 0,1
                ipoint = axis-2+3*pressure_case
                driving(axis,ipoint) = shear_input
                ! In perfect plasticity, pure Mandel shear returns to Y/sqrt(L), Y/sqrt(M) or Y/sqrt(N).
                expected_response(axis,ipoint) = yield_ref/root_coefficient
                expected_gradient(axis,axis+1,ipoint) = -yield_ref/(2.0_wp*coefficient*root_coefficient)
                expected_history(axis,ipoint) = shear_input-yield_ref/(2.0_wp*shear_ref*root_coefficient)
                expected_history(nstate,ipoint) = shear_input/root_coefficient - &
                    yield_ref/(2.0_wp*shear_ref*coefficient)
                expected_history_gradient(axis,axis+1,ipoint) = &
                    yield_ref/(4.0_wp*shear_ref*coefficient*root_coefficient)
                expected_history_gradient(nstate,axis+1,ipoint) = &
                    -shear_input/(2.0_wp*coefficient*root_coefficient) + &
                    yield_ref/(2.0_wp*shear_ref*coefficient**2)
            end do
        end do
        call material%eval(driving,response,jacobian,result,dsdp=gradient)
        if (result /= cf_ok) error stop "Hill hydrostatic pressure must not create plastic failure"
        do ipoint = 1,npoints
            do component = 1,nmandel
                call close(response(component,ipoint),expected_response(component,ipoint))
                do parameter = 1,nselected
                    call close(gradient(component,parameter,ipoint),expected_gradient(component,parameter,ipoint))
                end do
            end do
        end do
        do axis = 4,6
            call close(jacobian(axis,axis,axis-2),0.0_wp)
            call close(jacobian(axis,axis,axis+1),0.0_wp)
        end do
        call material%commit(result)
        if (result /= cf_ok) error stop "Hill pressure reference acceptance"
        call material%export_committed(history,result,sensitivity_state=history_gradient)
        if (result /= cf_ok) error stop "Hill pressure reference history export"
        do ipoint = 1,npoints
            do component = 1,nstate
                ! Express strain errors in stress units, matching the return solver residual.
                call close(2.0_wp*shear_ref*history(component,ipoint), &
                    2.0_wp*shear_ref*expected_history(component,ipoint))
                do parameter = 1,nselected
                    call close(2.0_wp*shear_ref*history_gradient(component,parameter,ipoint), &
                        2.0_wp*shear_ref*expected_history_gradient(component,parameter,ipoint))
                end do
            end do
        end do
    end subroutine test_pressure_and_directional_shear


    subroutine test_path_sensitivities()
        !! Perturb every parameter across the *same complete accepted strain
        !! history*, not just the final increment. This verifies analytical
        !! stress and plastic-history sensitivity propagation independently.
        integer, parameter :: nparameters = 10, npoints = 2, nsteps = 5
        character(len=18), parameter :: names(nparameters) = [character(len=18) :: &
            "hill_n", "yield_stress", "bulk_modulus", "hill_f", "shear_modulus", &
            "hill_g", "hardening_modulus", "hill_h", "hill_l", "hill_m"]
        real(wp), parameter :: factors(nsteps) = [0.10_wp,0.8_wp,1.3_wp,-0.9_wp,0.2_wp]
        real(wp), parameter :: reference_parameters(nparameters) = [bulk,shear,sy,hardening,anisotropy]
        type(material_t) :: tracked, ordinary, one_point, switch
        type(material_t) :: plus(nparameters,2), minus(nparameters,2)
        real(wp) :: xnew(nmandel,npoints), tracked_y(nmandel,npoints), ordinary_y(nmandel,npoints)
        real(wp) :: tangent(nmandel,nmandel,npoints), dsdp(nmandel,nparameters,npoints)
        real(wp) :: trial_sensitivity(nmandel,nparameters,npoints), trial_response(nmandel,npoints)
        real(wp) :: yp(nmandel,npoints), ym(nmandel,npoints)
        real(wp) :: plus_parameters(nparameters), minus_parameters(nparameters)
        real(wp) :: perturbations(nparameters,2), fd(nmandel,npoints,2), tolerance(nmandel,npoints)
        real(wp) :: accepted(nstate,npoints), accepted_sensitivity(nstate,nparameters,npoints)
        real(wp) :: accumulated_solve_scale(npoints)
        real(wp) :: one_x(nmandel,1), one_y(nmandel,1)
        integer :: stat, k, istep, id, refinement, ipoint

        call tracked%init(hill_plasticity_t(bulk,shear,sy,hardening,anisotropy),npoints,stat, &
            sensitivity_parameters=names)
        if (stat /= cf_ok) error stop "Hill analytical-sensitivity initialization"
        call ordinary%init(hill_plasticity_t(bulk,shear,sy,hardening,anisotropy),npoints,stat)
        if (stat /= cf_ok) error stop "Hill ordinary sensitivity reference"
        call one_point%init(hill_plasticity_t(bulk,shear,sy,hardening,anisotropy),1,stat)
        if (stat /= cf_ok) error stop "Hill single-point setup"

        do k = 1, nparameters
            select case (trim(names(k)))
            case ("bulk_modulus")
                id = 1
            case ("shear_modulus")
                id = 2
            case ("yield_stress")
                id = 3
            case ("hardening_modulus")
                id = 4
            case ("hill_f")
                id = 5
            case ("hill_g")
                id = 6
            case ("hill_h")
                id = 7
            case ("hill_l")
                id = 8
            case ("hill_m")
                id = 9
            case ("hill_n")
                id = 10
            case default
                error stop "Unrecognized test parameter"
            end select
            ! The solver residual is O(64*eps), so balance it with the O(h^2) secant error.
            perturbations(k,1) = (64.0_wp*epsilon(1.0_wp))**(1.0_wp/3.0_wp)* &
                abs(reference_parameters(id))
            perturbations(k,2) = 0.5_wp*perturbations(k,1)
            do refinement = 1,2
                plus_parameters = reference_parameters
                minus_parameters = reference_parameters
                plus_parameters(id) = plus_parameters(id)+perturbations(k,refinement)
                minus_parameters(id) = minus_parameters(id)-perturbations(k,refinement)
                call plus(k,refinement)%init(hill_plasticity_t(plus_parameters(1),plus_parameters(2),plus_parameters(3), &
                    plus_parameters(4),plus_parameters(5:10)),npoints,stat)
                if (stat /= cf_ok) error stop "Hill positive parameter perturbation"
                call minus(k,refinement)%init(hill_plasticity_t(minus_parameters(1),minus_parameters(2),minus_parameters(3), &
                    minus_parameters(4),minus_parameters(5:10)),npoints,stat)
                if (stat /= cf_ok) error stop "Hill negative parameter perturbation"
            end do
        end do

        accepted = 0.0_wp
        accumulated_solve_scale = 0.0_wp
        do istep = 1, nsteps
            xnew(:,1) = factors(istep)*[0.010_wp,-0.005_wp,-0.005_wp,0.001_wp,0.0008_wp,0.002_wp]
            xnew(:,2) = factors(istep)*[-0.004_wp,0.001_wp,0.003_wp,0.002_wp,0.001_wp,-0.001_wp]
            do ipoint = 1,npoints
                accumulated_solve_scale(ipoint) = accumulated_solve_scale(ipoint) + &
                    max(sy+hardening*accepted(7,ipoint), &
                    hill_equivalent(2.0_wp*shear*(xnew(:,ipoint)-accepted(1:6,ipoint)),anisotropy))
            end do
            call tracked%eval(xnew,tracked_y,tangent,stat,dsdp=dsdp)
            if (stat /= cf_ok) error stop "Hill tracked return through prescribed history"
            call ordinary%eval(xnew,ordinary_y,status=stat)
            if (stat /= cf_ok) error stop "Hill untracked return through prescribed history"
            if (.not. all(ieee_is_finite(tracked_y)) .or. .not. all(ieee_is_finite(ordinary_y)) .or. &
                .not. all(ieee_is_finite(dsdp))) error stop "Nonfinite Hill sensitivity reference"
            if (maxval(abs(tracked_y-ordinary_y)) > &
                reference_tolerance*max(sy,maxval(abs(ordinary_y)))) then
                error stop "Hill sensitivity path changed material stress"
            end if
            if (istep == 1) then
                call one_point%eval(xnew(:,1:1),one_y,status=stat)
                if (stat /= cf_ok) error stop "Hill single-point reference"
                if (.not. all(ieee_is_finite(one_y))) error stop "Nonfinite Hill single-point reference"
                if (maxval(abs(one_y(:,1)-tracked_y(:,1))) > &
                    reference_tolerance*max(sy,maxval(abs(one_y)))) then
                    error stop "Hill batch and single point differ"
                end if
            end if

            do k = 1, nparameters
                do refinement = 1,2
                    call plus(k,refinement)%eval(xnew,yp,status=stat)
                    if (stat /= cf_ok) error stop "Hill positive full-history evaluation"
                    call minus(k,refinement)%eval(xnew,ym,status=stat)
                    if (stat /= cf_ok) error stop "Hill negative full-history evaluation"
                    fd(:,:,refinement) = (yp-ym)/(2.0_wp*perturbations(k,refinement))
                end do
                ! Include accumulated solve residuals as well as the observed O(h^2) truncation.
                tolerance = abs(fd(:,:,2)-fd(:,:,1))/3.0_wp + &
                    192.0_wp*epsilon(1.0_wp)*spread(accumulated_solve_scale,1,nmandel)/perturbations(k,2)
                if (.not. all(ieee_is_finite(fd))) error stop "Nonfinite Hill parameter secant"
                if (any(abs(fd(:,:,2)-dsdp(:,k,:)) > tolerance)) then
                    print "(a,a,a,i0,a,es12.4)", "Parameter ",trim(names(k)), &
                        ", step ",istep,", discrepancy ",maxval(abs(fd(:,:,2)-dsdp(:,k,:)))
                    error stop "Hill exact path derivative disagrees with full-history perturbation"
                end if
            end do

            if (istep == 3) then
                trial_sensitivity = dsdp
                trial_response = tracked_y
                call tracked%revert()
                call tracked%eval(xnew,tracked_y,status=stat,dsdp=dsdp)
                if (stat /= cf_ok) error stop "Hill reverted derivative trial"
                if (.not. all(abs(dsdp-trial_sensitivity) <= &
                    reference_tolerance*max(1.0_wp,abs(trial_sensitivity)))) then
                    error stop "Hill rejected trial changed derivative history"
                end if
                if (.not. all(abs(tracked_y-trial_response) <= &
                    reference_tolerance*max(1.0_wp,abs(trial_response)))) then
                    error stop "Hill rejected trial changed stress history"
                end if
            end if

            call tracked%commit(stat)
            if (stat /= cf_ok) error stop "Hill tracked commit"
            call tracked%export_committed(accepted,stat,sensitivity_state=accepted_sensitivity)
            if (stat /= cf_ok) error stop "Hill accepted sensitivity history export"
            if (.not. all(ieee_is_finite(accepted)) .or. .not. all(ieee_is_finite(accepted_sensitivity))) then
                error stop "Nonfinite Hill accepted sensitivity history"
            end if
            call ordinary%commit(stat)
            if (stat /= cf_ok) error stop "Hill ordinary commit"
            if (istep == 1) then
                call one_point%commit(stat)
                if (stat /= cf_ok) error stop "Hill single-point commit"
            end if
            do k = 1, nparameters
                do refinement = 1,2
                    call plus(k,refinement)%commit(stat)
                    if (stat /= cf_ok) error stop "Hill perturbed accepted step"
                    call minus(k,refinement)%commit(stat)
                    if (stat /= cf_ok) error stop "Hill perturbed accepted step"
                end do
            end do
        end do

        call switch%init(hill_plasticity_t(bulk,shear,sy,hardening,anisotropy),1,stat, &
            sensitivity_parameters=[character(len=6) :: "hill_n"])
        if (stat /= cf_ok) error stop "Hill active-set setup"
        one_x = 0.0_wp
        one_x(nmandel,1) = sy/(2.0_wp*shear*sqrt(anisotropy(nmandel)))
        call switch%eval(one_x,one_y,status=stat,dsdp=trial_sensitivity(:,1:1,1:1))
        if (stat /= cf_nondifferentiable) error stop "Hill yield-surface sensitivity switch not reported"
    end subroutine test_path_sensitivities

    pure function hill_equivalent(value, coefficients) result(equivalent)
        real(wp), intent(in) :: value(nmandel), coefficients(nmandel)
        real(wp) :: equivalent

        equivalent = sqrt(coefficients(1)*(value(2)-value(3))**2 + &
            coefficients(2)*(value(3)-value(1))**2 + coefficients(3)*(value(1)-value(2))**2 + &
            sum(coefficients(4:6)*value(4:6)**2))
    end function hill_equivalent

    subroutine test_nearly_incompressible_pressure()
        !! A pressure-insensitive yield law must retain the virgin state under hydrostatic loading.
        type(material_t) :: material
        real(wp) :: driving(nmandel,1), response(nmandel,1), history(nstate,1)
        real(wp), parameter :: coefficients(nmandel) = [0.7_wp,0.4_wp,0.5_wp,1.6_wp,1.3_wp,1.8_wp]
        integer :: result

        call material%init(hill_plasticity_t(1.0e6_wp,100.0_wp,0.1_wp,0.0_wp,coefficients),1,result)
        if (result /= cf_ok) error stop "Hydrostatic Hill regression setup"
        driving = 0.0_wp
        driving(1:3,1) = 0.001_wp
        call material%eval(driving,response,status=result)
        if (result /= cf_ok) then
            print *, "Hydrostatic Hill status:",result
            error stop "Hill hydrostatic pressure must remain elastic"
        end if
        if (.not. all(ieee_is_finite(response))) error stop "Nonfinite hydrostatic Hill stress"
        if (maxval(abs(response(1:3,1)-3000.0_wp)) > 512.0_wp*epsilon(1.0_wp)*3000.0_wp) then
            error stop "Hill hydrostatic pressure response mismatch"
        end if
        call material%commit(result)
        if (result /= cf_ok) error stop "Hydrostatic Hill regression commit"
        call material%export_committed(history,result)
        if (result /= cf_ok) error stop "Hydrostatic Hill regression state"
        if (.not. all(ieee_is_finite(history)) .or. any(abs(history) > 0.0_wp)) then
            print *, "Hydrostatic Hill plastic history:",history(:,1)
            error stop "Hydrostatic pressure created Hill plastic flow"
        end if
    end subroutine test_nearly_incompressible_pressure

    subroutine test_exact_shear(modulus_scale)
        real(wp), intent(in) :: modulus_scale
        type(material_t) :: material
        real(wp), parameter :: coefficients(nmandel) = [0.6_wp,0.4_wp,0.5_wp,1.6_wp,1.3_wp,2.0_wp]
        real(wp) :: strain(nmandel,1), response(nmandel,1), jacobian(nmandel,nmandel,1)
        real(wp) :: derivative(nmandel,2,1), history(nstate,1), history_derivative(nstate,2,1)
        integer :: result, component

        ! Anisotropic pure Mandel shear has q=sqrt(N)*sigma6 and plastic flow6=sqrt(N).
        ! For N=2, q_trial=12, delta_alpha=3/85 and q_new=60/17 in reference stress units.
        call material%init(hill_plasticity_t(100.0_wp*modulus_scale,60.0_wp*modulus_scale, &
            3.0_wp*modulus_scale,15.0_wp*modulus_scale,coefficients),1,result, &
            sensitivity_parameters=[character(len=12) :: "yield_stress", "hill_n"])
        if (result /= cf_ok) error stop "Exact anisotropic Hill shear setup"
        strain = 0.0_wp
        strain(6,1) = 1.0_wp/(10.0_wp*sqrt(2.0_wp))
        call material%eval(strain,response,jacobian,result,dsdp=derivative)
        if (result /= cf_ok) error stop "Exact anisotropic Hill shear evaluation"
        call close(response(6,1)/modulus_scale,30.0_wp*sqrt(2.0_wp)/17.0_wp)
        call close(jacobian(6,6,1)/modulus_scale,120.0_wp/17.0_wp)
        call close(derivative(6,1,1),8.0_wp*sqrt(2.0_wp)/17.0_wp)
        call close(derivative(6,2,1)/modulus_scale,-138.0_wp*sqrt(2.0_wp)/289.0_wp)
        call material%commit(result)
        if (result /= cf_ok) error stop "Exact anisotropic Hill shear acceptance"
        call material%export_committed(history,result,sensitivity_state=history_derivative)
        if (result /= cf_ok) error stop "Exact anisotropic Hill shear history"
        do component = 1,5
            call close(response(component,1),0.0_wp)
            call close(history(component,1),0.0_wp)
            call close(derivative(component,1,1),0.0_wp)
            call close(derivative(component,2,1)/modulus_scale,0.0_wp)
            call close(history_derivative(component,1,1)*modulus_scale,0.0_wp)
            call close(history_derivative(component,2,1),0.0_wp)
        end do
        call close(history(6,1),3.0_wp*sqrt(2.0_wp)/85.0_wp)
        call close(history(7,1),3.0_wp/85.0_wp)
        call close(history_derivative(6,1,1)*modulus_scale,-sqrt(2.0_wp)/255.0_wp)
        call close(history_derivative(7,1,1)*modulus_scale,-1.0_wp/255.0_wp)
        call close(history_derivative(6,2,1),23.0_wp*sqrt(2.0_wp)/5780.0_wp)
        call close(history_derivative(7,2,1),-7.0_wp/1445.0_wp)
    end subroutine test_exact_shear

    subroutine test_hydrostatic_pressure()
        type(material_t) :: material
        real(wp), parameter :: large_bulk = 1.0_wp/epsilon(1.0_wp)
        real(wp) :: strain(nmandel,1), response(nmandel,1), derivative(nmandel,2,1)
        real(wp) :: history(nstate,1), history_derivative(nstate,2,1)
        integer :: result, component

        ! Hill's quadratic form annihilates pressure exactly, even when K greatly exceeds G.
        call material%init(hill_plasticity_t(large_bulk,60.0_wp,3.0_wp,15.0_wp,anisotropy),1,result, &
            sensitivity_parameters=[character(len=12) :: "bulk_modulus", "yield_stress"])
        if (result /= cf_ok) error stop "Hydrostatic Hill setup"
        strain = 0.0_wp
        strain(1:3,1) = 0.01_wp
        call material%eval(strain,response,status=result,dsdp=derivative)
        if (result /= cf_ok) error stop "Hill spuriously yields under hydrostatic pressure"
        do component = 1,3
            call close(response(component,1)/large_bulk,0.03_wp)
            call close(derivative(component,1,1),0.03_wp)
        end do
        do component = 4,nmandel
            call close(response(component,1),0.0_wp)
            call close(derivative(component,1,1),0.0_wp)
        end do
        do component = 1,nmandel
            call close(derivative(component,2,1),0.0_wp)
        end do
        call material%commit(result)
        if (result /= cf_ok) error stop "Hydrostatic Hill acceptance"
        call material%export_committed(history,result,sensitivity_state=history_derivative)
        if (result /= cf_ok) error stop "Hydrostatic Hill history"
        do component = 1,nstate
            call close(history(component,1),0.0_wp)
            call close(history_derivative(component,1,1),0.0_wp)
            call close(history_derivative(component,2,1),0.0_wp)
        end do
    end subroutine test_hydrostatic_pressure

    subroutine close(actual, expected)
        real(wp), intent(in) :: actual, expected

        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) then
            error stop "Nonfinite exact Hill comparison"
        end if
        if (abs(actual-expected) > reference_tolerance*max(1.0_wp,abs(expected))) then
            print *, "Expected",expected,"got",actual
            error stop "Exact anisotropic Hill shear reference mismatch"
        end if
    end subroutine close

end program test_hill_plasticity

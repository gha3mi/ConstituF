! SPDX-License-Identifier: Apache-2.0

program test_perzyna
    !! Test convention: moduli and stress in MPa, viscosity in MPa*s and
    !! dt in s. Explicit time-rescaling cases state their factors.
    !!
    !! Rate-dependent J2 radial return, incremental subdivision, exact
    !! endpoint tangent, state acceptance, rollback and error paths.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
    use constituf, only: wp, material_t, verify_tangent, cf_ok, cf_invalid_argument, cf_invalid_time, &
                        cf_evaluation_failed, cf_no_trial, cf_nondifferentiable
    use constituf_j2_plasticity, only: j2_plasticity_t
    use constituf_perzyna, only: perzyna_t
    implicit none (type, external)

    integer, parameter :: nmandel = 6, nstate = 7, npoints = 2, nsubsteps = 4
    real(wp), parameter :: bulk = 160000.0_wp, shear = 80000.0_wp
    real(wp), parameter :: sy = 350.0_wp, hardening = 2000.0_wp, eta = 5000.0_wp
    real(wp), parameter :: dt = 0.2_wp, upper = 0.003_wp
    type(material_t) :: material, steps, j2, copy
    real(wp) :: initial(nmandel,npoints), x(nmandel,npoints), previous(nmandel,npoints)
    real(wp) :: previous_steps(nmandel,npoints)
    real(wp) :: y(nmandel,npoints), y_steps(nmandel,npoints), tangent(nmandel,nmandel,npoints)
    real(wp) :: s(nstate,npoints), s_steps(nstate,npoints), s_before(nstate,npoints)
    real(wp) :: step_input(nmandel,npoints), normal(nmandel)
    real(wp) :: q, predicted
    real(wp) :: error, h, stress_scale, limit_tolerance
    integer :: status, step, i

    initial = 0.0_wp
    call material%init(perzyna_t(bulk,shear,sy,hardening,eta,0.0_wp),npoints,status)
    if (status /= cf_invalid_argument) error stop "Perzyna accepted absent initial strain"
    call material%init(perzyna_t(bulk,shear,sy,hardening,eta,upper,nsubsteps),npoints,status,initial, &
        sensitivity_parameters=[character(len=12) :: "bulk_modulus"])
    if (status /= cf_ok) error stop "Perzyna analytical sensitivity setup"
    call material%init(perzyna_t(bulk,shear,sy,hardening,eta,upper,nsubsteps),npoints,status,initial)
    if (status /= cf_ok) error stop "Perzyna initialization"

    x(:,1) = [0.010_wp,-0.004_wp,-0.006_wp,0.001_wp,0.0_wp,0.001_wp]
    x(:,2) = [-0.010_wp,0.003_wp,0.007_wp,0.001_wp,0.0_wp,0.0_wp]
    call material%eval(x,y,status=status)
    if (status /= cf_invalid_time) error stop "Perzyna allowed missing positive dt"
    call material%eval(x,y,status=status,dt=dt)
    if (status /= cf_ok) error stop "Perzyna subdivided response"
    call material%export_committed(s_before,status,previous_input=previous)
    if (status /= cf_ok .or. .not. all(ieee_is_finite(s_before)) .or. maxval(abs(s_before)) > 0.0_wp) then
        error stop "Perzyna trial changed committed history"
    end if
    if (.not. all(ieee_is_finite(previous)) .or. maxval(abs(previous-initial)) > 0.0_wp) then
        error stop "Perzyna did not retain initial strain"
    end if

    call steps%init(perzyna_t(bulk,shear,sy,hardening,eta),npoints,status,initial)
    if (status /= cf_ok) error stop "Perzyna manual substep setup"
    do step = 1, nsubsteps
        step_input = real(step,wp)/real(nsubsteps,wp)*x
        call steps%eval(step_input,y_steps,status=status,dt=dt/real(nsubsteps,wp))
        if (status /= cf_ok) error stop "Manual Perzyna substep failed"
        call steps%commit(status)
        if (status /= cf_ok) error stop "Manual Perzyna substep commit"
    end do
    stress_scale = 3.0_wp*shear*maxval(sqrt(sum(x*x,dim=1)))
    if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(y_steps)) .or. &
        maxval(abs(y-y_steps)) > 64.0_wp*epsilon(1.0_wp)*stress_scale) then
        error stop "Increment-controlled Perzyna differs from explicit substeps"
    end if

    call material%eval(x,y,tangent,status,dt=dt)
    if (status /= cf_ok) error stop "Perzyna exact tangent evaluation"
    call material%commit(status)
    if (status /= cf_ok) error stop "Perzyna commit"
    call material%export_committed(s,status,previous_input=previous)
    if (status /= cf_ok) error stop "Perzyna accepted history export"
    call steps%export_committed(s_steps,status,previous_input=previous_steps)
    if (status /= cf_ok) error stop "Manual Perzyna committed history export"
    if (.not. all(ieee_is_finite(previous_steps)) .or. maxval(abs(previous_steps-x)) > 0.0_wp) then
        error stop "Manual Perzyna previous strain mismatch"
    end if
    if (.not. all(ieee_is_finite(s)) .or. .not. all(ieee_is_finite(s_steps)) .or. &
        maxval(abs(s-s_steps)) > 64.0_wp*epsilon(1.0_wp)*maxval(abs(x))) then
        error stop "Substep plastic history mismatch"
    end if
    if (s(nstate,1) <= 0.0_wp) error stop "Perzyna plastic flow never activated"
    normal = x(:,1)
    normal(1:3) = normal(1:3)-sum(normal(1:3))/3.0_wp
    normal = 2.0_wp*shear*normal
    q = sqrt(1.5_wp*dot_product(normal,normal))
    predicted = q-sy
    if (predicted <= 0.0_wp) error stop "Perzyna trial must exceed yield"

    call copy%init(perzyna_t(bulk,shear,sy,hardening,eta,upper,nsubsteps),npoints,status,initial)
    if (status /= cf_ok) error stop "Perzyna clone setup"
    call copy%import_committed(s,status,previous_input=previous)
    if (status /= cf_ok) error stop "Perzyna restart import"
    step_input = 0.39_wp*x
    call copy%eval(step_input,y,status=status,dt=dt)
    if (status /= cf_ok .or. .not. all(ieee_is_finite(y))) error stop "Perzyna tangent reference response"
    h = maxval(abs(step_input))*epsilon(1.0_wp)**(1.0_wp/3.0_wp)
    call verify_tangent(copy,step_input,error,status,dt=dt, &
        absolute_tolerance=64.0_wp*epsilon(1.0_wp)*maxval(abs(y))/h,step=[(h,i=1,nmandel)])
    if (status /= cf_ok .or. .not. ieee_is_finite(error) .or. error > 1.0_wp) then
        error stop "Substep-exact global endpoint tangent"
    end if

    call copy%eval(step_input,y,status=status,dt=dt)
    if (status /= cf_ok) error stop "Perzyna trial after tangent check"
    call copy%revert()
    call copy%export_committed(s_before,status,previous_input=previous)
    if (status /= cf_ok .or. .not. all(ieee_is_finite(s_before)) .or. maxval(abs(s_before-s)) > 0.0_wp) then
        error stop "Perzyna revert changed committed history"
    end if

    call j2%init(j2_plasticity_t(bulk,shear,sy,hardening),npoints,status)
    if (status /= cf_ok) error stop "J2 rate-independent reference setup"
    call copy%init(perzyna_t(bulk,shear,sy,hardening,1.0e-8_wp),npoints,status,initial)
    if (status /= cf_ok) error stop "Small viscosity Perzyna setup"
    call copy%eval(x,y,status=status,dt=dt)
    if (status /= cf_ok) error stop "Small viscosity Perzyna trial"
    call j2%eval(x,y_steps,status=status)
    if (status /= cf_ok) error stop "Rate-independent J2 trial"
    ! The viscosity error is O(eta/[dt*(3G+H)]), with a trial-stress roundoff floor.
    stress_scale = 3.0_wp*shear*maxval(sqrt(sum(x*x,dim=1)))
    limit_tolerance = (1.0e-8_wp/(dt*(3.0_wp*shear+hardening)) + &
        64.0_wp*epsilon(1.0_wp))*stress_scale
    if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(y_steps)) .or. &
        maxval(abs(y-y_steps)) > limit_tolerance) then
        error stop "Perzyna did not approach J2 limit"
    end if

    call copy%init(perzyna_t(bulk,shear,sy,hardening,eta,1.0e-8_wp,1),npoints,status,initial)
    if (status /= cf_ok) error stop "Substep overflow setup"
    call copy%eval(x,y,status=status,dt=dt)
    if (status /= cf_evaluation_failed) error stop "Substep maximum was not enforced"
    call copy%commit(status)
    if (status /= cf_no_trial) error stop "Failed substep trial committed"

    call copy%init(perzyna_t(bulk,shear,sy,hardening,eta,upper,0),npoints,status,initial)
    if (status /= cf_invalid_argument) error stop "Zero maximum substeps accepted"
    call copy%init(perzyna_t(bulk,shear,sy,hardening,0.0_wp),npoints,status,initial)
    if (status /= cf_invalid_argument) error stop "Zero viscosity accepted"
    call copy%init(perzyna_t(bulk,shear,sy,hardening,eta),npoints,status,initial)
    if (status /= cf_ok) error stop "Perzyna invalid input test setup"
    x(2,2) = ieee_value(0.0_wp,ieee_quiet_nan)
    call copy%eval(x,y,status=status,dt=dt)
    if (status /= cf_evaluation_failed) error stop "Nonfinite Perzyna input accepted"

    call test_exact_relaxation(1.0_wp,1.0_wp)
    call test_exact_relaxation(epsilon(1.0_wp)**2,1.0_wp)
    call test_exact_relaxation(1.0_wp,huge(1.0_wp)/64.0_wp)
    call test_exact_substeps()
    call test_large_moduli()
    call test_branch_switches()
    call test_path_sensitivities()

contains

    subroutine test_exact_relaxation(modulus_scale,time_scale)
        real(wp), intent(in) :: modulus_scale,time_scale
        type(material_t) :: model
        real(wp), parameter :: direction(nmandel) = [1.0_wp,-0.5_wp,-0.5_wp,0.0_wp,0.0_wp,0.0_wp]
        real(wp), parameter :: expected_alpha(2) = [7.0_wp/390.0_wp,49.0_wp/2340.0_wp]
        real(wp), parameter :: expected_stress(2) = [172.0_wp/65.0_wp,446.0_wp/195.0_wp]
        real(wp), parameter :: expected_yield_derivative(2) = [-1.0_wp/234.0_wp,-7.0_wp/1404.0_wp]
        real(wp), parameter :: expected_viscosity_derivative(2) = [-7.0_wp/91260.0_wp,-7.0_wp/273780.0_wp]
        real(wp), parameter :: expected_shear_tangent(2) = [860.0_wp/13.0_wp,4460.0_wp/43.0_wp]
        real(wp), parameter :: expected_overstress(2) = [7.0_wp/10.0_wp,7.0_wp/60.0_wp]
        real(wp) :: initial_input(nmandel,1), input(nmandel,1), output(nmandel,1)
        real(wp) :: jacobian(nmandel,nmandel,1), derivative(nmandel,2,1)
        real(wp) :: state(nstate,1), state_derivative(nstate,2,1), saved_input(nmandel,1)
        real(wp) :: derivative_roundoff, computed_overstress
        integer :: result, load_step, component

        ! For eta/(3G+H)=1/5 and dt=1, backward Euler reduces overstress by 1/6 per strain hold.
        ! The extreme time-unit case has finite eta and dt but eta+dt*(3G+H) would overflow.
        initial_input = 0.0_wp
        call model%init(perzyna_t(100.0_wp*modulus_scale,60.0_wp*modulus_scale, &
            3.0_wp*modulus_scale,15.0_wp*modulus_scale,39.0_wp*modulus_scale*time_scale),1,result, &
            initial_input,sensitivity_parameters=[character(len=12) :: "yield_stress", "viscosity"])
        call expect(result)
        input(:,1) = 0.04_wp*direction
        do load_step = 1,2
            call model%eval(input,output,jacobian,result,dt=time_scale,dsdp=derivative)
            call expect(result)
            ! The eta derivatives can become subnormal when time units are extremely large.
            derivative_roundoff = 8.0_wp*nearest(0.0_wp,1.0_wp)*time_scale
            do component = 1,nmandel
                call close(output(component,1)/modulus_scale,expected_stress(load_step)*direction(component))
                call close(derivative(component,1,1), &
                    -120.0_wp*expected_yield_derivative(load_step)*direction(component))
                call close(derivative(component,2,1)*time_scale, &
                    -120.0_wp*expected_viscosity_derivative(load_step)*direction(component),derivative_roundoff)
            end do
            call close(jacobian(1,1,1)/modulus_scale,1540.0_wp/13.0_wp)
            call close(jacobian(1,2,1)/modulus_scale,1180.0_wp/13.0_wp)
            call close(jacobian(4,4,1)/modulus_scale,expected_shear_tangent(load_step))
            call model%commit(result)
            call expect(result)
            call model%export_committed(state,result,previous_input=saved_input,sensitivity_state=state_derivative)
            call expect(result)
            do component = 1,nmandel
                call close(state(component,1),expected_alpha(load_step)*direction(component))
                call close(state_derivative(component,1,1)*modulus_scale, &
                    expected_yield_derivative(load_step)*direction(component))
                call close(state_derivative(component,2,1)*modulus_scale*time_scale, &
                    expected_viscosity_derivative(load_step)*direction(component),modulus_scale*derivative_roundoff)
                call close(saved_input(component,1),input(component,1))
            end do
            call close(state(nstate,1),expected_alpha(load_step))
            call close(state_derivative(nstate,1,1)*modulus_scale,expected_yield_derivative(load_step))
            call close(state_derivative(nstate,2,1)*modulus_scale*time_scale/expected_viscosity_derivative(load_step), &
                1.0_wp,modulus_scale*derivative_roundoff/abs(expected_viscosity_derivative(load_step)))
            computed_overstress = 1.5_wp*output(1,1)/modulus_scale-3.0_wp-15.0_wp*state(nstate,1)
            call close(computed_overstress,expected_overstress(load_step))
        end do
    end subroutine test_exact_relaxation

    subroutine test_exact_substeps()
        type(material_t) :: model
        real(wp), parameter :: direction(nmandel) = [1.0_wp,-0.5_wp,-0.5_wp,0.0_wp,0.0_wp,0.0_wp]
        real(wp) :: initial_input(nmandel,1), input(nmandel,1), output(nmandel,1)
        real(wp) :: jacobian(nmandel,nmandel,1), derivative(nmandel,2,1)
        real(wp) :: state(nstate,1), state_derivative(nstate,2,1)
        integer :: result, component

        ! Two equal ramp substeps have coefficient 1/273 and history factor 2/7.
        ! Differentiating their closed recurrences gives the global endpoint tangent.
        initial_input = 0.0_wp
        call model%init(perzyna_t(100.0_wp,60.0_wp,3.0_wp,15.0_wp,39.0_wp,0.03_wp,2),1,result, &
            initial_input,sensitivity_parameters=[character(len=12) :: "yield_stress", "viscosity"])
        call expect(result)
        input(:,1) = 0.04_wp*direction
        call model%eval(input,output,jacobian,result,dt=1.0_wp,dsdp=derivative)
        call expect(result)
        do component = 1,nmandel
            call close(output(component,1),(9168.0_wp/3185.0_wp)*direction(component))
            call close(derivative(component,1,1),(360.0_wp/637.0_wp)*direction(component))
            call close(derivative(component,2,1),(736.0_wp/57967.0_wp)*direction(component))
        end do
        call close(jacobian(1,1,1),76260.0_wp/637.0_wp)
        call close(jacobian(1,2,1),57420.0_wp/637.0_wp)
        call close(jacobian(4,4,1),45840.0_wp/637.0_wp)
        call model%commit(result)
        call expect(result)
        call model%export_committed(state,result,previous_input=initial_input,sensitivity_state=state_derivative)
        call expect(result)
        do component = 1,nmandel
            call close(state(component,1),(51.0_wp/3185.0_wp)*direction(component))
            call close(state_derivative(component,1,1),-(3.0_wp/637.0_wp)*direction(component))
            call close(state_derivative(component,2,1),-(92.0_wp/869505.0_wp)*direction(component))
        end do
        call close(state(nstate,1),51.0_wp/3185.0_wp)
        call close(state_derivative(nstate,1,1),-3.0_wp/637.0_wp)
        call close(state_derivative(nstate,2,1),-92.0_wp/869505.0_wp)
    end subroutine test_exact_substeps


    subroutine test_large_moduli()
        type(material_t) :: model
        real(wp), parameter :: large_scale = huge(1.0_wp)
        real(wp) :: initial_input(nmandel,1), output(nmandel,1), jacobian(nmandel,nmandel,1)
        real(wp) :: state(nstate,1), saved_input(nmandel,1), expected(nmandel,nmandel)
        integer :: result, component, column

        ! All elastic entries are finite although the intermediate product 4*G would overflow.
        initial_input = 0.0_wp
        call model%init(perzyna_t(large_scale/32.0_wp,large_scale/3.5_wp,1.0_wp,0.0_wp,1.0_wp),1,result,initial_input)
        call expect(result)
        call model%eval(initial_input,output,jacobian,result,dt=1.0_wp)
        call expect(result)
        expected = 0.0_wp
        expected(1:3,1:3) = 1.0_wp/32.0_wp-4.0_wp/21.0_wp
        do component = 1,3
            expected(component,component) = 1.0_wp/32.0_wp+8.0_wp/21.0_wp
        end do
        do component = 4,nmandel
            expected(component,component) = 4.0_wp/7.0_wp
        end do
        do component = 1,nmandel
            call close(output(component,1),0.0_wp)
            do column = 1,nmandel
                call close(jacobian(component,column,1)/large_scale,expected(component,column))
            end do
        end do
        call model%commit(result)
        call expect(result)
        call model%export_committed(state,result,previous_input=saved_input)
        call expect(result)
        do component = 1,nstate
            call close(state(component,1),0.0_wp)
        end do
    end subroutine test_large_moduli

    subroutine test_branch_switches()
        type(material_t) :: ordinary, tracked
        real(wp), parameter :: direction(nmandel) = [1.0_wp,-0.5_wp,-0.5_wp,0.0_wp,0.0_wp,0.0_wp]
        real(wp) :: initial_input(nmandel,1), input(nmandel,1), output(nmandel,1)
        real(wp) :: jacobian(nmandel,nmandel,1), derivative(nmandel,1,1)
        real(wp) :: expected_stress
        integer :: result, branch, component

        initial_input = 0.0_wp
        call ordinary%init(perzyna_t(128.0_wp,64.0_wp,3.0_wp,0.0_wp,48.0_wp,3.0_wp/128.0_wp,2), &
            1,result,initial_input)
        call expect(result)
        call tracked%init(perzyna_t(128.0_wp,64.0_wp,3.0_wp,0.0_wp,48.0_wp,3.0_wp/128.0_wp,2), &
            1,result,initial_input,sensitivity_parameters=[character(len=12) :: "yield_stress"])
        call expect(result)
        do branch = 1,2
            ! One-step endpoint yield, then yield at the first of two internal substeps.
            input(:,1) = real(branch,wp)*direction/64.0_wp
            call ordinary%eval(input,output,jacobian,result,dt=1.0_wp)
            call expect(result)
            expected_stress = 2.0_wp
            if (branch==2) expected_stress = 8.0_wp/3.0_wp
            do component = 1,nmandel
                call close(output(component,1),expected_stress*direction(component))
            end do
            if (branch==1) then
                call close(jacobian(1,1,1),640.0_wp/3.0_wp)
                call close(jacobian(1,2,1),256.0_wp/3.0_wp)
                call close(jacobian(4,4,1),128.0_wp)
            else
                call close(jacobian(1,1,1),1408.0_wp/9.0_wp)
                call close(jacobian(1,2,1),1024.0_wp/9.0_wp)
                call close(jacobian(4,4,1),256.0_wp/3.0_wp)
            end if
            call tracked%eval(input,output,status=result,dt=1.0_wp,dsdp=derivative)
            if (result /= cf_nondifferentiable) error stop "Perzyna sensitivity accepted an internal yield switch"
        end do
    end subroutine test_branch_switches

    subroutine test_path_sensitivities()
        integer, parameter :: nparameters = 5, nloads = 4
        integer, parameter :: subdivisions(nloads) = [3,1,2,4]
        character(len=18), parameter :: names(nparameters) = [character(len=18) :: "bulk_modulus", &
            "shear_modulus","yield_stress","hardening_modulus","viscosity"]
        real(wp), parameter :: parameters(nparameters) = [100.0_wp,60.0_wp,3.0_wp,15.0_wp,39.0_wp]
        real(wp), parameter :: elapsed(nloads) = [0.8_wp,1.25_wp,0.6_wp,1.0_wp]
        type(material_t) :: tracked
        type(material_t) :: plus(nparameters,2), minus(nparameters,2)
        real(wp) :: path(nmandel,nloads), initial_input(nmandel,1), input(nmandel,1), output(nmandel,1)
        real(wp) :: derivative(nmandel,nparameters,1), state(nstate,1), state_derivative(nstate,nparameters,1)
        real(wp) :: saved_input(nmandel,1), perturbed_input(nmandel,1)
        real(wp) :: pplus(nparameters), pminus(nparameters), perturbations(nparameters,2), elastic_dev(nmandel)
        real(wp) :: yp(nmandel,1), ym(nmandel,1), hp(nstate,1), hm(nstate,1)
        real(wp) :: stress_difference(nmandel,2), state_difference(nstate,2)
        real(wp) :: stress_tolerance(nmandel), state_tolerance(nstate)
        real(wp) :: stress_scale, strain_scale, accumulated_stress_scale, accumulated_strain_scale, trace_e
        integer :: result, load_step, parameter, refinement

        ! A subdivided load, strain hold, elastic unloading, and nonproportional reversal.
        ! The prescribed increments select 3,1,2,4 substeps independently of physical parameters.
        path(:,1) = [0.044_wp,-0.018_wp,-0.014_wp,0.012_wp,0.004_wp,0.0_wp]
        path(:,2) = path(:,1)
        path(:,3) = 0.55_wp*path(:,1)
        path(:,4) = [-0.031_wp,0.014_wp,0.011_wp,-0.015_wp,0.008_wp,0.013_wp]
        initial_input = 0.0_wp
        call tracked%init(perzyna_t(parameters(1),parameters(2),parameters(3),parameters(4),parameters(5), &
            0.015_wp,maxval(subdivisions)),1,result,initial_input,sensitivity_parameters=names)
        call expect(result)
        do parameter = 1,nparameters
            perturbations(parameter,1) = (32.0_wp*epsilon(1.0_wp))**(1.0_wp/3.0_wp)*parameters(parameter)
            ! Bulk stress is affine in K and does not affect effective J2 history.
            if (parameter==1) perturbations(parameter,1) = parameters(parameter)/4.0_wp
            perturbations(parameter,2) = 0.5_wp*perturbations(parameter,1)
            do refinement = 1,2
                pplus = parameters
                pminus = parameters
                pplus(parameter) = pplus(parameter)+perturbations(parameter,refinement)
                pminus(parameter) = pminus(parameter)-perturbations(parameter,refinement)
                call plus(parameter,refinement)%init(perzyna_t(pplus(1),pplus(2),pplus(3),pplus(4),pplus(5), &
                    0.015_wp,maxval(subdivisions)),1,result,initial_input)
                call expect(result)
                call minus(parameter,refinement)%init(perzyna_t(pminus(1),pminus(2),pminus(3),pminus(4),pminus(5), &
                    0.015_wp,maxval(subdivisions)),1,result,initial_input)
                call expect(result)
            end do
        end do
        state = 0.0_wp
        accumulated_stress_scale = 0.0_wp
        accumulated_strain_scale = 0.0_wp
        do load_step = 1,nloads
            input(:,1) = path(:,load_step)
            trace_e = sum(input(1:3,1))
            elastic_dev = input(:,1)-state(1:nmandel,1)
            elastic_dev(1:3) = elastic_dev(1:3)-trace_e/3.0_wp
            stress_scale = max(abs(parameters(1)*trace_e),parameters(3)+parameters(4)*state(nstate,1), &
                2.0_wp*parameters(2)*sqrt(1.5_wp*dot_product(elastic_dev,elastic_dev)))
            call tracked%eval(input,output,status=result,dt=elapsed(load_step),dsdp=derivative)
            call expect(result)
            call tracked%commit(result)
            call expect(result)
            call tracked%export_committed(state,result,previous_input=saved_input,sensitivity_state=state_derivative)
            call expect(result)
            if (.not. all(ieee_is_finite(output)) .or. .not. all(ieee_is_finite(state)) .or. &
                .not. all(ieee_is_finite(derivative)) .or. .not. all(ieee_is_finite(state_derivative)) .or. &
                .not. all(ieee_is_finite(saved_input))) error stop "Nonfinite Perzyna history comparison"
            if (maxval(abs(saved_input-input)) > 0.0_wp) error stop "Perzyna accepted input history mismatch"
            strain_scale = max(maxval(abs(input)),maxval(abs(state)))
            accumulated_stress_scale = accumulated_stress_scale+real(subdivisions(load_step),wp)*stress_scale
            accumulated_strain_scale = accumulated_strain_scale+real(subdivisions(load_step),wp)*strain_scale
            if (.not. all(ieee_is_finite([accumulated_stress_scale,accumulated_strain_scale]))) then
                error stop "Nonfinite Perzyna history error scale"
            end if
            do parameter = 1,nparameters
                do refinement = 1,2
                    call plus(parameter,refinement)%eval(input,yp,status=result,dt=elapsed(load_step))
                    call expect(result)
                    call minus(parameter,refinement)%eval(input,ym,status=result,dt=elapsed(load_step))
                    call expect(result)
                    call plus(parameter,refinement)%commit(result)
                    call expect(result)
                    call minus(parameter,refinement)%commit(result)
                    call expect(result)
                    call plus(parameter,refinement)%export_committed(hp,result,previous_input=perturbed_input)
                    call expect(result)
                    call minus(parameter,refinement)%export_committed(hm,result,previous_input=perturbed_input)
                    call expect(result)
                    stress_difference(:,refinement) = (yp(:,1)-ym(:,1))/(2.0_wp*perturbations(parameter,refinement))
                    state_difference(:,refinement) = (hp(:,1)-hm(:,1))/(2.0_wp*perturbations(parameter,refinement))
                end do
                if (.not. all(ieee_is_finite(stress_difference)) .or. .not. all(ieee_is_finite(state_difference))) then
                    error stop "Nonfinite Perzyna history secants"
                end if
                ! The roundoff allowance accounts for all accepted backward-Euler substeps.
                stress_tolerance = abs(stress_difference(:,2)-stress_difference(:,1))/3.0_wp &
                    + 32.0_wp*epsilon(1.0_wp)*accumulated_stress_scale/perturbations(parameter,2)
                state_tolerance = abs(state_difference(:,2)-state_difference(:,1))/3.0_wp &
                    + 32.0_wp*epsilon(1.0_wp)*accumulated_strain_scale/perturbations(parameter,2)
                if (any(abs(stress_difference(:,2)-derivative(:,parameter,1)) > stress_tolerance)) then
                    print *, "Perzyna response history derivative",load_step,parameter
                    error stop "Perzyna parameter derivative mismatch"
                end if
                if (any(abs(state_difference(:,2)-state_derivative(:,parameter,1)) > state_tolerance)) then
                    print *, "Perzyna internal history derivative",load_step,parameter
                    error stop "Perzyna state derivative mismatch"
                end if
            end do
        end do
    end subroutine test_path_sensitivities

    subroutine expect(result)
        integer, intent(in) :: result

        if (result /= cf_ok) then
            print *, "Perzyna status",result
            error stop "Exact Perzyna reference evaluation failed"
        end if
    end subroutine expect

    subroutine close(actual, expected, subnormal_roundoff)
        real(wp), intent(in) :: actual, expected
        real(wp), intent(in), optional :: subnormal_roundoff
        real(wp) :: bound

        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) then
            error stop "Nonfinite exact Perzyna comparison"
        end if
        bound = 512.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(expected))
        if (present(subnormal_roundoff)) bound = bound+subnormal_roundoff
        if (abs(actual-expected) > bound) then
            print *, "Expected",expected,"got",actual
            error stop "Exact Perzyna reference mismatch"
        end if
    end subroutine close
end program test_perzyna

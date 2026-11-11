! SPDX-License-Identifier: Apache-2.0

program test_maxwell_shear
    !! Test convention: moduli, stress and overstress in Pa. Relaxation times
    !! and dt are in s. Explicit rescaling cases state their factors.
    !!
    !! Independent Maxwell relaxation, tensor tangent and accepted-path tests.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use constituf, only: wp, lk, material_t, cf_ok, cf_invalid_argument, &
                        cf_invalid_shape, cf_invalid_time, cf_evaluation_failed, cf_no_trial
    use constituf_maxwell_shear, only: maxwell_shear_t
    implicit none (type, external)

    integer, parameter :: nmandel = 6, nparameters = 4, npoints = 2
    real(wp), parameter :: bulk = 100.0_wp, equilibrium = 10.0_wp
    real(wp), parameter :: branch = 30.0_wp, relaxation = 2.0_wp

    call test_spec_and_limits()
    call test_relaxation_and_tangent()
    call test_short_relaxation_reference(0.25_wp*sqrt(tiny(1.0_wp)))
    call test_short_relaxation_reference(1.0_wp/huge(1.0_wp),huge(1.0_wp)/8.0_wp)
    call test_subnormal_attenuation_reference()
    call test_hydrostatic_history()
    call test_path_sensitivities()
    call test_unit_conventions()
    call test_invalid_trials()

contains

    subroutine require_status(actual, expected)
        integer, intent(in) :: actual, expected

        if (actual /= expected) then
            print *, "Expected status", expected, "received", actual
            error stop "Maxwell material status mismatch"
        end if
    end subroutine require_status

    subroutine close(actual, expected)
        real(wp), intent(in) :: actual, expected

        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) then
            error stop "Nonfinite Maxwell comparison"
        end if
        if (abs(actual-expected) > 512.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(expected))) then
            print *, "Expected", expected, "received", actual
            error stop "Maxwell reference mismatch"
        end if
    end subroutine close

    subroutine test_unit_conventions()
        !! One physical ramp/hold/unload history in Pa,s and MPa,ms.
        !! Stress, tangent and overstress scale by the stress factor.
        !! Their relaxation-time derivatives scale by stress/time.
        real(wp), parameter :: stress_factors(2) = [1.0_wp,1.0e-6_wp]
        real(wp), parameter :: time_factors(2) = [1.0_wp,1000.0_wp]
        real(wp), parameter :: amplitudes(3) = [1.0_wp,1.0_wp,-0.5_wp]
        type(material_t) :: materials(2)
        real(wp) :: initial(nmandel,1), x(nmandel,1)
        real(wp) :: stress(nmandel,1,2), tangent(nmandel,nmandel,1,2)
        real(wp) :: derivative(nmandel,1,1,2), history(nmandel,1,2)
        real(wp) :: history_derivative(nmandel,1,1,2), previous(nmandel,1,2)
        integer :: status, system, load_step, row, column

        initial = 0.0_wp
        do system = 1,2
            call materials(system)%init(maxwell_shear_t(100.0e9_wp*stress_factors(system), &
                10.0e9_wp*stress_factors(system),30.0e9_wp*stress_factors(system), &
                2.0_wp*time_factors(system)),1,status,initial, &
                sensitivity_parameters=[character(len=15) :: "relaxation_time"])
            call require_status(status,cf_ok)
        end do
        do load_step = 1,size(amplitudes)
            x(:,1) = amplitudes(load_step)*[0.010_wp,-0.005_wp,0.0_wp,0.002_wp,0.0_wp,0.003_wp]
            do system = 1,2
                call materials(system)%eval(x,stress(:,:,system),tangent(:,:,:,system),status, &
                    dt=0.5_wp*time_factors(system),dsdp=derivative(:,:,:,system))
                call require_status(status,cf_ok)
                call materials(system)%commit(status)
                call require_status(status,cf_ok)
                call materials(system)%export_committed(history(:,:,system),status, &
                    previous_input=previous(:,:,system),sensitivity_state=history_derivative(:,:,:,system))
                call require_status(status,cf_ok)
            end do
            do row = 1,nmandel
                call close(stress(row,1,2)/stress_factors(2),stress(row,1,1))
                call close(history(row,1,2)/stress_factors(2),history(row,1,1))
                call close(derivative(row,1,1,2)*time_factors(2)/stress_factors(2),derivative(row,1,1,1))
                call close(history_derivative(row,1,1,2)*time_factors(2)/stress_factors(2), &
                    history_derivative(row,1,1,1))
                call close(previous(row,1,1),x(row,1))
                call close(previous(row,1,2),x(row,1))
                do column = 1,nmandel
                    call close(tangent(row,column,1,2)/stress_factors(2),tangent(row,column,1,1))
                end do
            end do
        end do
    end subroutine test_unit_conventions

    subroutine test_spec_and_limits()
        type(material_t) :: material
        real(wp) :: initial(nmandel,npoints), x(nmandel,npoints), y(nmandel,npoints)
        integer(lk) :: history_bytes
        integer :: status

        initial = 0.0_wp
        call material%init(maxwell_shear_t(bulk,equilibrium,branch,relaxation), &
                           npoints, status)
        call require_status(status, cf_invalid_argument)
        call material%init(maxwell_shear_t(bulk,equilibrium,branch,relaxation), &
                           npoints, status, initial_input=initial)
        call require_status(status, cf_ok)
        if (material%input_size() /= nmandel .or. material%output_size() /= nmandel) then
            error stop "Maxwell tensor dimensions"
        end if
        if (material%state_size() /= nmandel) error stop "Maxwell history dimension"
        if (material%input_offset("strain.small") /= 1) error stop "Strain offset"
        if (material%output_offset("stress.cauchy") /= 1) error stop "Stress offset"
        if (.not. material%field_matches("state", "overstress.deviatoric", &
            "mandel", "internal_stress", "spatial", [3,3])) then
            error stop "Maxwell history field convention"
        end if
        history_bytes = 2_lk*int(npoints,lk)*int(nmandel,lk)* &
                        int(storage_size(0.0_wp)/8,lk)
        if (material%state_storage_bytes() /= history_bytes) then
            error stop "Unexpected committed/trial history memory"
        end if

        call material%init(maxwell_shear_t(0.0_wp,equilibrium,branch,relaxation), &
                           npoints, status, initial)
        call require_status(status, cf_invalid_argument)
        call material%init(maxwell_shear_t(bulk,-1.0_wp,branch,relaxation), &
                           npoints, status, initial)
        call require_status(status, cf_invalid_argument)
        call material%init(maxwell_shear_t(bulk,0.0_wp,0.0_wp,relaxation), &
                           npoints, status, initial)
        call require_status(status, cf_invalid_argument)
        call material%init(maxwell_shear_t(bulk,equilibrium,branch,0.0_wp), &
                           npoints, status, initial)
        call require_status(status, cf_invalid_argument)

        ! G1=0 is the elastic limiting case; dt remains mandatory by schema.
        call material%init(maxwell_shear_t(bulk,equilibrium,0.0_wp,relaxation), &
                           npoints, status, initial)
        call require_status(status, cf_ok)
        x = 0.0_wp
        x(1,:) = 0.02_wp
        call material%eval(x, y, status=status, dt=1.0_wp)
        call require_status(status, cf_ok)
        call close(y(1,1), bulk*0.02_wp + &
                   2.0_wp*equilibrium*(0.02_wp-0.02_wp/3.0_wp))
        call material%commit(status)
        call require_status(status, cf_ok)
        call material%eval(x, y, status=status, dt=1.0_wp)
        call require_status(status, cf_ok)
        call close(y(1,1), bulk*0.02_wp + &
                   2.0_wp*equilibrium*(0.02_wp-0.02_wp/3.0_wp))
    end subroutine test_spec_and_limits

    subroutine test_relaxation_and_tangent()
        type(material_t) :: material
        real(wp), parameter :: overstress_factor(2) = [40.0_wp,80.0_wp/3.0_wp]
        real(wp), parameter :: branch_derivative(2) = [4.0_wp/3.0_wp,8.0_wp/9.0_wp]
        real(wp), parameter :: time_derivative(2) = [20.0_wp/3.0_wp,80.0_wp/9.0_wp]
        real(wp) :: initial(nmandel,npoints), x(nmandel,npoints), y(nmandel,npoints)
        real(wp) :: tangent(nmandel,nmandel,npoints), positive(nmandel,npoints)
        real(wp) :: negative(nmandel,npoints), perturb(nmandel,npoints)
        real(wp) :: dsdp(nmandel,nparameters,npoints), state(nmandel,npoints)
        real(wp) :: state_derivative(nmandel,nparameters,npoints), previous(nmandel,npoints)
        real(wp) :: h, derivative, roundoff
        integer :: i, j, status, step

        initial = 0.0_wp
        call material%init(maxwell_shear_t(bulk,equilibrium,branch,relaxation), &
            npoints,status,initial,sensitivity_parameters=[character(len=26) :: &
                "bulk_modulus","equilibrium_shear_modulus", &
                "branch_shear_modulus","relaxation_time"])
        call require_status(status,cf_ok)
        x(:,1) = [0.02_wp,-0.01_wp,-0.01_wp,0.008_wp,-0.006_wp,0.005_wp]
        x(:,2) = [0.01_wp,0.01_wp,0.01_wp,0.0_wp,0.0_wp,0.0_wp]
        ! For dt=1 and tau=2, each accepted hold multiplies q by 2/3.
        ! Differentiating the ramp and hold gives 60*a' and 120*a*a',
        ! with a=2/3 and a'=1/9. These are total history derivatives.
        do step = 1,2
            call material%eval(x,y,tangent,status,dt=1.0_wp,dsdp=dsdp)
            call require_status(status,cf_ok)
            if (step == 2) then
                call material%revert()
                call material%eval(x,y,tangent,status,dt=1.0_wp,dsdp=dsdp)
                call require_status(status,cf_ok)
            end if
            if (.not. all(ieee_is_finite(dsdp))) error stop "Nonfinite Maxwell relaxation derivatives"
            do i = 1,nmandel
                call close(y(i,1),(20.0_wp+overstress_factor(step))*x(i,1))
                call close(dsdp(i,1,1),0.0_wp)
                call close(dsdp(i,2,1),2.0_wp*x(i,1))
                call close(dsdp(i,3,1),branch_derivative(step)*x(i,1))
                call close(dsdp(i,4,1),time_derivative(step)*x(i,1))
                call close(y(i,2),300.0_wp*x(i,2))
                call close(dsdp(i,1,2),3.0_wp*x(i,2))
                call close(maxval(abs(dsdp(i,2:nparameters,2))),0.0_wp)
            end do
            call close(tangent(1,1,1),140.0_wp)
            call close(tangent(1,2,1),80.0_wp)
            call close(tangent(4,4,1),60.0_wp)
            call material%commit(status)
            call require_status(status,cf_ok)
            call material%export_committed(state,status,previous_input=previous, &
                sensitivity_state=state_derivative)
            call require_status(status,cf_ok)
            if (.not. all(ieee_is_finite(previous)) .or. .not. all(ieee_is_finite(state_derivative))) then
                error stop "Nonfinite Maxwell relaxation history"
            end if
            call close(maxval(abs(previous-x)),0.0_wp)
            do i = 1,nmandel
                call close(state(i,1),overstress_factor(step)*x(i,1))
                call close(state_derivative(i,1,1),0.0_wp)
                call close(state_derivative(i,2,1),0.0_wp)
                call close(state_derivative(i,3,1),branch_derivative(step)*x(i,1))
                call close(state_derivative(i,4,1),time_derivative(step)*x(i,1))
                call close(state(i,2),0.0_wp)
                call close(maxval(abs(state_derivative(i,:,2))),0.0_wp)
            end do
        end do

        ! The discrete tangent is checked at a nonzero committed state.
        ! Stress is affine in strain, so a finite step has no truncation error.
        x(1,1) = x(1,1)+0.005_wp
        x(4,1) = x(4,1)-0.004_wp
        h = 0.25_wp*maxval(abs(x))
        call material%eval(x,y,tangent,status,dt=0.35_wp,dsdp=dsdp)
        call require_status(status,cf_ok)
        if (.not. all(ieee_is_finite(tangent))) error stop "Nonfinite Maxwell tangent"
        do j = 1,nmandel
            perturb = x
            perturb(j,1) = perturb(j,1)+h
            call material%eval(perturb,positive,status=status,dt=0.35_wp,dsdp=dsdp)
            call require_status(status,cf_ok)
            perturb(j,1) = perturb(j,1)-2.0_wp*h
            call material%eval(perturb,negative,status=status,dt=0.35_wp,dsdp=dsdp)
            call require_status(status,cf_ok)
            if (.not. all(ieee_is_finite(positive)) .or. .not. all(ieee_is_finite(negative))) then
                error stop "Nonfinite Maxwell tangent secant response"
            end if
            roundoff = 32.0_wp*epsilon(1.0_wp)*maxval(abs(positive)+abs(negative))/(2.0_wp*h)
            if (.not. ieee_is_finite(roundoff)) error stop "Nonfinite Maxwell tangent error bound"
            do i = 1,nmandel
                derivative = (positive(i,1)-negative(i,1))/(2.0_wp*h)
                if (.not. ieee_is_finite(derivative) .or. abs(derivative-tangent(i,j,1)) > roundoff) then
                    error stop "Inconsistent backward-Euler Maxwell tangent"
                end if
            end do
        end do
    end subroutine test_relaxation_and_tangent

    subroutine test_short_relaxation_reference(time_constant,large_time_step)
        !! Pure Maxwell shear retains a small, representable overstress even
        !! when dt/tau or dt/tau**2 would overflow. Normalization exposes it.
        real(wp), intent(in) :: time_constant
        real(wp), intent(in), optional :: large_time_step
        type(material_t) :: material
        real(wp) :: initial(nmandel,1), x(nmandel,1), y(nmandel,1)
        real(wp) :: tangent(nmandel,nmandel,1), dsdp(nmandel,nparameters,1)
        real(wp) :: state(nmandel,1), state_derivative(nmandel,nparameters,1)
        real(wp) :: previous(nmandel,1)
        real(wp) :: subnormal_bound
        integer :: status

        initial = 0.0_wp
        call material%init(maxwell_shear_t(bulk,0.0_wp,branch,time_constant), &
            1,status,initial,sensitivity_parameters=[character(len=26) :: &
                "bulk_modulus","equilibrium_shear_modulus", &
                "branch_shear_modulus","relaxation_time"])
        call require_status(status,cf_ok)
        x = 0.0_wp
        x(nmandel,1) = 0.125_wp
        call material%eval(x,y,tangent,status,dt=2.0_wp,dsdp=dsdp)
        call require_status(status,cf_ok)
        call close(y(nmandel,1)/time_constant,15.0_wp/4.0_wp)
        call close(tangent(nmandel,nmandel,1)/time_constant,30.0_wp)
        call close(dsdp(nmandel,1,1),0.0_wp)
        call close(dsdp(nmandel,2,1),0.25_wp)
        call close(dsdp(nmandel,3,1)/time_constant,0.125_wp)
        call close(dsdp(nmandel,4,1),15.0_wp/4.0_wp)
        call material%commit(status)
        call require_status(status,cf_ok)
        call material%export_committed(state,status,previous_input=previous, &
            sensitivity_state=state_derivative)
        call require_status(status,cf_ok)
        call close(state(nmandel,1)/time_constant,15.0_wp/4.0_wp)
        call close(state_derivative(nmandel,3,1)/time_constant,0.125_wp)
        call close(state_derivative(nmandel,4,1),15.0_wp/4.0_wp)

        if (.not. present(large_time_step)) return
        ! G1=dt=L and tau=t: tau/dt rounds to zero, but stress=t/4 and C66=2*t are representable.
        ! The normalized tolerance includes the subnormal quantum instead of treating zero as the limit.
        subnormal_bound = 8.0_wp*nearest(0.0_wp,1.0_wp)/time_constant
        call material%init(maxwell_shear_t(bulk,0.0_wp,large_time_step,time_constant), &
            1,status,initial,sensitivity_parameters=[character(len=26) :: &
                "bulk_modulus","equilibrium_shear_modulus","branch_shear_modulus","relaxation_time"])
        call require_status(status,cf_ok)
        call material%eval(x,y,tangent,status,dt=large_time_step,dsdp=dsdp)
        call require_status(status,cf_ok)
        if (.not. all(ieee_is_finite([y(nmandel,1),tangent(nmandel,nmandel,1),subnormal_bound]))) then
            error stop "Nonfinite coupled-scale Maxwell reference"
        end if
        if (abs(y(nmandel,1)/time_constant-0.25_wp) > subnormal_bound .or. &
            abs(tangent(nmandel,nmandel,1)/time_constant-2.0_wp) > subnormal_bound) then
            error stop "Maxwell discarded representable weighted stress or tangent"
        end if
        call close(dsdp(nmandel,4,1),0.25_wp)
        call material%commit(status)
        call require_status(status,cf_ok)
        call material%export_committed(state,status,previous_input=previous,sensitivity_state=state_derivative)
        call require_status(status,cf_ok)
        if (.not. ieee_is_finite(state(nmandel,1)) .or. &
            abs(state(nmandel,1)/time_constant-0.25_wp) > subnormal_bound) then
            error stop "Maxwell discarded representable weighted history"
        end if
        call close(state_derivative(nmandel,4,1),0.25_wp)

        ! A short load creates q=L*t/8 and q_tau=L/16. Holding for L then gives q_tau/t=3/16:
        ! propagated history contributes 1/16 and the direct derivative contributes 1/8.
        ! Prescribed strain and dt stay fixed when differentiating the parameter tau.
        call material%init(maxwell_shear_t(bulk,0.0_wp,large_time_step,time_constant), &
            1,status,initial,sensitivity_parameters=[character(len=26) :: &
                "bulk_modulus","equilibrium_shear_modulus","branch_shear_modulus","relaxation_time"])
        call require_status(status,cf_ok)
        x(nmandel,1) = time_constant/8.0_wp
        call material%eval(x,y,status=status,dt=time_constant,dsdp=dsdp)
        call require_status(status,cf_ok)
        call close(y(nmandel,1)/(large_time_step*time_constant),0.125_wp)
        call close(dsdp(nmandel,4,1)/large_time_step,1.0_wp/16.0_wp)
        call material%commit(status)
        call require_status(status,cf_ok)
        call material%eval(x,y,status=status,dt=large_time_step,dsdp=dsdp)
        call require_status(status,cf_ok)
        call close(y(nmandel,1),0.0_wp)
        if (.not. ieee_is_finite(dsdp(nmandel,4,1)) .or. &
            abs(dsdp(nmandel,4,1)/time_constant-3.0_wp/16.0_wp) > subnormal_bound) then
            error stop "Maxwell lost propagated short-relaxation sensitivity"
        end if
        call material%commit(status)
        call require_status(status,cf_ok)
        call material%export_committed(state,status,previous_input=previous,sensitivity_state=state_derivative)
        call require_status(status,cf_ok)
        call close(state(nmandel,1),0.0_wp)
        if (.not. ieee_is_finite(state_derivative(nmandel,4,1)) .or. &
            abs(state_derivative(nmandel,4,1)/time_constant-3.0_wp/16.0_wp) > subnormal_bound) then
            error stop "Maxwell lost committed short-relaxation sensitivity"
        end if
    end subroutine test_short_relaxation_reference

    subroutine test_subnormal_attenuation_reference()
        !! A nonzero rounded attenuation must not quantize normal weighted
        !! stress, tangent or the derivative of an accepted branch history.
        type(material_t) :: material
        real(wp) :: initial(nmandel,1), x(nmandel,1), y(nmandel,1)
        real(wp) :: tangent(nmandel,nmandel,1), dsdp(nmandel,1,1)
        real(wp) :: state(nmandel,1), state_derivative(nmandel,1,1)
        real(wp) :: previous(nmandel,1)
        real(wp) :: quantum, branch_modulus, stress_scale, time_constant
        integer :: status, step

        quantum = nearest(0.0_wp,1.0_wp)
        time_constant = 5.0_wp*quantum
        branch_modulus = huge(1.0_wp)/8.0_wp
        stress_scale = branch_modulus*quantum
        initial = 0.0_wp
        x = 0.0_wp
        x(nmandel,1) = 0.125_wp
        call material%init(maxwell_shear_t(bulk,0.0_wp,branch_modulus,time_constant), &
            1,status,initial,sensitivity_parameters=[character(len=15) :: "relaxation_time"])
        call require_status(status,cf_ok)

        ! tau/(tau+4) is 5*quantum/4 to working accuracy, but rounds to
        ! quantum before weighting. All normalized nonzero references are normal.
        do step = 1,2
            call material%eval(x,y,tangent,status,dt=4.0_wp,dsdp=dsdp)
            call require_status(status,cf_ok)
            if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(tangent)) .or. &
                .not. all(ieee_is_finite(dsdp))) error stop "Nonfinite subnormal-attenuation response"
            call close(tangent(nmandel,nmandel,1)/stress_scale,5.0_wp/2.0_wp)
            if (step == 1) then
                call close(y(nmandel,1)/stress_scale,5.0_wp/16.0_wp)
                call close(dsdp(nmandel,1,1)/branch_modulus,1.0_wp/16.0_wp)
            else
                if (abs(y(nmandel,1)) > 0.0_wp) error stop "Unrepresentable relaxed overstress"
                ! Each tau derivative contribution is 5/64 of stress_scale:
                ! propagation of the old derivative and differentiation of attenuation.
                call close(dsdp(nmandel,1,1)/stress_scale,5.0_wp/32.0_wp)
            end if
            call material%commit(status)
            call require_status(status,cf_ok)
            call material%export_committed(state,status,previous_input=previous,sensitivity_state=state_derivative)
            call require_status(status,cf_ok)
            if (.not. all(ieee_is_finite(state)) .or. .not. all(ieee_is_finite(state_derivative)) .or. &
                .not. all(ieee_is_finite(previous))) then
                error stop "Nonfinite subnormal-attenuation history"
            end if
            call close(maxval(abs(previous-x)),0.0_wp)
            if (step == 1) then
                call close(state(nmandel,1)/stress_scale,5.0_wp/16.0_wp)
                call close(state_derivative(nmandel,1,1)/branch_modulus,1.0_wp/16.0_wp)
            else
                if (abs(state(nmandel,1)) > 0.0_wp) error stop "Unrepresentable committed overstress"
                call close(state_derivative(nmandel,1,1)/stress_scale,5.0_wp/32.0_wp)
            end if
        end do

        ! A subnormal attenuation can coexist with an overflowing G1*tau.
        ! The weighted shear modulus tends to 3/2, so C66=3 and stress6=3/8.
        call material%init(maxwell_shear_t(bulk,0.0_wp,huge(1.0_wp)/2.0_wp,3.0_wp), &
            1,status,initial,sensitivity_parameters=[character(len=15) :: "relaxation_time"])
        call require_status(status,cf_ok)
        call material%eval(x,y,tangent,status,dt=huge(1.0_wp),dsdp=dsdp)
        call require_status(status,cf_ok)
        call close(y(nmandel,1),3.0_wp/8.0_wp)
        call close(tangent(nmandel,nmandel,1),3.0_wp)
        call close(dsdp(nmandel,1,1),1.0_wp/8.0_wp)

        ! The larger shear makes the unattenuated forcing 2*huge, while
        ! its attenuated stress and relaxation-time derivative remain finite.
        x(nmandel,1) = 2.0_wp
        call material%eval(x,y,tangent,status,dt=huge(1.0_wp),dsdp=dsdp)
        call require_status(status,cf_ok)
        if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(tangent)) .or. &
            .not. all(ieee_is_finite(dsdp))) error stop "Nonfinite attenuated large-forcing response"
        call close(y(nmandel,1),6.0_wp)
        call close(tangent(nmandel,nmandel,1),3.0_wp)
        call close(dsdp(nmandel,1,1),2.0_wp)

        ! With tau=dt=1 the attenuation is exactly one half: the same
        ! unrepresentable forcing returns the representable endpoint huge.
        call material%init(maxwell_shear_t(bulk,0.0_wp,huge(1.0_wp)/2.0_wp,1.0_wp), &
            1,status,initial,sensitivity_parameters=[character(len=15) :: "relaxation_time"])
        call require_status(status,cf_ok)
        call material%eval(x,y,tangent,status,dt=1.0_wp,dsdp=dsdp)
        call require_status(status,cf_ok)
        if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(tangent)) .or. &
            .not. all(ieee_is_finite(dsdp))) error stop "Nonfinite normal-attenuation large-forcing response"
        call close(y(nmandel,1)/huge(1.0_wp),1.0_wp)
        call close(tangent(nmandel,nmandel,1)/huge(1.0_wp),0.5_wp)
        call close(dsdp(nmandel,1,1)/huge(1.0_wp),0.5_wp)
        x(nmandel,1) = 0.125_wp

        ! The complementary short-step fraction dt/(tau+dt) can also be
        ! subnormal. Its unrounded value controls a normal tau derivative.
        call material%init(maxwell_shear_t(bulk,0.0_wp,branch_modulus,4.0_wp), &
            1,status,initial,sensitivity_parameters=[character(len=15) :: "relaxation_time"])
        call require_status(status,cf_ok)
        call material%eval(x,y,tangent,status,dt=time_constant,dsdp=dsdp)
        call require_status(status,cf_ok)
        call close(y(nmandel,1)/branch_modulus,0.25_wp)
        call close(tangent(nmandel,nmandel,1)/branch_modulus,2.0_wp)
        call close(dsdp(nmandel,1,1)/stress_scale,5.0_wp/64.0_wp)

        ! A zero branch modulus still has a nonzero branch-modulus derivative.
        call material%init(maxwell_shear_t(bulk,1.0_wp,0.0_wp,time_constant), &
            1,status,initial,sensitivity_parameters=[character(len=20) :: "branch_shear_modulus"])
        call require_status(status,cf_ok)
        x(nmandel,1) = huge(1.0_wp)/8.0_wp
        call material%eval(x,y,tangent,status,dt=4.0_wp,dsdp=dsdp)
        call require_status(status,cf_ok)
        call close(y(nmandel,1)/huge(1.0_wp),0.25_wp)
        call close(tangent(nmandel,nmandel,1),2.0_wp)
        call close(dsdp(nmandel,1,1)/(huge(1.0_wp)*quantum),5.0_wp/16.0_wp)
        call material%commit(status)
        call require_status(status,cf_ok)
        call material%export_committed(state,status,previous_input=previous,sensitivity_state=state_derivative)
        call require_status(status,cf_ok)
        if (.not. all(ieee_is_finite(state)) .or. .not. all(ieee_is_finite(state_derivative))) then
            error stop "Nonfinite zero-branch history"
        end if
        call close(maxval(abs(state)),0.0_wp)
        call close(state_derivative(nmandel,1,1)/(huge(1.0_wp)*quantum),5.0_wp/16.0_wp)

        ! With trial overstress=quantum, even normal time fractions can
        ! quantize intermediate products before division restores a normal derivative.
        ! tau=2*t,dt=4*t gives quantum/(9*t); swapping them gives quantum/(18*t).
        x(nmandel,1) = 0.125_wp
        do step = 1,2
            time_constant = real(2*step,wp)*tiny(1.0_wp)
            call material%init(maxwell_shear_t(bulk,1.0_wp,4.0_wp*quantum,time_constant), &
                1,status,initial,sensitivity_parameters=[character(len=15) :: "relaxation_time"])
            call require_status(status,cf_ok)
            call material%eval(x,y,tangent,status,dt=real(6-2*step,wp)*tiny(1.0_wp),dsdp=dsdp)
            call require_status(status,cf_ok)
            call close(y(nmandel,1),0.25_wp)
            call close(tangent(nmandel,nmandel,1),2.0_wp)
            call close(dsdp(nmandel,1,1)/epsilon(1.0_wp),1.0_wp/real(9*step,wp))
        end do

        ! Even the unweighted branch forcing can round to zero: G1=quantum
        ! gives forcing=quantum/4, but tau=dt=quantum gives its tau derivative 1/16.
        call material%init(maxwell_shear_t(bulk,1.0_wp,quantum,quantum), &
            1,status,initial,sensitivity_parameters=[character(len=15) :: "relaxation_time"])
        call require_status(status,cf_ok)
        call material%eval(x,y,tangent,status,dt=quantum,dsdp=dsdp)
        call require_status(status,cf_ok)
        call close(y(nmandel,1),0.25_wp)
        call close(tangent(nmandel,nmandel,1),2.0_wp)
        call close(dsdp(nmandel,1,1),1.0_wp/16.0_wp)
    end subroutine test_subnormal_attenuation_reference

    subroutine test_hydrostatic_history()
        !! Uniform dilation has no deviatoric forcing during loading, hold or reversal.
        real(wp), parameter :: dilations(3) = [0.1_wp,0.1_wp,-0.1_wp]
        type(material_t) :: material
        real(wp) :: initial(nmandel,1), x(nmandel,1), y(nmandel,1)
        real(wp) :: tangent(nmandel,nmandel,1), dsdp(nmandel,nparameters,1)
        real(wp) :: state(nmandel,1), state_derivative(nmandel,nparameters,1), previous(nmandel,1)
        real(wp) :: small_bulk, dilation
        integer :: status, step, i

        small_bulk = epsilon(1.0_wp)**2
        initial = 0.0_wp
        call material%init(maxwell_shear_t(small_bulk,1.0_wp,1.0_wp,1.0_wp), &
            1,status,initial,sensitivity_parameters=[character(len=26) :: &
                "bulk_modulus","equilibrium_shear_modulus","branch_shear_modulus","relaxation_time"])
        call require_status(status,cf_ok)
        do step = 1,size(dilations)
            dilation = dilations(step)
            x = 0.0_wp
            x(1:3,1) = dilation
            call material%eval(x,y,tangent,status,dt=1.0_wp,dsdp=dsdp)
            call require_status(status,cf_ok)
            if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(tangent)) .or. &
                .not. all(ieee_is_finite(dsdp))) error stop "Nonfinite hydrostatic Maxwell response"
            do i = 1,3
                call close(y(i,1)/(small_bulk*dilation),3.0_wp)
                call close(dsdp(i,1,1),3.0_wp*dilation)
            end do
            if (maxval(abs(y(4:nmandel,1))) > 0.0_wp .or. &
                maxval(abs(dsdp(:,2:nparameters,1))) > 0.0_wp .or. &
                maxval(abs(dsdp(4:nmandel,1,1))) > 0.0_wp) then
                error stop "Spurious hydrostatic Maxwell shear response"
            end if
            call material%commit(status)
            call require_status(status,cf_ok)
            call material%export_committed(state,status,previous_input=previous,sensitivity_state=state_derivative)
            call require_status(status,cf_ok)
            if (.not. all(ieee_is_finite(state)) .or. .not. all(ieee_is_finite(state_derivative)) .or. &
                .not. all(ieee_is_finite(previous))) error stop "Nonfinite hydrostatic Maxwell history"
            if (maxval(abs(state)) > 0.0_wp .or. maxval(abs(state_derivative)) > 0.0_wp) then
                error stop "Spurious hydrostatic Maxwell branch history"
            end if
            call close(maxval(abs(previous-x)),0.0_wp)
        end do
    end subroutine test_hydrostatic_history

    subroutine test_path_sensitivities()
        !! Full accepted histories, not single-step parameter perturbations.
        type(material_t) :: tracked, ordinary, single
        type(material_t) :: positive(nparameters), negative(nparameters)
        real(wp) :: initial(nmandel,npoints), x(nmandel,npoints), y(nmandel,npoints)
        real(wp) :: tangent(nmandel,nmandel,npoints)
        real(wp) :: dsdp(nmandel,nparameters,npoints), reference(nmandel,npoints)
        real(wp) :: previous_dsdp(nmandel,nparameters,npoints)
        real(wp) :: response_plus(nmandel,npoints), response_minus(nmandel,npoints)
        real(wp) :: single_stress(nmandel,1), initial_stress(nmandel,npoints)
        real(wp) :: initial_state(nmandel,npoints), initial_state_derivative(nmandel,nparameters,npoints)
        real(wp) :: saved_initial(nmandel,npoints)
        real(wp) :: parameters(nparameters), plus_parameters(nparameters)
        real(wp) :: minus_parameters(nparameters), perturbations(nparameters)
        real(wp) :: dt, derivative, relative_step, error_bound
        integer(lk) :: expected_extra
        integer :: status, k, j, i, step, parameter_id

        parameters = [bulk,equilibrium,branch,relaxation]
        initial = 0.0_wp
        initial(1,1) = 0.01_wp
        initial(2,2) = -0.02_wp
        call tracked%init(maxwell_shear_t(bulk,equilibrium,branch,relaxation), &
            npoints, status, initial, &
            sensitivity_parameters=[character(len=26) :: "relaxation_time", &
                "branch_shear_modulus", "equilibrium_shear_modulus", "bulk_modulus"])
        call require_status(status, cf_ok)
        call ordinary%init(maxwell_shear_t(bulk,equilibrium,branch,relaxation), &
            npoints, status, initial)
        call require_status(status, cf_ok)
        call single%init(maxwell_shear_t(bulk,equilibrium,branch,relaxation), &
            1, status, initial(:,1:1))
        call require_status(status, cf_ok)

        ! Nonzero initial strain is relaxed: q and its initial derivatives are zero.
        call tracked%export_committed(initial_state,status,previous_input=saved_initial, &
            sensitivity_state=initial_state_derivative)
        call require_status(status,cf_ok)
        if (.not. all(ieee_is_finite(initial_state)) .or. &
            .not. all(ieee_is_finite(initial_state_derivative)) .or. .not. all(ieee_is_finite(saved_initial))) then
            error stop "Nonfinite initially relaxed Maxwell history"
        end if
        call close(maxval(abs(initial_state)),0.0_wp)
        call close(maxval(abs(initial_state_derivative)),0.0_wp)
        call close(maxval(abs(saved_initial-initial)),0.0_wp)
        initial_stress = 0.0_wp
        initial_stress(1:3,1) = [17.0_wp/15.0_wp,14.0_wp/15.0_wp,14.0_wp/15.0_wp]
        initial_stress(1:3,2) = [-28.0_wp/15.0_wp,-34.0_wp/15.0_wp,-28.0_wp/15.0_wp]
        call tracked%eval(initial,y,status=status,dt=1.0_wp,dsdp=dsdp)
        call require_status(status,cf_ok)
        if (.not. all(ieee_is_finite(dsdp))) error stop "Nonfinite initially relaxed Maxwell derivatives"
        do i = 1,npoints
            do j = 1,nmandel
                call close(y(j,i),initial_stress(j,i))
            end do
        end do
        call close(maxval(abs(dsdp(:,1:2,:))),0.0_wp)
        call tracked%revert()

        expected_extra = 2_lk*int(nmandel,lk)*int(npoints,lk)* &
                         int(nparameters,lk)*int(storage_size(0.0_wp)/8,lk)
        if (tracked%state_storage_bytes()-ordinary%state_storage_bytes() /= &
            expected_extra) error stop "Incorrect path sensitivity memory"
        if (ordinary%sensitivity_size() /= 0) error stop "Unexpected sensitivity buffers"

        do k = 1, nparameters
            parameter_id = nparameters+1-k
            if (parameter_id == nparameters) then
                relative_step = epsilon(1.0_wp)**(1.0_wp/3.0_wp)
            else
                ! The accepted response is affine in each modulus.
                relative_step = 0.25_wp
            end if
            perturbations(k) = relative_step*parameters(parameter_id)
            plus_parameters = parameters
            minus_parameters = parameters
            plus_parameters(parameter_id) = plus_parameters(parameter_id)+perturbations(k)
            minus_parameters(parameter_id) = minus_parameters(parameter_id)-perturbations(k)
            call positive(k)%init(maxwell_shear_t(plus_parameters(1),plus_parameters(2), &
                plus_parameters(3),plus_parameters(4)), npoints, status, initial)
            call require_status(status, cf_ok)
            call negative(k)%init(maxwell_shear_t(minus_parameters(1),minus_parameters(2), &
                minus_parameters(3),minus_parameters(4)), npoints, status, initial)
            call require_status(status, cf_ok)
        end do

        do step = 1, 4
            dt = 0.2_wp + 0.15_wp*real(step,wp)
            x = initial
            x(1,1) = x(1,1) + 0.018_wp*real(step,wp)
            x(2,1) = x(2,1) - 0.015_wp*real(step,wp)
            x(4,1) = 0.006_wp*real(3-step,wp)
            x(1,2) = x(1,2) + 0.025_wp*real(3-step,wp)
            x(2,2) = x(2,2) + 0.011_wp*real(step,wp)
            x(5,2) = 0.004_wp*real(step-2,wp)
            call tracked%eval(x, y, tangent, status, dt=dt, dsdp=dsdp)
            call require_status(status, cf_ok)
            if (.not. all(ieee_is_finite(dsdp))) error stop "Nonfinite Maxwell path derivatives"
            call ordinary%eval(x, reference, status=status, dt=dt)
            call require_status(status, cf_ok)
            do i = 1, npoints
                do j = 1, nmandel
                    call close(y(j,i), reference(j,i))
                end do
            end do
            if (step == 1) then
                call single%eval(x(:,1:1), single_stress, status=status, dt=dt)
                call require_status(status, cf_ok)
                do j = 1, nmandel
                    call close(y(j,1), single_stress(j,1))
                end do
            end if

            do k = 1, nparameters
                call positive(k)%eval(x, response_plus, status=status, dt=dt)
                call require_status(status, cf_ok)
                call negative(k)%eval(x, response_minus, status=status, dt=dt)
                call require_status(status, cf_ok)
                if (.not. all(ieee_is_finite(response_plus)) .or. .not. all(ieee_is_finite(response_minus))) then
                    error stop "Nonfinite Maxwell parameter secant response"
                end if
                do i = 1, npoints
                    do j = 1, nmandel
                        derivative = (response_plus(j,i)-response_minus(j,i))/ &
                                     (2.0_wp*perturbations(k))
                        error_bound = 32.0_wp*epsilon(1.0_wp)* &
                            maxval(abs(response_plus(:,i))+abs(response_minus(:,i)))/ &
                            (2.0_wp*perturbations(k))
                        if (k == 1) then
                            ! Only tau is nonlinear; four accepted rational
                            ! updates contribute O((h/tau)**2) truncation.
                            error_bound = error_bound+16.0_wp* &
                                (perturbations(k)/relaxation)**2* &
                                maxval(abs(y(:,i)))/relaxation
                        end if
                        if (.not. all(ieee_is_finite([derivative,error_bound]))) then
                            error stop "Nonfinite Maxwell parameter secant or error bound"
                        end if
                        if (abs(derivative-dsdp(j,k,i)) > error_bound) then
                            error stop "Incorrect full-history parameter sensitivity"
                        end if
                    end do
                end do
            end do

            if (step == 2) then
                previous_dsdp = dsdp
                call tracked%revert()
                call ordinary%revert()
                do k = 1, nparameters
                    call positive(k)%revert()
                    call negative(k)%revert()
                end do
                call tracked%eval(x, y, tangent, status, dt=dt, dsdp=dsdp)
                call require_status(status, cf_ok)
                do k = 1, nparameters
                    do i = 1, npoints
                        do j = 1, nmandel
                            call close(dsdp(j,k,i), previous_dsdp(j,k,i))
                        end do
                    end do
                end do
                call ordinary%eval(x, reference, status=status, dt=dt)
                call require_status(status, cf_ok)
                do k = 1, nparameters
                    call positive(k)%eval(x, response_plus, status=status, dt=dt)
                    call require_status(status, cf_ok)
                    call negative(k)%eval(x, response_minus, status=status, dt=dt)
                    call require_status(status, cf_ok)
                end do
            end if

            call tracked%commit(status)
            call require_status(status, cf_ok)
            call ordinary%commit(status)
            call require_status(status, cf_ok)
            do k = 1, nparameters
                call positive(k)%commit(status)
                call require_status(status, cf_ok)
                call negative(k)%commit(status)
                call require_status(status, cf_ok)
            end do
        end do
        call tracked%commit(status)
        call require_status(status, cf_no_trial)
    end subroutine test_path_sensitivities

    subroutine test_invalid_trials()
        type(material_t) :: material
        real(wp) :: initial(nmandel,npoints), x(nmandel,npoints)
        real(wp) :: y(nmandel,npoints), invalid_x(nmandel-1,npoints)
        integer :: status

        initial = 0.0_wp
        call material%init(maxwell_shear_t(bulk,equilibrium,branch,relaxation), &
                           npoints, status, initial)
        call require_status(status, cf_ok)
        x = 0.0_wp
        x(1,:) = 0.02_wp
        call material%eval(x, y, status=status)
        call require_status(status, cf_invalid_time)
        call material%eval(x, y, status=status, dt=0.0_wp)
        call require_status(status, cf_invalid_time)
        x(1,2) = ieee_value(0.0_wp,ieee_quiet_nan)
        call material%eval(x, y, status=status, dt=0.3_wp)
        call require_status(status, cf_evaluation_failed)
        if (material%last_failure_point() /= 2) error stop "Wrong failing point"
        call material%commit(status)
        call require_status(status, cf_no_trial)

        x(1,2) = 0.02_wp
        call material%eval(x, y, status=status, dt=0.3_wp)
        call require_status(status, cf_ok)
        invalid_x = 0.0_wp
        call material%eval(invalid_x, y, status=status, dt=0.3_wp)
        call require_status(status, cf_invalid_shape)
        call material%init(maxwell_shear_t(-1.0_wp,equilibrium,branch,relaxation), &
                           npoints, status, initial)
        call require_status(status, cf_invalid_argument)
        call material%eval(x, y, status=status, dt=0.3_wp)
        call require_status(status, cf_ok)
        call material%commit(status)
        call require_status(status, cf_ok)
    end subroutine test_invalid_trials
end program test_maxwell_shear

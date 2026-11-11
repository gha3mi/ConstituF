! SPDX-License-Identifier: Apache-2.0

program test_maxwell_multi
    !! Test convention: moduli, stress and overstress in Pa. Relaxation times
    !! and dt are in s. Explicit rescaling cases state their factors.
    !!
    !! Independent Prony-branch recurrence and accepted-path derivatives.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf, only: wp, material_t, verify_tangent, cf_ok, cf_invalid_argument
    use constituf_maxwell_multi, only: maxwell_multi_t
    use constituf_maxwell_shear, only: maxwell_shear_t
    implicit none (type, external)

    integer, parameter :: nmandel = 6, nbranches = 2, nreference = 4
    integer, parameter :: nstate = nmandel*nbranches, nparameters = 2+2*nbranches
    real(wp), parameter :: bulk = 80.0_wp, equilibrium = 15.0_wp
    real(wp), parameter :: gs(nbranches) = [30.0_wp,12.0_wp], taus(nbranches) = [0.5_wp,2.0_wp]
    real(wp), parameter :: q_factor(nbranches,2) = reshape( &
        [40.0_wp,64.0_wp/3.0_wp,80.0_wp/3.0_wp,512.0_wp/27.0_wp],[nbranches,2])
    real(wp), parameter :: q_modulus(nbranches,2) = reshape( &
        [4.0_wp/3.0_wp,16.0_wp/9.0_wp,8.0_wp/9.0_wp,128.0_wp/81.0_wp],[nbranches,2])
    real(wp), parameter :: q_time(nbranches,2) = reshape( &
        [80.0_wp/3.0_wp,32.0_wp/27.0_wp,320.0_wp/9.0_wp,512.0_wp/243.0_wp],[nbranches,2])
    type(material_t) :: material, scalar, tracked, positive, negative
    real(wp) :: initial(nmandel,2), x(nmandel,2), y(nmandel,2), jac(nmandel,nmandel,2), dsdp(nmandel,nparameters,2)
    real(wp) :: single(nmandel,2), yp(nmandel,2), ym(nmandel,2)
    real(wp) :: gp(nbranches), gm(nbranches), tp(nbranches), tm(nbranches)
    real(wp) :: reference_dsdp(nmandel,nreference,2), state(nstate,2), previous(nmandel,2)
    real(wp) :: state_derivative(nstate,nreference,2)
    real(wp) :: tracked_state_derivative(nstate,nparameters,2), state_plus(nstate,2), state_minus(nstate,2)
    real(wp) :: numerical(nmandel,2), state_numerical(nstate,2), deviation(nmandel)
    real(wp) :: error, h, dt, error_bound, state_error_bound, trace_strain, bulk_derivative
    real(wp) :: bulk_plus, bulk_minus, equilibrium_plus, equilibrium_minus
    integer :: status, step, k, i, component, branch_index
    character(len=26) :: names(nparameters)

    names = [character(len=26) :: "bulk_modulus","equilibrium_shear_modulus", &
        "branch_shear_modulus_1","relaxation_time_1","branch_shear_modulus_2","relaxation_time_2"]
    initial = 0.0_wp
    x(:,1) = [0.02_wp,-0.01_wp,-0.01_wp,0.005_wp,0.0_wp,0.0_wp]
    x(:,2) = [-0.006_wp,0.004_wp,0.002_wp,0.0_wp,0.003_wp,-0.002_wp]
    call material%init(maxwell_multi_t(bulk,equilibrium,gs,taus),2,status)
    if (status /= cf_invalid_argument) error stop "Missing initial strain accepted"
    call material%init(maxwell_multi_t(bulk,equilibrium,gs,taus),2,status,initial, &
        sensitivity_parameters=[character(len=24) :: "relaxation_time_2", &
            "branch_shear_modulus_1","relaxation_time_1","branch_shear_modulus_2"])
    if (status /= cf_ok) error stop "Multi-branch Maxwell initialization"
    if (material%state_size() /= nstate) error stop "Incorrect branch state size"
    if (.not. material%field_matches("state","overstress.branches",&
        "custom.mandel_branches","internal_stress","spatial",[nmandel,nbranches])) then
        error stop "Incorrect Maxwell packed state contract"
    end if
    ! Unequal branch attenuation is [2/3,8/9]. Exact ramp/hold references
    ! check both packed branch histories and deliberately reordered indices.
    do step = 1,2
        call material%eval(x,y,jac,status,dt=0.25_wp,dsdp=reference_dsdp)
        if (status /= cf_ok) error stop "Multi-branch Maxwell reference step"
        if (step == 2) then
            call material%revert()
            call material%eval(x,y,jac,status,dt=0.25_wp,dsdp=reference_dsdp)
            if (status /= cf_ok) error stop "Maxwell branch revert"
        end if
        do i = 1,2
            do component = 1,nmandel
                call close(y(component,i),(30.0_wp+sum(q_factor(:,step)))*x(component,i))
                call close(reference_dsdp(component,1,i),q_time(2,step)*x(component,i))
                call close(reference_dsdp(component,2,i),q_modulus(1,step)*x(component,i))
                call close(reference_dsdp(component,3,i),q_time(1,step)*x(component,i))
                call close(reference_dsdp(component,4,i),q_modulus(2,step)*x(component,i))
            end do
        end do
        call close(jac(1,1,1),1268.0_wp/9.0_wp)
        call close(jac(1,2,1),446.0_wp/9.0_wp)
        call close(jac(4,4,1),274.0_wp/3.0_wp)
        call material%commit(status)
        if (status /= cf_ok) error stop "Maxwell branch commit"
        call material%export_committed(state,status,previous_input=previous, &
            sensitivity_state=state_derivative)
        if (status /= cf_ok) error stop "Maxwell branch state export"
        if (.not. all(ieee_is_finite(previous))) error stop "Nonfinite Maxwell previous input"
        call close(maxval(abs(previous-x)),0.0_wp)
        do i = 1,2
            do component = 1,nmandel
                call close(state(component,i),q_factor(1,step)*x(component,i))
                call close(state(nmandel+component,i),q_factor(2,step)*x(component,i))
                call close(state_derivative(component,1,i),0.0_wp)
                call close(state_derivative(component,2,i),q_modulus(1,step)*x(component,i))
                call close(state_derivative(component,3,i),q_time(1,step)*x(component,i))
                call close(state_derivative(component,4,i),0.0_wp)
                call close(state_derivative(nmandel+component,1,i),q_time(2,step)*x(component,i))
                call close(state_derivative(nmandel+component,2,i),0.0_wp)
                call close(state_derivative(nmandel+component,3,i),0.0_wp)
                call close(state_derivative(nmandel+component,4,i),q_modulus(2,step)*x(component,i))
            end do
        end do
    end do

    call test_shifted_branch_bounds()
    call test_mixed_relaxation(0.25_wp*sqrt(tiny(1.0_wp)))
    call test_mixed_relaxation(1.0_wp/huge(1.0_wp),huge(1.0_wp)/8.0_wp)
    call test_subnormal_branch_reference()
    call test_hydrostatic_history()

    ! A one-branch multi model must reproduce the existing single-branch law.
    call scalar%init(maxwell_shear_t(bulk,equilibrium,gs(1),taus(1)),2,status,initial)
    if (status /= cf_ok) error stop "Single branch setup"
    call positive%init(maxwell_multi_t(bulk,equilibrium,gs(1:1),taus(1:1)),2,status,initial)
    if (status /= cf_ok) error stop "One branch multi setup"
    call scalar%eval(x,single,status=status,dt=0.4_wp)
    if (status /= cf_ok) error stop "Single branch response"
    call positive%eval(x,yp,status=status,dt=0.4_wp)
    if (status /= cf_ok) error stop "One branch multi response"
    if (.not. all(ieee_is_finite(single)) .or. .not. all(ieee_is_finite(yp)) .or. &
        maxval(abs(single-yp)) > 2000.0_wp*epsilon(1.0_wp)) error stop "Maxwell N=1 reduction"

    call tracked%init(maxwell_multi_t(bulk,equilibrium,gs,taus),2,status,initial,sensitivity_parameters=names)
    if (status /= cf_ok) error stop "Maxwell selected sensitivity setup"
    call material%init(maxwell_multi_t(bulk,equilibrium,gs,taus),2,status,initial)
    if (status /= cf_ok) error stop "Maxwell verification setup"
    do step = 1, 3
        dt = 0.15_wp+0.1_wp*real(step,wp)
        x(1,:) = x(1,:)+0.003_wp*real(step,wp)
        x(4,:) = x(4,:)-0.001_wp*real(step,wp)
        call tracked%eval(x,y,jac,status,dt=dt,dsdp=dsdp)
        if (status /= cf_ok) error stop "Maxwell tracked evaluation"

        call tracked%commit(status)
        if (status /= cf_ok) error stop "Maxwell tracked commit"
        call tracked%export_committed(state,status,previous_input=previous, &
            sensitivity_state=tracked_state_derivative)
        if (status /= cf_ok) error stop "Maxwell complete history export"
        if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(dsdp)) .or. &
            .not. all(ieee_is_finite(state)) .or. .not. all(ieee_is_finite(tracked_state_derivative)) .or. &
            .not. all(ieee_is_finite(previous))) error stop "Nonfinite Maxwell path reference"
        call close(maxval(abs(previous-x)),0.0_wp)
        do i = 1,2
            trace_strain = sum(x(1:3,i))
            deviation = x(:,i)
            deviation(1:3) = deviation(1:3)-trace_strain/3.0_wp
            call close(sum(y(1:3,i))/3.0_wp,bulk*trace_strain)
            do component = 1,nmandel
                bulk_derivative = 0.0_wp
                if (component<=3) bulk_derivative = trace_strain
                call close(dsdp(component,1,i),bulk_derivative)
                call close(dsdp(component,2,i),2.0_wp*deviation(component))
            end do
        end do

        do k = 1,nparameters
            bulk_plus = bulk
            bulk_minus = bulk
            equilibrium_plus = equilibrium
            equilibrium_minus = equilibrium
            gp = gs
            gm = gs
            tp = taus
            tm = taus
            branch_index = 0
            select case(k)
            case(1)
                h = 0.25_wp*bulk
                bulk_plus = bulk+h
                bulk_minus = bulk-h
            case(2)
                h = 0.25_wp*equilibrium
                equilibrium_plus = equilibrium+h
                equilibrium_minus = equilibrium-h
            case default
                branch_index = (k-1)/2
                if (mod(k,2)==1) then
                    ! Accepted histories are affine in each modulus.
                    h = 0.25_wp*gs(branch_index)
                    gp(branch_index) = gp(branch_index)+h
                    gm(branch_index) = gm(branch_index)-h
                else
                    h = epsilon(1.0_wp)**(1.0_wp/3.0_wp)*taus(branch_index)
                    tp(branch_index) = tp(branch_index)+h
                    tm(branch_index) = tm(branch_index)-h
                end if
            end select
            call positive%init(maxwell_multi_t(bulk_plus,equilibrium_plus,gp,tp),2,status,initial)
            if (status /= cf_ok) error stop "Maxwell positive parameter setup"
            call negative%init(maxwell_multi_t(bulk_minus,equilibrium_minus,gm,tm),2,status,initial)
            if (status /= cf_ok) error stop "Maxwell negative parameter setup"
            ! Replay the complete accepted path for both stress and branch-state derivatives.
            call replay(positive,step,yp,state_plus)
            call replay(negative,step,ym,state_minus)
            if (.not. all(ieee_is_finite(yp)) .or. .not. all(ieee_is_finite(ym)) .or. &
                .not. all(ieee_is_finite(state_plus)) .or. .not. all(ieee_is_finite(state_minus))) then
                error stop "Nonfinite Maxwell perturbed history"
            end if
            numerical = (yp-ym)/(2.0_wp*h)
            state_numerical = (state_plus-state_minus)/(2.0_wp*h)
            error_bound = 32.0_wp*epsilon(1.0_wp)*maxval(abs(yp)+abs(ym))/(2.0_wp*h)
            state_error_bound = 32.0_wp*epsilon(1.0_wp)*maxval(abs(state_plus)+abs(state_minus))/(2.0_wp*h)
            if (k>2 .and. mod(k,2)==0) then
                error_bound = error_bound+16.0_wp*(h/taus(branch_index))**2*maxval(abs(y))/taus(branch_index)
                state_error_bound = state_error_bound+16.0_wp*(h/taus(branch_index))**2* &
                    maxval(abs(state))/taus(branch_index)
            end if
            if (.not. all(ieee_is_finite(numerical)) .or. .not. all(ieee_is_finite(state_numerical)) .or. &
                .not. all(ieee_is_finite([error_bound,state_error_bound]))) then
                error stop "Nonfinite Maxwell parameter secant or bound"
            end if
            if (maxval(abs(numerical-dsdp(:,k,:))) > error_bound) then
                error stop "Multi-branch Maxwell full-history sensitivity mismatch"
            end if
            if (maxval(abs(state_numerical-tracked_state_derivative(:,k,:))) > state_error_bound) then
                error stop "Multi-branch Maxwell state sensitivity mismatch"
            end if
        end do
        call material%eval(x,yp,status=status,dt=dt)
        if (status /= cf_ok) error stop "Maxwell ordinary evaluation"
        if (.not. all(ieee_is_finite(yp)) .or. maxval(abs(yp-y)) > 5000.0_wp*epsilon(1.0_wp)) then
            error stop "Maxwell tracked stress differs"
        end if
        call material%commit(status)
        if (status /= cf_ok) error stop "Maxwell ordinary commit"
    end do
    h = 0.25_wp*maxval(abs(x))
    call verify_tangent(material,x,error,status,dt=0.3_wp,step=[(h,i=1,nmandel)], &
        absolute_tolerance=32.0_wp*epsilon(1.0_wp)*maxval(abs(y))/h)
    if (status /= cf_ok .or. .not. ieee_is_finite(error) .or. error > 1.0_wp) then
        error stop "Maxwell multi tangent verification"
    end if

    call positive%init(maxwell_multi_t(bulk,equilibrium,[gs(1),-1.0_wp],taus),2,status,initial)
    if (status /= cf_invalid_argument) error stop "Negative branch modulus accepted"
    call positive%init(maxwell_multi_t(bulk,equilibrium,gs,[taus(1),0.0_wp]),2,status,initial)
    if (status /= cf_invalid_argument) error stop "Zero branch relaxation time accepted"
    call positive%init(maxwell_multi_t(bulk,equilibrium,gs,taus(1:1)),2,status,initial)
    if (status /= cf_invalid_argument) error stop "Mismatched branch array sizes accepted"

contains

    subroutine close(actual, expected)
        real(wp), intent(in) :: actual, expected

        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) then
            error stop "Nonfinite Maxwell multi comparison"
        end if
        if (abs(actual-expected) > 512.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(expected))) then
            print *, "Expected",expected,"received",actual
            error stop "Maxwell independent recurrence discrepancy"
        end if
    end subroutine close

    subroutine test_shifted_branch_bounds()
        !! Parameter indices count branches, independently of allocatable bounds.
        type(maxwell_multi_t) :: law
        type(material_t) :: model
        real(wp) :: initial_input(nmandel,1), input(nmandel,1), response(nmandel,1)
        real(wp) :: derivative(nmandel,2,1), history(nstate,1), saved_input(nmandel,1)
        real(wp) :: history_derivative(nstate,2,1)
        integer :: result

        law%bulk_modulus = bulk
        law%equilibrium_shear_modulus = equilibrium
        allocate(law%branch_shear_moduli(-1:0),law%relaxation_times(3:4))
        law%branch_shear_moduli = gs
        law%relaxation_times = taus
        initial_input = 0.0_wp
        call model%init(law,1,result,initial_input, &
            sensitivity_parameters=[character(len=24) :: "relaxation_time_2","branch_shear_modulus_1"])
        if (result /= cf_ok) error stop "Shifted branch arrays initialization"
        input = 0.0_wp
        input(nmandel,1) = 0.125_wp
        call model%eval(input,response,status=result,dt=0.25_wp,dsdp=derivative)
        if (result /= cf_ok) error stop "Shifted branch arrays evaluation"
        call close(response(nmandel,1),137.0_wp/12.0_wp)
        call close(derivative(nmandel,1,1),4.0_wp/27.0_wp)
        call close(derivative(nmandel,2,1),1.0_wp/6.0_wp)
        call model%commit(result)
        if (result /= cf_ok) error stop "Shifted branch arrays commit"
        call model%export_committed(history,result,previous_input=saved_input, &
            sensitivity_state=history_derivative)
        if (result /= cf_ok) error stop "Shifted branch arrays export"
        call close(history(nmandel,1),5.0_wp)
        call close(history(nstate,1),8.0_wp/3.0_wp)
        call close(history_derivative(nmandel,1,1),0.0_wp)
        call close(history_derivative(nmandel,2,1),1.0_wp/6.0_wp)
        call close(history_derivative(nstate,1,1),4.0_wp/27.0_wp)
        call close(history_derivative(nstate,2,1),0.0_wp)
    end subroutine test_shifted_branch_bounds

    subroutine test_mixed_relaxation(time_constant,large_time_step)
        !! A zero-modulus ordinary branch has zero state but a nonzero modulus
        !! derivative. The second branch retains tiny, finite Maxwell stress.
        real(wp), intent(in) :: time_constant
        real(wp), intent(in), optional :: large_time_step
        type(material_t) :: model
        type(maxwell_multi_t) :: law
        real(wp) :: input(nmandel,1), initial_input(nmandel,1), response(nmandel,1)
        real(wp) :: tangent(nmandel,nmandel,1), derivative(nmandel,nreference,1)
        real(wp) :: history(nstate,1), history_derivative(nstate,nreference,1)
        real(wp) :: saved_input(nmandel,1)
        real(wp) :: subnormal_bound
        integer :: result

        initial_input = 0.0_wp
        call model%init(maxwell_multi_t(bulk,0.0_wp,[0.0_wp,30.0_wp],[0.5_wp,time_constant]), &
            1,result,initial_input,sensitivity_parameters=[character(len=24) :: &
                "relaxation_time_2","branch_shear_modulus_1", &
                "relaxation_time_1","branch_shear_modulus_2"])
        if (result /= cf_ok) error stop "Mixed relaxation initialization"
        input = 0.0_wp
        input(nmandel,1) = 0.125_wp
        call model%eval(input,response,tangent,result,dt=2.0_wp,dsdp=derivative)
        if (result /= cf_ok) error stop "Mixed finite relaxation response"
        call close(response(nmandel,1)/time_constant,15.0_wp/4.0_wp)
        call close(tangent(nmandel,nmandel,1)/time_constant,30.0_wp)
        call close(derivative(nmandel,1,1),15.0_wp/4.0_wp)
        call close(derivative(nmandel,2,1),0.05_wp)
        call close(derivative(nmandel,3,1),0.0_wp)
        call close(derivative(nmandel,4,1)/time_constant,0.125_wp)
        call model%commit(result)
        if (result /= cf_ok) error stop "Mixed relaxation commit"
        call model%export_committed(history,result,previous_input=saved_input, &
            sensitivity_state=history_derivative)
        if (result /= cf_ok) error stop "Mixed relaxation export"
        if (.not. all(ieee_is_finite(history)) .or. .not. all(ieee_is_finite(history_derivative))) then
            error stop "Nonfinite mixed relaxation history"
        end if
        call close(maxval(abs(history(1:nmandel,1))),0.0_wp)
        call close(history(nstate,1)/time_constant,15.0_wp/4.0_wp)
        call close(history_derivative(nmandel,2,1),0.05_wp)
        call close(maxval(abs(history_derivative(1:nmandel,1,1))),0.0_wp)
        call close(maxval(abs(history_derivative(1:nmandel,3:4,1))),0.0_wp)
        call close(history_derivative(nstate,1,1),15.0_wp/4.0_wp)
        call close(history_derivative(nstate,4,1)/time_constant,0.125_wp)
        call close(maxval(abs(history_derivative(nmandel+1:nstate,2:3,1))),0.0_wp)

        if (.not. present(large_time_step)) return
        ! G2=dt=L and tau2=t: tau2/dt rounds to zero, but sigma6=t/4 and C66=2*t are representable.
        ! Independent array bounds exercise the second branch's weighted update and parameter lookup.
        law%bulk_modulus = bulk
        law%equilibrium_shear_modulus = 0.0_wp
        allocate(law%branch_shear_moduli(-1:0),law%relaxation_times(3:4))
        law%branch_shear_moduli = [0.0_wp,large_time_step]
        law%relaxation_times = [0.5_wp,time_constant]
        subnormal_bound = 8.0_wp*nearest(0.0_wp,1.0_wp)/time_constant
        call model%init(law,1,result,initial_input,sensitivity_parameters=[character(len=24) :: &
            "relaxation_time_2","branch_shear_modulus_1","relaxation_time_1","branch_shear_modulus_2"])
        if (result /= cf_ok) error stop "Weighted multi-branch initialization"
        call model%eval(input,response,tangent,result,dt=large_time_step,dsdp=derivative)
        if (result /= cf_ok) error stop "Weighted multi-branch response"
        if (.not. all(ieee_is_finite(response)) .or. .not. all(ieee_is_finite(tangent)) .or. &
            .not. all(ieee_is_finite(derivative)) .or. .not. ieee_is_finite(subnormal_bound)) then
            error stop "Nonfinite weighted multi-branch reference"
        end if
        if (abs(response(nmandel,1)/time_constant-0.25_wp) > subnormal_bound .or. &
            abs(tangent(nmandel,nmandel,1)/time_constant-2.0_wp) > subnormal_bound) then
            error stop "Maxwell multi discarded representable weighted stress or tangent"
        end if
        call close(derivative(nmandel,1,1),0.25_wp)
        call model%commit(result)
        if (result /= cf_ok) error stop "Weighted multi-branch commit"
        call model%export_committed(history,result,previous_input=saved_input,sensitivity_state=history_derivative)
        if (result /= cf_ok) error stop "Weighted multi-branch export"
        if (.not. all(ieee_is_finite(history)) .or. .not. all(ieee_is_finite(history_derivative))) then
            error stop "Nonfinite weighted multi-branch history"
        end if
        call close(maxval(abs(history(1:nmandel,1))),0.0_wp)
        if (abs(history(nstate,1)/time_constant-0.25_wp) > subnormal_bound) then
            error stop "Maxwell multi discarded representable weighted history"
        end if
        call close(history_derivative(nstate,1,1),0.25_wp)
        call close(maxval(abs(history_derivative(1:nmandel,1,1))),0.0_wp)

        ! A short load creates q2=L*t/8 and q2_tau=L/16. A long hold then gives q2_tau/t=3/16:
        ! propagated history contributes 1/16 and the direct derivative contributes 1/8.
        ! Prescribed strain and dt remain fixed when differentiating tau2.
        call model%init(law,1,result,initial_input,sensitivity_parameters=[character(len=24) :: &
            "relaxation_time_2","branch_shear_modulus_1","relaxation_time_1","branch_shear_modulus_2"])
        if (result /= cf_ok) error stop "Weighted multi-branch path initialization"
        input(nmandel,1) = time_constant/8.0_wp
        call model%eval(input,response,status=result,dt=time_constant,dsdp=derivative)
        if (result /= cf_ok) error stop "Weighted multi-branch short load"
        call close(response(nmandel,1)/(large_time_step*time_constant),0.125_wp)
        call close(derivative(nmandel,1,1)/large_time_step,1.0_wp/16.0_wp)
        call model%commit(result)
        if (result /= cf_ok) error stop "Weighted multi-branch short-load commit"
        call model%eval(input,response,status=result,dt=large_time_step,dsdp=derivative)
        if (result /= cf_ok) error stop "Weighted multi-branch long hold"
        call close(response(nmandel,1),0.0_wp)
        if (.not. ieee_is_finite(derivative(nmandel,1,1)) .or. &
            abs(derivative(nmandel,1,1)/time_constant-3.0_wp/16.0_wp) > subnormal_bound) then
            error stop "Maxwell multi lost propagated short-relaxation sensitivity"
        end if
        call model%commit(result)
        if (result /= cf_ok) error stop "Weighted multi-branch hold commit"
        call model%export_committed(history,result,previous_input=saved_input,sensitivity_state=history_derivative)
        if (result /= cf_ok) error stop "Weighted multi-branch hold export"
        if (.not. all(ieee_is_finite(history)) .or. .not. all(ieee_is_finite(history_derivative))) then
            error stop "Nonfinite weighted multi-branch hold history"
        end if
        call close(maxval(abs(history)),0.0_wp)
        if (abs(history_derivative(nstate,1,1)/time_constant-3.0_wp/16.0_wp) > subnormal_bound) then
            error stop "Maxwell multi lost committed short-relaxation sensitivity"
        end if
    end subroutine test_mixed_relaxation

    subroutine test_subnormal_branch_reference()
        !! Independent bounds and reordered parameters exercise the second
        !! branch without letting the zero first branch contaminate its history.
        type(material_t) :: model
        type(maxwell_multi_t) :: law
        real(wp) :: input(nmandel,1), initial_input(nmandel,1), response(nmandel,1)
        real(wp) :: tangent(nmandel,nmandel,1), derivative(nmandel,nreference,1)
        real(wp) :: history(nstate,1), history_derivative(nstate,nreference,1), saved_input(nmandel,1)
        real(wp) :: quantum, branch_modulus, stress_scale, step_size
        integer :: result, load_step
        character(len=24), parameter :: selected(nreference) = [character(len=24) :: &
            "relaxation_time_2","branch_shear_modulus_1","relaxation_time_1","branch_shear_modulus_2"]

        quantum = nearest(0.0_wp,1.0_wp)
        branch_modulus = huge(1.0_wp)/8.0_wp
        stress_scale = branch_modulus*quantum
        law%bulk_modulus = bulk
        law%equilibrium_shear_modulus = 0.0_wp
        allocate(law%branch_shear_moduli(-1:0),law%relaxation_times(3:4))
        law%branch_shear_moduli = [0.0_wp,branch_modulus]
        law%relaxation_times = [0.5_wp,5.0_wp*quantum]
        initial_input = 0.0_wp
        input = 0.0_wp
        input(nmandel,1) = 0.125_wp
        call model%init(law,1,result,initial_input,sensitivity_parameters=selected)
        if (result /= cf_ok) error stop "Subnormal branch initialization"

        ! A2=5*quantum/4 to working accuracy. Rounded A2=quantum must not
        ! quantize its normal weighted stress, tangent or accepted derivative.
        do load_step = 1,2
            call model%eval(input,response,tangent,result,dt=4.0_wp,dsdp=derivative)
            if (result /= cf_ok) error stop "Subnormal branch evaluation"
            if (.not. all(ieee_is_finite(response)) .or. .not. all(ieee_is_finite(tangent)) .or. &
                .not. all(ieee_is_finite(derivative))) error stop "Nonfinite subnormal branch response"
            call close(tangent(nmandel,nmandel,1)/stress_scale,5.0_wp/2.0_wp)
            call close(derivative(nmandel,2,1),0.25_wp/9.0_wp**load_step)
            call close(derivative(nmandel,3,1),0.0_wp)
            if (load_step == 1) then
                call close(response(nmandel,1)/stress_scale,5.0_wp/16.0_wp)
                call close(derivative(nmandel,1,1)/branch_modulus,1.0_wp/16.0_wp)
            else
                if (abs(response(nmandel,1)) > 0.0_wp) error stop "Unrepresentable relaxed branch stress"
                call close(derivative(nmandel,1,1)/stress_scale,5.0_wp/32.0_wp)
            end if
            call model%commit(result)
            if (result /= cf_ok) error stop "Subnormal branch commit"
            call model%export_committed(history,result,previous_input=saved_input,sensitivity_state=history_derivative)
            if (result /= cf_ok) error stop "Subnormal branch export"
            if (.not. all(ieee_is_finite(history)) .or. .not. all(ieee_is_finite(history_derivative)) .or. &
                .not. all(ieee_is_finite(saved_input))) error stop "Nonfinite subnormal branch history"
            call close(maxval(abs(saved_input-input)),0.0_wp)
            call close(maxval(abs(history(1:nmandel,1))),0.0_wp)
            call close(maxval(abs(history_derivative(1:nmandel,1,1))),0.0_wp)
            call close(maxval(abs(history_derivative(nmandel+1:nstate,2:3,1))),0.0_wp)
            call close(history_derivative(nmandel,2,1),0.25_wp/9.0_wp**load_step)
            if (load_step == 1) then
                call close(history(nstate,1)/stress_scale,5.0_wp/16.0_wp)
                call close(history_derivative(nstate,1,1)/branch_modulus,1.0_wp/16.0_wp)
            else
                if (abs(history(nstate,1)) > 0.0_wp) error stop "Unrepresentable committed branch stress"
                call close(history_derivative(nstate,1,1)/stress_scale,5.0_wp/32.0_wp)
            end if
        end do

        ! A subnormal attenuation does not imply that G2*tau2 is representable.
        law%branch_shear_moduli(0) = huge(1.0_wp)/2.0_wp
        law%relaxation_times(4) = 3.0_wp
        call model%init(law,1,result,initial_input,sensitivity_parameters=selected)
        if (result /= cf_ok) error stop "Balanced branch initialization"
        call model%eval(input,response,tangent,result,dt=huge(1.0_wp),dsdp=derivative)
        if (result /= cf_ok) error stop "Balanced branch evaluation"
        call close(response(nmandel,1),3.0_wp/8.0_wp)
        call close(tangent(nmandel,nmandel,1),3.0_wp)
        call close(derivative(nmandel,1,1),1.0_wp/8.0_wp)

        ! The second branch forcing is 2*huge. Both small and normal
        ! attenuation recover finite response and accepted history derivatives.
        input(nmandel,1) = 2.0_wp
        do load_step = 1,2
            if (load_step == 1) then
                law%relaxation_times(4) = 3.0_wp
                step_size = huge(1.0_wp)
            else
                law%relaxation_times(4) = 1.0_wp
                step_size = 1.0_wp
            end if
            call model%init(law,1,result,initial_input,sensitivity_parameters=selected)
            if (result /= cf_ok) error stop "Large-forcing branch initialization"
            call model%eval(input,response,tangent,result,dt=step_size,dsdp=derivative)
            if (result /= cf_ok) error stop "Large-forcing branch evaluation"
            if (.not. all(ieee_is_finite(response)) .or. .not. all(ieee_is_finite(tangent)) .or. &
                .not. all(ieee_is_finite(derivative))) error stop "Nonfinite large-forcing branch response"
            if (load_step == 1) then
                call close(response(nmandel,1),6.0_wp)
                call close(tangent(nmandel,nmandel,1),3.0_wp)
                call close(derivative(nmandel,1,1),2.0_wp)
                call close(derivative(nmandel,2,1)*huge(1.0_wp),2.0_wp)
                call close(derivative(nmandel,4,1)*huge(1.0_wp),12.0_wp)
            else
                call close(response(nmandel,1)/huge(1.0_wp),1.0_wp)
                call close(tangent(nmandel,nmandel,1)/huge(1.0_wp),0.5_wp)
                call close(derivative(nmandel,1,1)/huge(1.0_wp),0.5_wp)
                call close(derivative(nmandel,2,1),4.0_wp/3.0_wp)
                call close(derivative(nmandel,4,1),2.0_wp)
            end if
            call close(derivative(nmandel,3,1),0.0_wp)
            call model%commit(result)
            if (result /= cf_ok) error stop "Large-forcing branch commit"
            call model%export_committed(history,result,previous_input=saved_input,sensitivity_state=history_derivative)
            if (result /= cf_ok) error stop "Large-forcing branch export"
            if (.not. all(ieee_is_finite(history)) .or. .not. all(ieee_is_finite(history_derivative)) .or. &
                .not. all(ieee_is_finite(saved_input))) error stop "Nonfinite large-forcing branch history"
            if (maxval(abs(history(1:nmandel,1))) > 0.0_wp .or. &
                maxval(abs(history_derivative(1:nmandel,1,1))) > 0.0_wp .or. &
                maxval(abs(history_derivative(nmandel+1:nstate,2:3,1))) > 0.0_wp) then
                error stop "Large forcing contaminated another Maxwell branch"
            end if
            call close(history(nstate,1)/response(nmandel,1),1.0_wp)
            call close(history_derivative(nstate,1,1)/derivative(nmandel,1,1),1.0_wp)
            call close(history_derivative(nmandel,2,1)/derivative(nmandel,2,1),1.0_wp)
            call close(history_derivative(nstate,4,1)/derivative(nmandel,4,1),1.0_wp)
            call close(maxval(abs(saved_input-input)),0.0_wp)
        end do
        input(nmandel,1) = 0.125_wp

        ! The short-step complementary fraction also has to retain its weight.
        law%branch_shear_moduli(0) = branch_modulus
        law%relaxation_times(4) = 4.0_wp
        call model%init(law,1,result,initial_input,sensitivity_parameters=selected)
        if (result /= cf_ok) error stop "Complementary branch initialization"
        call model%eval(input,response,tangent,result,dt=5.0_wp*quantum,dsdp=derivative)
        if (result /= cf_ok) error stop "Complementary branch evaluation"
        call close(response(nmandel,1)/branch_modulus,0.25_wp)
        call close(tangent(nmandel,nmandel,1)/branch_modulus,2.0_wp)
        call close(derivative(nmandel,1,1)/stress_scale,5.0_wp/64.0_wp)

        ! At G2=0, the branch state is zero but its modulus derivative is not.
        law%branch_shear_moduli(0) = 0.0_wp
        law%equilibrium_shear_modulus = 1.0_wp
        law%relaxation_times(4) = 5.0_wp*quantum
        input(nmandel,1) = huge(1.0_wp)/8.0_wp
        call model%init(law,1,result,initial_input,sensitivity_parameters=selected)
        if (result /= cf_ok) error stop "Zero weighted branch initialization"
        call model%eval(input,response,tangent,result,dt=4.0_wp,dsdp=derivative)
        if (result /= cf_ok) error stop "Zero weighted branch evaluation"
        call close(response(nmandel,1)/huge(1.0_wp),0.25_wp)
        call close(tangent(nmandel,nmandel,1),2.0_wp)
        call close(derivative(nmandel,4,1)/(huge(1.0_wp)*quantum),5.0_wp/16.0_wp)

        ! Normal time fractions can still quantize a subnormal forcing before
        ! division. With forcing=quantum, the two tau derivatives are eps/9,eps/18.
        law%branch_shear_moduli(0) = 4.0_wp*quantum
        input(nmandel,1) = 0.125_wp
        do load_step = 1,2
            law%relaxation_times(4) = real(2*load_step,wp)*tiny(1.0_wp)
            call model%init(law,1,result,initial_input,sensitivity_parameters=selected)
            if (result /= cf_ok) error stop "Subnormal forcing branch initialization"
            call model%eval(input,response,tangent,result,dt=real(6-2*load_step,wp)*tiny(1.0_wp),dsdp=derivative)
            if (result /= cf_ok) error stop "Subnormal forcing branch evaluation"
            call close(response(nmandel,1),0.25_wp)
            call close(tangent(nmandel,nmandel,1),2.0_wp)
            call close(derivative(nmandel,1,1)/epsilon(1.0_wp),1.0_wp/real(9*load_step,wp))
        end do

        ! Unrounded forcing=quantum/4, tau2=dt=quantum: the derivative is 1/16
        ! even though both the forcing and branch stress round to zero.
        law%branch_shear_moduli(0) = quantum
        law%relaxation_times(4) = quantum
        call model%init(law,1,result,initial_input,sensitivity_parameters=selected)
        if (result /= cf_ok) error stop "Underflowed forcing branch initialization"
        call model%eval(input,response,tangent,result,dt=quantum,dsdp=derivative)
        if (result /= cf_ok) error stop "Underflowed forcing branch evaluation"
        call close(response(nmandel,1),0.25_wp)
        call close(tangent(nmandel,nmandel,1),2.0_wp)
        call close(derivative(nmandel,1,1),1.0_wp/16.0_wp)
    end subroutine test_subnormal_branch_reference

    subroutine test_hydrostatic_history()
        !! Uniform loading, hold and reversal excite no branch or shear parameter.
        real(wp), parameter :: dilations(3) = [0.1_wp,0.1_wp,-0.1_wp]
        type(material_t) :: model
        real(wp) :: initial_input(nmandel,1), input(nmandel,1), response(nmandel,1)
        real(wp) :: tangent(nmandel,nmandel,1), derivative(nmandel,nparameters,1)
        real(wp) :: history(nstate,1), history_derivative(nstate,nparameters,1), saved_input(nmandel,1)
        real(wp) :: small_bulk, dilation
        integer :: result, load_step, component

        small_bulk = epsilon(1.0_wp)**2
        initial_input = 0.0_wp
        call model%init(maxwell_multi_t(small_bulk,1.0_wp,[1.0_wp,2.0_wp],[1.0_wp,2.0_wp]), &
            1,result,initial_input,sensitivity_parameters=names)
        if (result /= cf_ok) error stop "Hydrostatic Maxwell multi initialization"
        do load_step = 1,size(dilations)
            dilation = dilations(load_step)
            input = 0.0_wp
            input(1:3,1) = dilation
            call model%eval(input,response,tangent,result,dt=1.0_wp,dsdp=derivative)
            if (result /= cf_ok) error stop "Hydrostatic Maxwell multi evaluation"
            if (.not. all(ieee_is_finite(response)) .or. .not. all(ieee_is_finite(tangent)) .or. &
                .not. all(ieee_is_finite(derivative))) error stop "Nonfinite hydrostatic Maxwell multi response"
            do component = 1,3
                call close(response(component,1)/(small_bulk*dilation),3.0_wp)
                call close(derivative(component,1,1),3.0_wp*dilation)
            end do
            if (maxval(abs(response(4:nmandel,1))) > 0.0_wp .or. &
                maxval(abs(derivative(:,2:nparameters,1))) > 0.0_wp .or. &
                maxval(abs(derivative(4:nmandel,1,1))) > 0.0_wp) then
                error stop "Spurious hydrostatic Maxwell multi shear response"
            end if
            call model%commit(result)
            if (result /= cf_ok) error stop "Hydrostatic Maxwell multi commit"
            call model%export_committed(history,result,previous_input=saved_input,sensitivity_state=history_derivative)
            if (result /= cf_ok) error stop "Hydrostatic Maxwell multi export"
            if (.not. all(ieee_is_finite(history)) .or. .not. all(ieee_is_finite(history_derivative)) .or. &
                .not. all(ieee_is_finite(saved_input))) error stop "Nonfinite hydrostatic Maxwell multi history"
            if (maxval(abs(history)) > 0.0_wp .or. maxval(abs(history_derivative)) > 0.0_wp) then
                error stop "Spurious hydrostatic Maxwell multi branch history"
            end if
            call close(maxval(abs(saved_input-input)),0.0_wp)
        end do
    end subroutine test_hydrostatic_history

    subroutine replay(mat, current_step, out, history)
        type(material_t), intent(inout) :: mat
        integer, intent(in) :: current_step
        real(wp), intent(out) :: out(nmandel,2), history(nstate,2)
        real(wp) :: path(nmandel,2)
        real(wp) :: tstep
        integer :: k, stat

        path(:,1) = [0.02_wp,-0.01_wp,-0.01_wp,0.005_wp,0.0_wp,0.0_wp]
        path(:,2) = [-0.006_wp,0.004_wp,0.002_wp,0.0_wp,0.003_wp,-0.002_wp]
        do k = 1, current_step
            tstep = 0.15_wp+0.1_wp*real(k,wp)
            path(1,:) = path(1,:)+0.003_wp*real(k,wp)
            path(4,:) = path(4,:)-0.001_wp*real(k,wp)
            call mat%eval(path,out,status=stat,dt=tstep)
            if (stat /= cf_ok) error stop "Perturbed Maxwell history evaluation"
            call mat%commit(stat)
            if (stat /= cf_ok) error stop "Perturbed Maxwell history commit"
        end do
        call mat%export_committed(history,stat,previous_input=path)
        if (stat /= cf_ok) error stop "Perturbed Maxwell history export"
    end subroutine replay

end program test_maxwell_multi

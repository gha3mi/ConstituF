! SPDX-License-Identifier: Apache-2.0

program test_cohesive_mixed_mode
    !! Test convention: opening in m, traction in Pa and stiffness in Pa/m.
    !! Interface energy per area is J/m2. Explicit rescaling cases state their factors.
    !!
    !! Mixed opening, bilinear envelope, contact, tangent and history.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
    use constituf, only: wp,material_t,cf_ok,cf_invalid_argument,cf_evaluation_failed, &
                         cf_nondifferentiable,cf_no_trial
    use constituf_cohesive_mixed_mode, only: cohesive_mixed_mode_t
    implicit none (type, external)
    integer,parameter :: nvec=3,nparameters=4
    type(material_t) :: mat,fresh
    real(wp) :: x(3,1),y(3,1),t(3,3,1),xp(3,1),yp(3,1),ym(3,1)
    real(wp) :: h,numerical,error_bound
    real(wp),parameter :: stiffness=1000.0_wp,onset=0.01_wp,failure=0.04_wp,beta=1.5_wp
    integer :: status,i,j

    call mat%init(cohesive_mixed_mode_t(stiffness,onset,failure,beta),1,status)
    call expect(status,cf_ok)
    if (mat%input_size()/=3.or.mat%output_size()/=3.or.mat%state_size()/=1) error stop "Mixed schema"
    if (.not. mat%field_matches("input","opening.jump","vector", &
        "displacement_jump","interface",[3])) error stop "Mixed field semantics"

    x(:,1)=[0.008_wp,0.0_wp,0.0_wp]
    call mat%eval(x,y,t,status)
    call expect(status,cf_ok)
    call close(y(1,1),8.0_wp)
    call close(t(1,1,1),stiffness)
    x(1,1)=0.02_wp
    call mat%eval(x,y,t,status)
    call expect(status,cf_ok)
    call close(y(1,1),stiffness*onset*(failure-x(1,1))/(failure-onset))
    call close(t(1,1,1),-stiffness*onset/(failure-onset))
    call mat%commit(status)
    call expect(status,cf_ok)

    x(1,1)=0.005_wp
    call mat%eval(x,y,t,status)
    call expect(status,cf_ok)
    call close(y(1,1),stiffness*x(1,1)/3.0_wp)
    call close(t(1,1,1),stiffness/3.0_wp)
    call mat%revert()
    x(1,1)=0.03_wp
    call mat%eval(x,y,t,status)
    call expect(status,cf_ok)
    call mat%revert()
    x(:,1)=[-0.004_wp,0.002_wp,0.0_wp]
    call mat%eval(x,y,t,status)
    call expect(status,cf_ok)
    call close(y(1,1),-4.0_wp)
    call close(y(2,1),stiffness*beta**2*x(2,1)/3.0_wp)
    call close(t(1,1,1),stiffness)
    call mat%revert()
    call mat%commit(status)
    call expect(status,cf_no_trial)

    call fresh%init(cohesive_mixed_mode_t(stiffness,onset,failure,beta),1,status)
    call expect(status,cf_ok)
    x(:,1)=[0.015_wp,0.003_wp,0.002_wp]
    call fresh%eval(x,y,t,status)
    call expect(status,cf_ok)
    if (.not. all(ieee_is_finite(t))) error stop "Nonfinite mixed cohesive analytical tangent"
    h=maxval(abs(x))*epsilon(1.0_wp)**(1.0_wp/3.0_wp)
    do j=1,3
        xp=x
        xp(j,1)=xp(j,1)+h
        call fresh%eval(xp,yp,status=status)
        call expect(status,cf_ok)
        xp(j,1)=xp(j,1)-2.0_wp*h
        call fresh%eval(xp,ym,status=status)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(yp)).or..not. all(ieee_is_finite(ym))) then
            error stop "Nonfinite mixed cohesive tangent response"
        end if
        do i=1,3
            numerical=(yp(i,1)-ym(i,1))/(2.0_wp*h)
            error_bound=32.0_wp*epsilon(1.0_wp)*maxval(abs(yp)+abs(ym))/(2.0_wp*h)+ &
                16.0_wp*(h/maxval(abs(x)))**2*maxval(abs(t))
            if (.not. all(ieee_is_finite([numerical,error_bound]))) then
                error stop "Nonfinite mixed cohesive tangent comparison"
            end if
            if (abs(numerical-t(i,j,1))>error_bound) then
                error stop "Mixed cohesive centered tangent mismatch"
            end if
        end do
    end do
    do j=1,3
        do i=1,3
            call close(t(i,j,1),t(j,i,1))
        end do
    end do

    x(:,1)=[onset,0.0_wp,0.0_wp]
    call fresh%eval(x,y,t,status)
    call expect(status,cf_nondifferentiable)
    call fresh%eval(x,y,status=status)
    call expect(status,cf_ok)
    call fresh%revert()

    call fresh%init(cohesive_mixed_mode_t(stiffness,onset,onset,beta),1,status)
    call expect(status,cf_invalid_argument)
    call fresh%init(cohesive_mixed_mode_t(stiffness,onset,failure,beta),1,status)
    call expect(status,cf_ok)
    x(2,1)=ieee_value(0.0_wp,ieee_quiet_nan)
    call fresh%eval(x,y,status=status)
    call expect(status,cf_evaluation_failed)

    call test_exact_mixed_path(1.0_wp)
    call test_exact_mixed_path(sqrt(tiny(1.0_wp))*epsilon(1.0_wp))
    call test_exact_mixed_path(sqrt(huge(1.0_wp))/epsilon(1.0_wp))
    call test_proportional_work()
    call test_history_switches(1.0_wp)
    call test_history_switches(sqrt(tiny(1.0_wp))*epsilon(1.0_wp))
    call test_history_switches(sqrt(huge(1.0_wp))/epsilon(1.0_wp))
    call test_threshold_unloading()
    call test_separated_thresholds()
    call test_small_shear_weight()

contains

    subroutine expect(actual,expected)
        integer,intent(in) :: actual,expected
        if (actual/=expected) then
            print *, "Mixed cohesive status",actual,expected
            error stop "Mixed cohesive status mismatch"
        end if
    end subroutine expect

    subroutine close(actual,expected)
        real(wp),intent(in) :: actual,expected
        if (.not. ieee_is_finite(actual).or..not. ieee_is_finite(expected)) then
            error stop "Nonfinite mixed cohesive comparison"
        end if
        if (abs(actual-expected)>512.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(expected))) then
            print *, "Mixed cohesive values",actual,expected
            error stop "Mixed cohesive reference mismatch"
        end if
    end subroutine close

    subroutine test_exact_mixed_path(length_scale)
        !! A 3-4-5 effective direction with a signed rotation of the shears.
        !! beta-dependent maximum history must survive unloading and contact.
        real(wp),intent(in) :: length_scale
        integer,parameter :: npath=5
        real(wp),parameter :: jump(nvec)=[3.0_wp/250.0_wp,4.0_wp/625.0_wp,-16.0_wp/1875.0_wp]
        real(wp),parameter :: load_traction(nvec)=[4.0_wp,24.0_wp/5.0_wp,-32.0_wp/5.0_wp]
        real(wp),parameter :: load_beta(nvec)=[-256.0_wp/75.0_wp,288.0_wp/125.0_wp,-384.0_wp/125.0_wp]
        real(wp),parameter :: load_failure(nvec)=[200.0_wp/3.0_wp,80.0_wp,-320.0_wp/3.0_wp]
        real(wp),parameter :: load_onset(nvec)=[1600.0_wp/3.0_wp,640.0_wp,-2560.0_wp/3.0_wp]
        real(wp),parameter :: load_tangent(nvec,nvec)=reshape( &
            [280.0_wp/3.0_wp,-288.0_wp,384.0_wp, &
             -288.0_wp,2022.0_wp/5.0_wp,2304.0_wp/5.0_wp, &
             384.0_wp,2304.0_wp/5.0_wp,678.0_wp/5.0_wp],[nvec,nvec])
        type(material_t) :: model
        real(wp) :: input(nvec,1),response(nvec,1),jacobian(nvec,nvec,1)
        real(wp) :: dsdp(nvec,nparameters,1),state(1,1),state_derivative(1,nparameters,1)
        real(wp) :: expected(nvec),expected_derivative(nvec,nparameters),expected_tangent(nvec,nvec)
        real(wp) :: maximum,maximum_beta
        integer :: result,step,component,column

        call model%init(cohesive_mixed_mode_t(stiffness/length_scale,onset*length_scale, &
            failure*length_scale,beta),1,result,sensitivity_parameters=[character(len=16) :: &
                "shear_weight","failure_opening","stiffness","onset_opening"])
        call expect(result,cf_ok)
        do step=1,npath
            input(:,1)=jump*length_scale
            expected=load_traction
            expected_derivative(:,1)=load_beta
            expected_derivative(:,2)=load_failure
            expected_derivative(:,3)=load_traction/stiffness
            expected_derivative(:,4)=load_onset
            expected_tangent=load_tangent
            maximum=0.02_wp
            maximum_beta=16.0_wp/1875.0_wp
            if (step==2.or.step==3) then
                input=0.5_wp*input
                expected=0.5_wp*expected
                expected_derivative=0.5_wp*expected_derivative
                expected_tangent=0.0_wp
                expected_tangent(1,1)=1000.0_wp/3.0_wp
                expected_tangent(2,2)=750.0_wp
                expected_tangent(3,3)=750.0_wp
            else if (step>=4) then
                input=2.5_wp*input
                expected=0.0_wp
                expected_derivative=0.0_wp
                expected_tangent=0.0_wp
                maximum=0.05_wp
                maximum_beta=8.0_wp/375.0_wp
            end if
            if (step==3.or.step==5) then
                input(:,1)=0.5_wp*jump*length_scale
                input(1,1)=-0.006_wp*length_scale
                expected(1)=-6.0_wp
                expected_derivative(1,:)=0.0_wp
                expected_derivative(1,3)=-0.006_wp
                expected_tangent(1,1)=stiffness
            end if
            call model%eval(input,response,jacobian,result,dsdp=dsdp)
            call expect(result,cf_ok)
            do component=1,nvec
                call close(response(component,1),expected(component))
                call close(dsdp(component,1,1),expected_derivative(component,1))
                call close(dsdp(component,2,1)*length_scale,expected_derivative(component,2))
                call close(dsdp(component,3,1)/length_scale,expected_derivative(component,3))
                call close(dsdp(component,4,1)*length_scale,expected_derivative(component,4))
                do column=1,nvec
                    call close(jacobian(component,column,1)*length_scale,expected_tangent(component,column))
                end do
            end do
            call model%commit(result)
            call expect(result,cf_ok)
            call model%export_committed(state,result,sensitivity_state=state_derivative)
            call expect(result,cf_ok)
            call close(state(1,1)/length_scale,maximum)
            call close(state_derivative(1,1,1)/length_scale,maximum_beta)
            if (.not. all(ieee_is_finite(state_derivative))) error stop "Nonfinite mixed cohesive history derivative"
            call close(maxval(abs(state_derivative(1,2:nparameters,1))),0.0_wp)
            if (step==1) then
                ! Identical prescribed jumps have equal beta derivatives of
                ! both max arguments, even though their opening tangent kinks.
                call model%eval(input,response,status=result,dsdp=dsdp)
                call expect(result,cf_ok)
                do component=1,nvec
                    call close(response(component,1),load_traction(component))
                    call close(dsdp(component,1,1),load_beta(component))
                end do
                call model%revert()
                call model%eval(input,response,jacobian,result,dsdp=dsdp)
                call expect(result,cf_nondifferentiable)
                call model%commit(result)
                call expect(result,cf_no_trial)
            end if
        end do
    end subroutine test_exact_mixed_path

    subroutine test_proportional_work()
        !! Traction dot displacement gives exact trapezoidal work on each
        !! straight segment of this proportional bilinear loading path.
        type(material_t) :: model
        real(wp),parameter :: direction(nvec)=[3.0_wp/5.0_wp,8.0_wp/25.0_wp,-32.0_wp/75.0_wp]
        real(wp) :: input(nvec,1),response(nvec,1),peak(nvec),loaded(nvec)
        real(wp) :: work,stored,recovered
        integer :: result

        call model%init(cohesive_mixed_mode_t(stiffness,onset,failure,beta),1,result)
        call expect(result,cf_ok)
        input(:,1)=onset*direction
        call model%eval(input,response,status=result)
        call expect(result,cf_ok)
        peak=response(:,1)
        call close(dot_product(peak,direction),10.0_wp)
        call model%commit(result)
        call expect(result,cf_ok)
        work=0.5_wp*dot_product(peak,input(:,1))
        input(:,1)=0.02_wp*direction
        call model%eval(input,response,status=result)
        call expect(result,cf_ok)
        loaded=response(:,1)
        work=work+0.5_wp*dot_product(peak+loaded,0.01_wp*direction)
        stored=0.5_wp*dot_product(loaded,input(:,1))
        call close(work,2.0_wp/15.0_wp)
        call close(stored,1.0_wp/15.0_wp)
        call close(work-stored,1.0_wp/15.0_wp)
        call model%commit(result)
        call expect(result,cf_ok)
        input(:,1)=0.01_wp*direction
        call model%eval(input,response,status=result)
        call expect(result,cf_ok)
        recovered=0.5_wp*dot_product(loaded+response(:,1),0.01_wp*direction)
        call close(recovered,1.0_wp/20.0_wp)
        call close(0.5_wp*dot_product(response(:,1),input(:,1)),1.0_wp/60.0_wp)
        call model%revert()
        input(:,1)=failure*direction
        call model%eval(input,response,status=result)
        call expect(result,cf_ok)
        work=work+0.5_wp*dot_product(loaded+response(:,1),0.02_wp*direction)
        if (.not. all(ieee_is_finite(response))) error stop "Nonfinite mixed cohesive failure traction"
        call close(maxval(abs(response)),0.0_wp)
        call close(work,1.0_wp/5.0_wp)
    end subroutine test_proportional_work

    subroutine test_history_switches(length_scale)
        real(wp),intent(in) :: length_scale
        type(material_t) :: model
        real(wp) :: input(nvec,1),response(nvec,1),jacobian(nvec,nvec,1)
        real(wp) :: dsdp(nvec,1,1),state(1,1),state_derivative(1,1,1)
        integer :: result

        call model%init(cohesive_mixed_mode_t(stiffness/length_scale,onset*length_scale,failure*length_scale,beta),1,result, &
            sensitivity_parameters=[character(len=12) :: "shear_weight"])
        call expect(result,cf_ok)
        input=0.0_wp
        call model%eval(input,response,jacobian,result,dsdp=dsdp)
        call expect(result,cf_ok)
        if (.not. all(ieee_is_finite(response)).or..not. all(ieee_is_finite(dsdp))) then
            error stop "Nonfinite mixed cohesive origin response"
        end if
        call close(maxval(abs(response)),0.0_wp)
        call close(maxval(abs(dsdp)),0.0_wp)
        call close(jacobian(1,1,1)*length_scale,stiffness)
        input(1,1)=-0.002_wp*length_scale
        call model%eval(input,response,jacobian,result,dsdp=dsdp)
        call expect(result,cf_ok)
        call close(response(1,1),-2.0_wp)
        if (.not. all(ieee_is_finite(dsdp))) error stop "Nonfinite mixed cohesive contact derivative"
        call close(maxval(abs(dsdp)),0.0_wp)
        call model%commit(result)
        call expect(result,cf_ok)
        call model%export_committed(state,result,sensitivity_state=state_derivative)
        call expect(result,cf_ok)
        call close(state(1,1),0.0_wp)
        call close(state_derivative(1,1,1),0.0_wp)

        input(1,1)=0.005_wp*length_scale
        call model%eval(input,response,status=result,dsdp=dsdp)
        call expect(result,cf_ok)
        call model%commit(result)
        call expect(result,cf_ok)
        input(:,1)=[0.0_wp,0.005_wp/beta,0.0_wp]*length_scale
        ! Equal effective separations have different beta derivatives here.
        ! The elastic traction is smooth but its future damage history is not.
        call model%eval(input,response,status=result,dsdp=dsdp)
        call expect(result,cf_nondifferentiable)
        call model%commit(result)
        call expect(result,cf_no_trial)
    end subroutine test_history_switches

    subroutine test_threshold_unloading()
        !! Import a physical accepted normal-opening history; its parameter
        !! derivatives vanish because the prescribed normal path is fixed.
        type(material_t) :: plain,tracked,selected
        real(wp) :: input(nvec,2),response(nvec,2),jacobian(nvec,nvec,2)
        real(wp) :: dsdp(nvec,nparameters,2),state(1,2),state_derivative(1,nparameters,2)
        real(wp) :: single_dsdp(nvec,1,2),single_state_derivative(1,1,2)
        real(wp) :: nearby_failure
        integer :: result,offset

        call plain%init(cohesive_mixed_mode_t(stiffness,onset,failure,beta),2,result)
        call expect(result,cf_ok)
        input=0.0_wp
        input(1,:)=[onset,failure]
        call plain%eval(input,response,status=result)
        call expect(result,cf_ok)
        call plain%commit(result)
        call expect(result,cf_ok)
        call plain%export_committed(state,result)
        call expect(result,cf_ok)
        call tracked%init(cohesive_mixed_mode_t(stiffness,onset,failure,beta),2,result, &
            sensitivity_parameters=[character(len=16) :: &
                "shear_weight","failure_opening","stiffness","onset_opening"])
        call expect(result,cf_ok)
        state_derivative=0.0_wp
        call tracked%import_committed(state,result,sensitivity_state=state_derivative)
        call expect(result,cf_ok)
        input(1,:)=0.5_wp*state(1,:)
        call plain%eval(input,response,jacobian,result)
        call expect(result,cf_ok)
        call close(jacobian(1,1,1),stiffness)
        call close(jacobian(1,1,2),0.0_wp)
        call tracked%eval(input,response,status=result,dsdp=dsdp)
        call expect(result,cf_nondifferentiable)
        call tracked%commit(result)
        call expect(result,cf_no_trial)
        input=0.0_wp
        call tracked%eval(input,response,status=result,dsdp=dsdp)
        call expect(result,cf_ok)
        if (.not. all(ieee_is_finite(response)).or..not. all(ieee_is_finite(dsdp))) then
            error stop "Nonfinite mixed cohesive threshold response"
        end if
        call close(maxval(abs(response)),0.0_wp)
        call close(maxval(abs(dsdp)),0.0_wp)

        call selected%init(cohesive_mixed_mode_t(stiffness,onset,failure,beta),2,result, &
            sensitivity_parameters=[character(len=9) :: "stiffness"])
        call expect(result,cf_ok)
        single_state_derivative=0.0_wp
        call selected%import_committed(state,result,sensitivity_state=single_state_derivative)
        call expect(result,cf_ok)
        input(1,:)=0.5_wp*state(1,:)
        call selected%eval(input,response,jacobian,result,dsdp=single_dsdp)
        call expect(result,cf_ok)
        call close(single_dsdp(1,1,1),0.5_wp*onset)
        call close(single_dsdp(1,1,2),0.0_wp)

        ! Overlapping proximity bands must still distinguish the two exact
        ! endpoints and the parameter derivative that is unique at each one.
        nearby_failure=1.0_wp+16.0_wp*epsilon(1.0_wp)
        call plain%init(cohesive_mixed_mode_t(1.0_wp,1.0_wp,nearby_failure,1.0_wp),2,result)
        call expect(result,cf_ok)
        call tracked%init(cohesive_mixed_mode_t(1.0_wp,1.0_wp,nearby_failure,1.0_wp),2,result, &
            sensitivity_parameters=[character(len=16) :: "onset_opening"])
        call expect(result,cf_ok)
        call selected%init(cohesive_mixed_mode_t(1.0_wp,1.0_wp,nearby_failure,1.0_wp),2,result, &
            sensitivity_parameters=[character(len=16) :: "failure_opening"])
        call expect(result,cf_ok)
        input=0.0_wp
        input(1,:)=1.0_wp
        call tracked%eval(input,response,status=result,dsdp=single_dsdp)
        call expect(result,cf_nondifferentiable)
        call selected%eval(input,response,status=result,dsdp=single_dsdp)
        call expect(result,cf_ok)
        if (.not. all(ieee_is_finite(single_dsdp))) error stop "Nonfinite mixed onset derivative"
        call close(maxval(abs(single_dsdp)),0.0_wp)
        call selected%revert()

        input(1,:)=nearby_failure
        call selected%eval(input,response,status=result,dsdp=single_dsdp)
        call expect(result,cf_nondifferentiable)
        call tracked%eval(input,response,status=result,dsdp=single_dsdp)
        call expect(result,cf_ok)
        if (.not. all(ieee_is_finite(response)).or..not. all(ieee_is_finite(single_dsdp))) then
            error stop "Nonfinite mixed endpoint response"
        end if
        call close(maxval(abs(response)),0.0_wp)
        call close(maxval(abs(single_dsdp)),0.0_wp)
        call tracked%revert()

        ! First arrival is a kink; after failure is accepted, positive-normal
        ! opening has zero tangent even inside the nearby onset switch band.
        call plain%eval(input,response,jacobian,result)
        call expect(result,cf_nondifferentiable)
        call plain%eval(input,response,status=result)
        call expect(result,cf_ok)
        call plain%commit(result)
        call expect(result,cf_ok)
        do offset=-1,1
            input(1,:)=nearby_failure
            if (offset/=0) input(1,2)=nearest(nearby_failure,real(offset,wp))
            call plain%eval(input,response,jacobian,result)
            call expect(result,cf_ok)
            if (.not. all(ieee_is_finite(response)).or..not. all(ieee_is_finite(jacobian))) then
                error stop "Nonfinite accepted mixed failure response"
            end if
            call close(maxval(abs(response)),0.0_wp)
            call close(maxval(abs(jacobian)),0.0_wp)
        end do
        call plain%revert()

        ! Unilateral contact remains nonsmooth at zero normal jump after failure.
        input(1,2)=0.0_wp
        call plain%eval(input,response,jacobian,result)
        call expect(result,cf_nondifferentiable)
        if (plain%last_failure_point()/=2) error stop "Missing failed-interface contact kink"
        call plain%commit(result)
        call expect(result,cf_no_trial)
        call plain%eval(input,response,status=result)
        call expect(result,cf_ok)
        if (.not. all(ieee_is_finite(response))) error stop "Nonfinite failed-interface contact traction"
        call close(maxval(abs(response)),0.0_wp)
    end subroutine test_threshold_unloading

    subroutine test_separated_thresholds()
        type(material_t) :: model
        real(wp) :: input(nvec,1),response(nvec,1),jacobian(nvec,nvec,1)
        real(wp) :: small_onset,normalized
        integer :: result

        small_onset=epsilon(1.0_wp)**2
        call model%init(cohesive_mixed_mode_t(1.0_wp/small_onset,small_onset,1.0_wp,1.0_wp),1,result)
        call expect(result,cf_ok)
        input=0.0_wp
        input(1,1)=2.0_wp*small_onset
        call model%eval(input,response,jacobian,result)
        call expect(result,cf_ok)
        call close(response(1,1),1.0_wp)
        call close(jacobian(1,1,1),-1.0_wp/(1.0_wp-small_onset))

        ! The near-normal transverse projector remains O(eps), even if the
        ! normal unit-vector component rounds to one. Its contribution is
        ! amplified by K, with eps*T11=1/8-(67/64)*eps+O(eps**2).
        input(2,1)=small_onset*sqrt(epsilon(1.0_wp))
        call model%eval(input,response,jacobian,result)
        call expect(result,cf_ok)
        normalized=epsilon(1.0_wp)*jacobian(1,1,1)
        if (.not. ieee_is_finite(normalized).or. &
            abs(normalized-0.125_wp)>16.0_wp*epsilon(1.0_wp)) then
            error stop "Mixed cohesive near-normal projector cancellation"
        end if
    end subroutine test_separated_thresholds

    subroutine test_small_shear_weight()
        !! beta**2 underflows, but K*beta**2 and every requested output remain
        !! representable. The large prescribed shear jump keeps beta*u_s=1.
        type(material_t) :: model
        real(wp) :: input(nvec,1),response(nvec,1),jacobian(nvec,nvec,1)
        real(wp) :: dsdp(nvec,2,1),state(1,1),state_derivative(1,2,1)
        real(wp) :: weight
        integer :: result

        weight=sqrt(tiny(1.0_wp))*epsilon(1.0_wp)
        call model%init(cohesive_mixed_mode_t(1.0_wp/weight,2.0_wp,4.0_wp,weight),1,result, &
            sensitivity_parameters=[character(len=12) :: "shear_weight","stiffness"])
        call expect(result,cf_ok)
        input(:,1)=[0.5_wp,1.0_wp/weight,0.0_wp]
        call model%eval(input,response,jacobian,result,dsdp=dsdp)
        call expect(result,cf_ok)
        call close(response(1,1)*weight,0.5_wp)
        call close(response(2,1),1.0_wp)
        call close(jacobian(1,1,1)*weight,1.0_wp)
        call close(jacobian(2,2,1)/weight,1.0_wp)
        call close(jacobian(3,3,1)/weight,1.0_wp)
        call close(dsdp(1,1,1),0.0_wp)
        call close(dsdp(2,1,1)*weight,2.0_wp)
        call close(dsdp(1,2,1),0.5_wp)
        call close(dsdp(2,2,1)/weight,1.0_wp)
        call model%commit(result)
        call expect(result,cf_ok)
        call model%export_committed(state,result,sensitivity_state=state_derivative)
        call expect(result,cf_ok)
        call close(state(1,1),sqrt(1.25_wp))
        call close(state_derivative(1,1,1)*weight,2.0_wp/sqrt(5.0_wp))
        call close(state_derivative(1,2,1),0.0_wp)
    end subroutine test_small_shear_weight

end program test_cohesive_mixed_mode

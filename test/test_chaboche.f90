! SPDX-License-Identifier: Apache-2.0

program test_chaboche
    !! Test convention: stress and stress-like material parameters in Pa.
    !! Strain and deformation gradients are dimensionless. Explicit rescaling
    !! cases state their factors.
    !!
    !! One-backstress nonlinear cyclic return, exact tangent, J2 limit.
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_is_finite
    use constituf, only: wp,material_t,verify_tangent,cf_ok,cf_invalid_argument,cf_evaluation_failed,cf_no_trial
    use constituf_chaboche, only: chaboche_t
    use constituf_j2_plasticity, only: j2_plasticity_t
    implicit none (type, external)

    integer, parameter :: nmandel=6,nstate=13,nparameters=6
    type(material_t) :: mat,j2,virgin
    real(wp) :: x(nmandel,2),y(nmandel,2),t(nmandel,nmandel,2),xp(nmandel,2),yp(nmandel,2),baseline(nmandel,2)
    real(wp) :: history(nstate,2),trial(nmandel,2)
    real(wp) :: h,error,solve_scale,difference_tolerance
    real(wp),parameter :: reference_tolerance=512.0_wp*epsilon(1.0_wp)
    real(wp),parameter :: bulk=1000.0_wp,shear=400.0_wp,sy=0.2_wp,hard=40.0_wp
    integer :: status,i,ipoint

    call test_small_increment()
    call test_exact_recovery(1.0_wp)
    call test_exact_recovery(epsilon(1.0_wp)**2)
    call test_small_plastic_increment()
    call test_cyclic_shear(1.0_wp)
    call test_cyclic_shear(epsilon(1.0_wp)**2)
    call test_path_sensitivities()

    x(:,1)=[0.01_wp,-0.005_wp,-0.005_wp,0.02_wp,0.005_wp,0.0_wp]
    x(:,2)=[0.006_wp,-0.003_wp,-0.003_wp,0.01_wp,0.0_wp,0.002_wp]
    call mat%init(chaboche_t(bulk,shear,sy,hard,0.0_wp,0.0_wp),2,status)
    call expect(status,cf_ok)
    if (mat%state_size()/=nstate) error stop "Chaboche history schema"
    call j2%init(j2_plasticity_t(bulk,shear,sy,hard),2,status)
    call expect(status,cf_ok)
    call mat%eval(x,y,t,status)
    call expect(status,cf_ok)
    call j2%eval(x,yp,status=status)
    call expect(status,cf_ok)
    solve_scale=max(sy,sqrt(1.5_wp*maxval(sum((2.0_wp*shear*x)**2,dim=1))))
    do ipoint=1,2
        do i=1,nmandel
            call close(y(i,ipoint)/solve_scale,yp(i,ipoint)/solve_scale,128.0_wp*epsilon(1.0_wp))
        end do
    end do
    call mat%commit(status)
    call expect(status,cf_ok)
    call j2%commit(status)
    call expect(status,cf_ok)
    call mat%export_committed(history,status)
    call expect(status,cf_ok)
    x=-0.65_wp*x
    trial=2.0_wp*shear*(x-history(1:nmandel,:))
    solve_scale=max(sy+hard*maxval(history(nstate,:)),sqrt(1.5_wp*maxval(sum(trial**2,dim=1))))
    call mat%eval(x,y,t,status)
    call expect(status,cf_ok)
    call j2%eval(x,yp,status=status)
    call expect(status,cf_ok)
    do ipoint=1,2
        do i=1,nmandel
            call close(y(i,ipoint)/solve_scale,yp(i,ipoint)/solve_scale,128.0_wp*epsilon(1.0_wp))
        end do
    end do

    ! Finite dynamic recovery plus backstress: accepted/reverted paths.
    call mat%init(chaboche_t(bulk,shear,sy,hard,150.0_wp,5.0_wp),2,status)
    call expect(status,cf_ok)
    x(:,1)=[0.01_wp,-0.005_wp,-0.005_wp,0.02_wp,0.005_wp,0.0_wp]
    x(:,2)=[0.006_wp,-0.003_wp,-0.003_wp,0.01_wp,0.0_wp,0.002_wp]
    call mat%eval(x,y,t,status)
    call expect(status,cf_ok)
    call mat%commit(status)
    call expect(status,cf_ok)
    call mat%export_committed(history,status)
    call expect(status,cf_ok)
    x=-0.8_wp*x
    trial=2.0_wp*shear*(x-history(1:nmandel,:))-history(nmandel+1:2*nmandel,:)
    solve_scale=max(sy+hard*maxval(history(nstate,:)),sqrt(1.5_wp*maxval(sum(trial**2,dim=1))))
    call mat%eval(x,y,t,status)
    call expect(status,cf_ok)
    baseline=y
    call mat%revert()
    call mat%eval(x,y,t,status)
    call expect(status,cf_ok)
    if (.not. all(ieee_is_finite(y)).or..not. all(ieee_is_finite(baseline))) error stop "Nonfinite Chaboche revert"
    if (maxval(abs(y-baseline))>128.0_wp*epsilon(1.0_wp)*solve_scale) error stop "Chaboche revert"

    h=maxval(abs(x))*epsilon(1.0_wp)**(1.0_wp/3.0_wp)
    ! The 32*eps relative return residual is amplified by the difference quotient.
    difference_tolerance=96.0_wp*epsilon(1.0_wp)*solve_scale/h
    call verify_tangent(mat,x,error,status,step=[(h,i=1,nmandel)],absolute_tolerance=difference_tolerance)
    call expect(status,cf_ok)
    if (.not. ieee_is_finite(error).or.error>1.0_wp) error stop "Chaboche consistent tangent"
    call virgin%init(chaboche_t(bulk,shear,sy,hard,150.0_wp,5.0_wp),2,status)
    call expect(status,cf_ok)
    call virgin%eval(x,yp,status=status)
    call expect(status,cf_ok)
    if (.not. all(ieee_is_finite(yp))) error stop "Nonfinite Chaboche virgin response"
    if (maxval(abs(y-yp))<0.01_wp) error stop "Chaboche backstress had no cyclic effect"
    call mat%revert()
    call mat%commit(status)
    call expect(status,cf_no_trial)

    call virgin%init(chaboche_t(bulk,shear,sy,hard,150.0_wp,-1.0_wp),2,status)
    call expect(status,cf_invalid_argument)
    call virgin%init(chaboche_t(bulk,shear,sy,hard,150.0_wp,5.0_wp),2,status)
    call expect(status,cf_ok)
    xp=x
    xp(3,2)=ieee_value(0.0_wp,ieee_quiet_nan)
    call virgin%eval(xp,yp,status=status)
    call expect(status,cf_evaluation_failed)

contains

    subroutine test_exact_recovery(modulus_scale)
        real(wp),intent(in) :: modulus_scale
        type(material_t) :: material
        real(wp),parameter :: direction(nmandel)=[1.0_wp,-0.5_wp,-0.5_wp,0.0_wp,0.0_wp,0.0_wp]
        real(wp) :: strain(nmandel,1),response(nmandel,1),jacobian(nmandel,nmandel,1)
        real(wp) :: derivative(nmandel,2,1),state(nstate,1),state_derivative(nstate,2,1)
        real(wp) :: expected_stress
        integer :: result,load_step,component

        ! The backward-Euler residual has the exact root r=1/20 and 1+gamma*r=3/2.
        call material%init(chaboche_t(100.0_wp*modulus_scale,60.0_wp*modulus_scale, &
            3.0_wp*modulus_scale,15.0_wp*modulus_scale,30.0_wp*modulus_scale,10.0_wp),1,result, &
            sensitivity_parameters=[character(len=12) :: "yield_stress", "recovery"])
        call expect(result,cf_ok)
        strain(:,1)=(11.0_wp/144.0_wp)*direction
        do load_step=1,2
            expected_stress=19.0_wp/6.0_wp
            if (load_step==2) then
                ! Elastic unloading retains both plastic strain and nonlinear backstress derivatives.
                strain=0.5_wp*strain
                expected_stress=-17.0_wp/12.0_wp
            end if
            call material%eval(strain,response,jacobian,result,dsdp=derivative)
            call expect(result,cf_ok)
            do component=1,nmandel
                call close(response(component,1)/modulus_scale,expected_stress*direction(component),reference_tolerance)
                call close(derivative(component,1,1),(72.0_wp/125.0_wp)*direction(component),reference_tolerance)
                call close(derivative(component,2,1)/modulus_scale, &
                    -(12.0_wp/625.0_wp)*direction(component),reference_tolerance)
            end do
            if (load_step==1) then
                call close(jacobian(1,1,1)/modulus_scale,2772.0_wp/25.0_wp,reference_tolerance)
                call close(jacobian(1,2,1)/modulus_scale,2364.0_wp/25.0_wp,reference_tolerance)
                call close(jacobian(4,4,1)/modulus_scale,456.0_wp/11.0_wp,reference_tolerance)
            else
                call close(jacobian(1,1,1)/modulus_scale,180.0_wp,reference_tolerance)
                call close(jacobian(1,2,1)/modulus_scale,60.0_wp,reference_tolerance)
                call close(jacobian(4,4,1)/modulus_scale,120.0_wp,reference_tolerance)
            end if
            call material%commit(result)
            call expect(result,cf_ok)
            call material%export_committed(state,result,sensitivity_state=state_derivative)
            call expect(result,cf_ok)
            do component=1,nmandel
                call close(state(component,1),direction(component)/20.0_wp,reference_tolerance)
                call close(state(nmandel+component,1)/modulus_scale,(2.0_wp/3.0_wp)*direction(component),reference_tolerance)
                call close(state_derivative(component,1,1)*modulus_scale, &
                    -(3.0_wp/625.0_wp)*direction(component),reference_tolerance)
                call close(state_derivative(component,2,1),direction(component)/6250.0_wp,reference_tolerance)
                call close(state_derivative(nmandel+component,1,1), &
                    -(16.0_wp/375.0_wp)*direction(component),reference_tolerance)
                call close(state_derivative(nmandel+component,2,1)/modulus_scale, &
                    -(13.0_wp/625.0_wp)*direction(component),reference_tolerance)
            end do
            call close(state(nstate,1),1.0_wp/20.0_wp,reference_tolerance)
            call close(state_derivative(nstate,1,1)*modulus_scale,-3.0_wp/625.0_wp,reference_tolerance)
            call close(state_derivative(nstate,2,1),1.0_wp/6250.0_wp,reference_tolerance)
        end do
    end subroutine test_exact_recovery

    subroutine test_small_plastic_increment()
        type(material_t) :: material
        real(wp),parameter :: strain_scale=epsilon(1.0_wp)**2
        real(wp),parameter :: direction(nmandel)=[1.0_wp,-0.5_wp,-0.5_wp,0.0_wp,0.0_wp,0.0_wp]
        real(wp) :: strain(nmandel,1),response(nmandel,1),jacobian(nmandel,nmandel,1),state(nstate,1)
        integer :: result,component

        ! The exact plastic multiplier is 7*strain_scale/325. A fixed strain-width stopping
        ! criterion would accept the first midpoint before satisfying the yield condition.
        call material%init(chaboche_t(100.0_wp,60.0_wp,3.0_wp*strain_scale,15.0_wp,0.0_wp,0.0_wp),1,result)
        call expect(result,cf_ok)
        strain(:,1)=0.04_wp*strain_scale*direction
        call material%eval(strain,response,jacobian,result)
        call expect(result,cf_ok)
        call material%commit(result)
        call expect(result,cf_ok)
        call material%export_committed(state,result)
        call expect(result,cf_ok)
        do component=1,nmandel
            call close(response(component,1)/strain_scale,(144.0_wp/65.0_wp)*direction(component),reference_tolerance)
            call close(state(component,1)/strain_scale,(7.0_wp/325.0_wp)*direction(component),reference_tolerance)
            call close(state(nmandel+component,1),0.0_wp,reference_tolerance)
        end do
        call close(state(nstate,1)/strain_scale,7.0_wp/325.0_wp,reference_tolerance)
    end subroutine test_small_plastic_increment

    subroutine test_small_increment()
        type(material_t) :: model
        real(wp) :: input(nmandel,1),stress(nmandel,1),state(nstate,1),trial_stress(nmandel)
        real(wp) :: qtrial,increment,denom
        integer :: code

        ! With zero recovery, the small plastic increment has a linear closed form.
        input=0.0_wp
        input(1,1)=(250.0e6_wp+0.25e6_wp)/(2.0_wp*80.0e9_wp)
        trial_stress=input(:,1)
        trial_stress(1:3)=trial_stress(1:3)-sum(input(1:3,1))/3.0_wp
        trial_stress=2.0_wp*80.0e9_wp*trial_stress
        qtrial=sqrt(1.5_wp*dot_product(trial_stress,trial_stress))
        denom=3.0_wp*80.0e9_wp+10.0e9_wp+1.0e9_wp
        increment=(qtrial-250.0e6_wp)/denom
        call model%init(chaboche_t(160.0e9_wp,80.0e9_wp,250.0e6_wp,1.0e9_wp,10.0e9_wp,0.0_wp),1,code)
        call expect(code,cf_ok)
        call model%eval(input,stress,status=code)
        call expect(code,cf_ok)
        call model%commit(code)
        call expect(code,cf_ok)
        call model%export_committed(state,code)
        call expect(code,cf_ok)
        ! Normalize by the increment: an eps-sized absolute strain bound would hide this error.
        call close(state(nstate,1)/increment,1.0_wp, &
            64.0_wp*epsilon(1.0_wp)*qtrial/(denom*increment))
    end subroutine test_small_increment

    subroutine test_cyclic_shear(modulus_scale)
        real(wp),intent(in) :: modulus_scale
        type(material_t) :: model
        real(wp) :: input(nmandel,1),stress(nmandel,1),tangent(nmandel,nmandel,1),gradient(nmandel,1,1)
        real(wp) :: state(nstate,1),state_gradient(nstate,1,1)
        real(wp) :: expected_state(nstate),expected_gradient(nstate),state_units(nstate),gradient_units(nstate)
        real(wp) :: direction,a,expected_tangent,cyclic_tolerance
        integer :: code,step,component

        ! Finite recovery gives exact multipliers 1/100 and 1/50 for this shear reversal.
        ! Scaling all stress parameters leaves plastic strain and recovery unchanged.
        a=sqrt(1.5_wp)
        cyclic_tolerance=256.0_wp*epsilon(1.0_wp)
        call model%init(chaboche_t(100.0_wp*modulus_scale,60.0_wp*modulus_scale,3.0_wp*modulus_scale, &
            0.0_wp,120.0_wp*modulus_scale,100.0_wp),1,code, &
            sensitivity_parameters=[character(len=12) :: "yield_stress"])
        call expect(code,cf_ok)
        state_units=100.0_wp
        state_units(nmandel+1:2*nmandel)=1.0_wp/modulus_scale
        gradient_units=210.0_wp*modulus_scale
        gradient_units(nmandel+1:2*nmandel)=7.0_wp*a
        do step=1,2
            direction=real(3-2*step,wp)
            input=0.0_wp
            input(nmandel,1)=direction*9.0_wp/(200.0_wp*a)
            call model%eval(input,stress,tangent,code,dsdp=gradient)
            call expect(code,cf_ok)
            do component=1,nmandel-1
                call close(stress(component,1)/modulus_scale,0.0_wp,cyclic_tolerance)
                call close(gradient(component,1,1),0.0_wp,cyclic_tolerance)
            end do
            call close(stress(nmandel,1)/modulus_scale,direction*18.0_wp/(5.0_wp*a),cyclic_tolerance)
            call close(gradient(nmandel,1,1),direction*6.0_wp/(7.0_wp*a),cyclic_tolerance)
            expected_tangent=120.0_wp/7.0_wp
            if (step==2) expected_tangent=12.0_wp
            call close(tangent(nmandel,nmandel,1)/modulus_scale,expected_tangent,cyclic_tolerance)
            call model%commit(code)
            call expect(code,cf_ok)
            call model%export_committed(state,code,sensitivity_state=state_gradient)
            call expect(code,cf_ok)
            expected_state=0.0_wp
            expected_state(nmandel)=direction*a/100.0_wp
            expected_state(2*nmandel)=direction*3.0_wp*modulus_scale/(5.0_wp*a)
            expected_state(nstate)=real(2*step-1,wp)/100.0_wp
            expected_gradient=0.0_wp
            expected_gradient(nmandel)=-direction*a/(210.0_wp*modulus_scale)
            expected_gradient(2*nmandel)=-direction/(7.0_wp*a)
            expected_gradient(nstate)=-real(2*step-1,wp)/(210.0_wp*modulus_scale)
            do component=1,nstate
                call close(state(component,1)*state_units(component), &
                    expected_state(component)*state_units(component),cyclic_tolerance)
                call close(state_gradient(component,1,1)*gradient_units(component), &
                    expected_gradient(component)*gradient_units(component),cyclic_tolerance)
            end do
        end do
    end subroutine test_cyclic_shear

    subroutine test_path_sensitivities()
        integer,parameter :: nsteps=4
        character(len=18),parameter :: names(nparameters)=[character(len=18) :: "bulk_modulus", &
            "shear_modulus","yield_stress","hardening_modulus","kinematic_modulus","recovery"]
        real(wp),parameter :: parameters(nparameters)=[100.0_wp,60.0_wp,3.0_wp,15.0_wp,120.0_wp,100.0_wp]
        type(material_t) :: tracked
        type(material_t) :: plus(nparameters,2),minus(nparameters,2)
        real(wp) :: path(nmandel,nsteps),input(nmandel,1),stress(nmandel,1),gradient(nmandel,nparameters,1)
        real(wp) :: state(nstate,1),state_gradient(nstate,nparameters,1)
        real(wp) :: pplus(nparameters),pminus(nparameters),steps(nparameters,2),trial_stress(nmandel)
        real(wp) :: yp(nmandel,1),ym(nmandel,1),hp(nstate,1),hm(nstate,1)
        real(wp) :: stress_difference(nmandel,2),state_difference(nstate,2),stress_tolerance(nmandel)
        real(wp) :: state_tolerance(nstate),state_roundoff(nstate)
        real(wp) :: accumulated_scale,solve_scale,response_roundoff
        integer :: code,step,parameter,refinement

        ! Plastic load, elastic unloading, nonproportional reversal, then reloading.
        ! A nonzero trace also exercises the bulk-modulus derivative.
        path(:,1)=[0.03_wp,-0.005_wp,-0.01_wp,0.045_wp,0.01_wp,0.0_wp]
        path(:,2)=0.95_wp*path(:,1)
        path(:,3)=-0.6_wp*path(:,1)
        path(nmandel,3)=0.004_wp
        path(:,4)=[0.028_wp,-0.002_wp,-0.008_wp,0.036_wp,-0.014_wp,0.007_wp]
        call tracked%init(chaboche_t(parameters(1),parameters(2),parameters(3),parameters(4), &
            parameters(5),parameters(6)),1,code,sensitivity_parameters=names)
        call expect(code,cf_ok)
        do parameter=1,nparameters
            ! Balance the O(32*eps) return residual against the O(h^2) secant error.
            steps(parameter,1)=(32.0_wp*epsilon(1.0_wp))**(1.0_wp/3.0_wp)*parameters(parameter)
            steps(parameter,2)=0.5_wp*steps(parameter,1)
            do refinement=1,2
                pplus=parameters
                pminus=parameters
                pplus(parameter)=pplus(parameter)+steps(parameter,refinement)
                pminus(parameter)=pminus(parameter)-steps(parameter,refinement)
                call plus(parameter,refinement)%init(chaboche_t(pplus(1),pplus(2),pplus(3),pplus(4), &
                    pplus(5),pplus(6)),1,code)
                call expect(code,cf_ok)
                call minus(parameter,refinement)%init(chaboche_t(pminus(1),pminus(2),pminus(3),pminus(4), &
                    pminus(5),pminus(6)),1,code)
                call expect(code,cf_ok)
            end do
        end do
        state=0.0_wp
        accumulated_scale=0.0_wp
        do step=1,nsteps
            input(:,1)=path(:,step)
            trial_stress=input(:,1)
            trial_stress(1:3)=trial_stress(1:3)-sum(input(1:3,1))/3.0_wp
            trial_stress=2.0_wp*parameters(2)*(trial_stress-state(1:nmandel,1))-state(nmandel+1:2*nmandel,1)
            solve_scale=max(parameters(3)+parameters(4)*state(nstate,1), &
                sqrt(1.5_wp*dot_product(trial_stress,trial_stress)))
            accumulated_scale=accumulated_scale+solve_scale
            call tracked%eval(input,stress,status=code,dsdp=gradient)
            call expect(code,cf_ok)
            call tracked%commit(code)
            call expect(code,cf_ok)
            call tracked%export_committed(state,code,sensitivity_state=state_gradient)
            call expect(code,cf_ok)
            if (.not. all(ieee_is_finite(gradient)).or..not. all(ieee_is_finite(state_gradient))) then
                error stop "Nonfinite Chaboche path derivatives"
            end if
            do parameter=1,nparameters
                do refinement=1,2
                    call plus(parameter,refinement)%eval(input,yp,status=code)
                    call expect(code,cf_ok)
                    call minus(parameter,refinement)%eval(input,ym,status=code)
                    call expect(code,cf_ok)
                    call plus(parameter,refinement)%commit(code)
                    call expect(code,cf_ok)
                    call minus(parameter,refinement)%commit(code)
                    call expect(code,cf_ok)
                    call plus(parameter,refinement)%export_committed(hp,code)
                    call expect(code,cf_ok)
                    call minus(parameter,refinement)%export_committed(hm,code)
                    call expect(code,cf_ok)
                    stress_difference(:,refinement)=(yp(:,1)-ym(:,1))/(2.0_wp*steps(parameter,refinement))
                    state_difference(:,refinement)=(hp(:,1)-hm(:,1))/(2.0_wp*steps(parameter,refinement))
                end do
                if (.not. all(ieee_is_finite(stress_difference)).or..not. all(ieee_is_finite(state_difference))) then
                    error stop "Nonfinite Chaboche history secants"
                end if
                ! Accumulated consistency residuals contribute O(eps/h) to the accepted-history secants.
                response_roundoff=96.0_wp*epsilon(1.0_wp)*accumulated_scale/steps(parameter,2)
                state_roundoff=response_roundoff/(2.0_wp*parameters(2))
                state_roundoff(nmandel+1:2*nmandel)=parameters(5)*response_roundoff/(3.0_wp*parameters(2))
                stress_tolerance=abs(stress_difference(:,2)-stress_difference(:,1))/3.0_wp+response_roundoff
                state_tolerance=abs(state_difference(:,2)-state_difference(:,1))/3.0_wp+state_roundoff
                if (any(abs(stress_difference(:,2)-gradient(:,parameter,1))>stress_tolerance)) then
                    print *, "Chaboche response history derivative",step,parameter
                    error stop "Chaboche parameter derivative mismatch"
                end if
                if (any(abs(state_difference(:,2)-state_gradient(:,parameter,1))>state_tolerance)) then
                    print *, "Chaboche internal history derivative",step,parameter
                    error stop "Chaboche state derivative mismatch"
                end if
            end do
        end do
    end subroutine test_path_sensitivities

    subroutine expect(actual,expected)
        integer,intent(in) :: actual,expected
        if (actual/=expected) then
            print *, "Chaboche status",actual,expected
            error stop "Chaboche status mismatch"
        end if
    end subroutine expect

    subroutine close(actual,expected,relative)
        real(wp),intent(in) :: actual,expected,relative
        if (.not. all(ieee_is_finite([actual,expected]))) error stop "Nonfinite Chaboche comparison"
        if (abs(actual-expected)>relative*max(1.0_wp,abs(expected))) then
            print *, "Chaboche values",actual,expected
            error stop "Chaboche return/tangent mismatch"
        end if
    end subroutine close

end program test_chaboche

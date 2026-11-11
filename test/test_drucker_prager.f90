! SPDX-License-Identifier: Apache-2.0

program test_drucker_prager
    !! Test convention: stress and stress-like material parameters in Pa.
    !! Strain and deformation gradients are dimensionless. Explicit rescaling
    !! cases state their factors.
    !!
    !! Smooth/apex returns, finite-difference tangents, J2 reduction.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
    use constituf, only: wp,material_t,verify_tangent,cf_ok,cf_invalid_argument,cf_evaluation_failed,cf_no_trial, &
        cf_nondifferentiable
    use constituf_drucker_prager, only: drucker_prager_t
    use constituf_j2_plasticity, only: j2_plasticity_t
    implicit none (type, external)

    integer, parameter :: nmandel=6,nstate=7
    type(material_t) :: mat,j2,plain
    real(wp) :: x(nmandel,2),y(nmandel,2),t(nmandel,nmandel,2),xp(nmandel,2),yp(nmandel,2)
    real(wp) :: h,error,step,mean_final,stress_scale
    real(wp),parameter :: reference_tolerance=512.0_wp*epsilon(1.0_wp)
    real(wp),parameter :: bulk=1000.0_wp,shear=400.0_wp,sy=0.2_wp,beta=0.3_wp,hard=50.0_wp
    integer :: status,i,ipoint

    x(:,1)=[0.001_wp,-0.0005_wp,-0.0005_wp,0.02_wp,0.002_wp,0.0_wp]
    x(:,2)=[0.02_wp,0.02_wp,0.02_wp,0.0_wp,0.0_wp,0.0_wp]
    call mat%init(drucker_prager_t(bulk,shear,sy,beta,hard),2,status)
    call expect(status,cf_ok)
    if (mat%state_size()/=nstate.or.mat%input_size()/=nmandel) error stop "DP schema"
    call mat%eval(x,y,t,status)
    call expect(status,cf_ok)
    step=(beta*bulk*sum(x(1:3,2))-sy)/(bulk*beta**2+hard)
    mean_final=bulk*sum(x(1:3,2))-bulk*beta*step
    do i=1,3
        call close(y(i,2),mean_final,reference_tolerance)
    end do
    do i=4,nmandel
        call close(y(i,2),0.0_wp,reference_tolerance)
    end do
    if (.not. ieee_is_finite(y(4,1)).or.y(4,1)<=0.0_wp) error stop "DP smooth-cone shear sign"

    call plain%init(drucker_prager_t(bulk,shear,sy,beta,hard),2,status)
    call expect(status,cf_ok)
    h=maxval(abs(x))*epsilon(1.0_wp)**(1.0_wp/3.0_wp)
    stress_scale=max(sy,3.0_wp*shear*maxval(sqrt(sum(x*x,dim=1))),bulk*maxval(abs(sum(x(1:3,:),dim=1))))
    call verify_tangent(plain,x,error,status,step=[(h,i=1,nmandel)], &
        absolute_tolerance=64.0_wp*epsilon(1.0_wp)*stress_scale/h)
    call expect(status,cf_ok)
    if (.not. ieee_is_finite(error).or.error>1.0_wp) error stop "DP consistent tangent"

    ! beta=0 reproduces the original J2 radial return including hardening.
    call mat%init(drucker_prager_t(bulk,shear,sy,0.0_wp,hard),2,status)
    call expect(status,cf_ok)
    call j2%init(j2_plasticity_t(bulk,shear,sy,hard),2,status)
    call expect(status,cf_ok)
    call mat%eval(x,y,t,status)
    call expect(status,cf_ok)
    call j2%eval(x,yp,status=status)
    call expect(status,cf_ok)
    do ipoint=1,2
        do i=1,nmandel
            call close(y(i,ipoint),yp(i,ipoint),reference_tolerance)
        end do
    end do
    call mat%commit(status)
    call expect(status,cf_ok)
    call j2%commit(status)
    call expect(status,cf_ok)
    x=-0.5_wp*x
    call mat%eval(x,y,status=status)
    call expect(status,cf_ok)
    call j2%eval(x,yp,status=status)
    call expect(status,cf_ok)
    do ipoint=1,2
        do i=1,nmandel
            call close(y(i,ipoint),yp(i,ipoint),reference_tolerance)
        end do
    end do
    call mat%revert()
    call mat%commit(status)
    call expect(status,cf_no_trial)

    call plain%init(drucker_prager_t(bulk,shear,sy,-beta,hard),2,status)
    call expect(status,cf_invalid_argument)
    call plain%init(drucker_prager_t(bulk,shear,sy,beta,hard),2,status)
    call expect(status,cf_ok)
    xp=x
    xp(2,2)=ieee_value(0.0_wp,ieee_quiet_nan)
    call plain%eval(xp,yp,status=status)
    call expect(status,cf_evaluation_failed)

    call test_exact_branches(1.0_wp)
    call test_exact_branches(epsilon(1.0_wp)**2)
    call test_pressure_limits()
    call test_branch_switches()
    call test_path_sensitivities()

contains

    subroutine test_exact_branches(modulus_scale)
        real(wp),intent(in) :: modulus_scale
        type(material_t) :: material
        real(wp) :: strain(nmandel,2),response(nmandel,2),jacobian(nmandel,nmandel,2)
        real(wp) :: derivative(nmandel,2,2),state(nstate,2),state_derivative(nstate,2,2)
        real(wp) :: expected(nmandel,2),expected_state(nstate,2),expected_derivative(nmandel,2,2)
        real(wp) :: expected_state_derivative(nstate,2,2)
        integer :: result,load_step,component,point

        ! Point 1 returns to the smooth cone with delta_gamma=1/100 and q_new=6/25.
        ! Point 2 returns to the apex with delta_gamma=1/16 and p_new=105/8.
        call material%init(drucker_prager_t(100.0_wp*modulus_scale,60.0_wp*modulus_scale, &
            3.0_wp*modulus_scale,0.3_wp,15.0_wp*modulus_scale),2,result, &
            sensitivity_parameters=[character(len=14) :: "yield_stress", "pressure_slope"])
        call expect(result,cf_ok)
        strain=0.0_wp
        strain(1:3,1)=[67.0_wp/1500.0_wp,83.0_wp/3000.0_wp,83.0_wp/3000.0_wp]
        strain(1:3,2)=0.05_wp
        expected=0.0_wp
        expected(1:3,1)=[493.0_wp/50.0_wp,481.0_wp/50.0_wp,481.0_wp/50.0_wp]
        expected(1:3,2)=105.0_wp/8.0_wp
        expected_state=0.0_wp
        expected_state(1:3,1)=[11.0_wp/1000.0_wp,-1.0_wp/250.0_wp,-1.0_wp/250.0_wp]
        expected_state(1:3,2)=1.0_wp/160.0_wp
        expected_state(nstate,:)=[1.0_wp/100.0_wp,1.0_wp/16.0_wp]
        expected_derivative=0.0_wp
        expected_derivative(1:3,1,1)=[25.0_wp/34.0_wp,-5.0_wp/34.0_wp,-5.0_wp/34.0_wp]
        expected_derivative(1:3,2,1)=[-269.0_wp/34.0_wp,13.0_wp/34.0_wp,13.0_wp/34.0_wp]
        expected_derivative(1:3,1,2)=5.0_wp/4.0_wp
        expected_derivative(1:3,2,2)=-325.0_wp/16.0_wp
        expected_state_derivative=0.0_wp
        expected_state_derivative(1:3,1,1)=[-11.0_wp/2040.0_wp,1.0_wp/510.0_wp,1.0_wp/510.0_wp]
        expected_state_derivative(1:3,2,1)=[551.0_wp/10200.0_wp,-77.0_wp/5100.0_wp,-77.0_wp/5100.0_wp]
        expected_state_derivative(1:3,1,2)=-1.0_wp/240.0_wp
        expected_state_derivative(1:3,2,2)=13.0_wp/192.0_wp
        expected_state_derivative(nstate,1,:)=[-1.0_wp/204.0_wp,-1.0_wp/24.0_wp]
        expected_state_derivative(nstate,2,:)=[47.0_wp/1020.0_wp,15.0_wp/32.0_wp]
        do load_step=1,2
            if (load_step==2) then
                ! Both points unload elastically with accepted plastic histories unchanged.
                strain=0.5_wp*strain
                expected(1:3,1)=[209.0_wp/50.0_wp,124.0_wp/25.0_wp,124.0_wp/25.0_wp]
                expected(1:3,2)=45.0_wp/8.0_wp
            end if
            call material%eval(strain,response,jacobian,result,dsdp=derivative)
            call expect(result,cf_ok)
            do point=1,2
                do component=1,nmandel
                    call close(response(component,point)/modulus_scale,expected(component,point),reference_tolerance)
                    call close(derivative(component,1,point),expected_derivative(component,1,point),reference_tolerance)
                    call close(derivative(component,2,point)/modulus_scale, &
                        expected_derivative(component,2,point),reference_tolerance)
                end do
                if (load_step==2) then
                    call close(jacobian(1,1,point)/modulus_scale,180.0_wp,reference_tolerance)
                    call close(jacobian(1,2,point)/modulus_scale,60.0_wp,reference_tolerance)
                    call close(jacobian(4,4,point)/modulus_scale,120.0_wp,reference_tolerance)
                end if
            end do
            if (load_step==1) then
                call close(jacobian(1,1,1)/modulus_scale,1185.0_wp/17.0_wp,reference_tolerance)
                call close(jacobian(1,2,1)/modulus_scale,1395.0_wp/17.0_wp,reference_tolerance)
                call close(jacobian(4,4,1)/modulus_scale,240.0_wp/17.0_wp,reference_tolerance)
                call close(jacobian(1,1,2)/modulus_scale,125.0_wp/2.0_wp,reference_tolerance)
                call close(jacobian(1,2,2)/modulus_scale,125.0_wp/2.0_wp,reference_tolerance)
                call close(jacobian(4,4,2)/modulus_scale,0.0_wp,reference_tolerance)
            end if
            call material%commit(result)
            call expect(result,cf_ok)
            call material%export_committed(state,result,sensitivity_state=state_derivative)
            call expect(result,cf_ok)
            do point=1,2
                do component=1,nstate
                    call close(state(component,point),expected_state(component,point),reference_tolerance)
                    call close(state_derivative(component,1,point)*modulus_scale, &
                        expected_state_derivative(component,1,point),reference_tolerance)
                    call close(state_derivative(component,2,point), &
                        expected_state_derivative(component,2,point),reference_tolerance)
                end do
            end do
        end do
    end subroutine test_exact_branches

    subroutine test_pressure_limits()
        type(material_t) :: material
        real(wp) :: strain(nmandel,1),response(nmandel,1),jacobian(nmandel,nmandel,1)
        real(wp) :: derivative(nmandel,3,1)
        real(wp),parameter :: large_bulk=100.0_wp/epsilon(1.0_wp)
        integer :: result,component

        ! With beta=0, hydrostatic pressure cannot change yield or derivative availability.
        call material%init(drucker_prager_t(large_bulk,60.0_wp,3.0_wp,0.0_wp,15.0_wp),1,result, &
            sensitivity_parameters=[character(len=12) :: "yield_stress"])
        call expect(result,cf_ok)
        strain=0.0_wp
        strain(1:3,1)=0.01_wp
        strain(6,1)=sqrt(2.0_wp/3.0_wp)/10.0_wp
        call material%eval(strain,response,jacobian,result,dsdp=derivative(:,1:1,:))
        call expect(result,cf_ok)
        do component=1,3
            call close(response(component,1)/large_bulk,0.03_wp,reference_tolerance)
            call close(derivative(component,1,1),0.0_wp,reference_tolerance)
        end do
        call close(response(6,1),(48.0_wp/13.0_wp)*sqrt(2.0_wp/3.0_wp),reference_tolerance)
        call close(jacobian(6,6,1),120.0_wp/13.0_wp,reference_tolerance)
        call close(derivative(6,1,1),(12.0_wp/13.0_wp)*sqrt(2.0_wp/3.0_wp),reference_tolerance)

        ! At the perfectly plastic apex, p=Y/beta is independent of the large elastic trial pressure.
        call material%init(drucker_prager_t(large_bulk,60.0_wp,0.01_wp,0.5_wp,0.0_wp),1,result, &
            sensitivity_parameters=[character(len=14) :: "yield_stress", "pressure_slope", "bulk_modulus"])
        call expect(result,cf_ok)
        strain=0.0_wp
        strain(1:3,1)=0.001_wp
        call material%eval(strain,response,jacobian,result,dsdp=derivative)
        call expect(result,cf_ok)
        do component=1,3
            call close(response(component,1)/0.02_wp,1.0_wp,reference_tolerance)
            call close(derivative(component,1,1),2.0_wp,reference_tolerance)
            call close(derivative(component,2,1)/0.04_wp,-1.0_wp,reference_tolerance)
            call close(derivative(component,3,1),0.0_wp,reference_tolerance)
        end do
        if (.not. all(ieee_is_finite(jacobian)).or.maxval(abs(jacobian))>reference_tolerance) then
            error stop "Perfectly plastic DP apex tangent"
        end if
    end subroutine test_pressure_limits


    subroutine test_branch_switches()
        type(material_t) :: ordinary,tracked
        real(wp) :: strain(nmandel,1),response(nmandel,1),jacobian(nmandel,nmandel,1),derivative(nmandel,1,1)
        real(wp) :: expected_mean
        integer :: result,branch,component

        call ordinary%init(drucker_prager_t(100.0_wp,60.0_wp,3.0_wp,0.3_wp,15.0_wp),1,result)
        call expect(result,cf_ok)
        call tracked%init(drucker_prager_t(100.0_wp,60.0_wp,3.0_wp,0.3_wp,15.0_wp),1,result, &
            sensitivity_parameters=[character(len=12) :: "yield_stress"])
        call expect(result,cf_ok)
        do branch=1,2
            strain=0.0_wp
            ! Elastic/plastic switch, then the smooth-cone/apex switch.
            strain(1:3,1)=1.0_wp/30.0_wp
            if (branch==2) strain(1:3,1)=[0.046_wp,0.031_wp,0.031_wp]
            expected_mean=10.0_wp+0.5_wp*real(branch-1,wp)
            call ordinary%eval(strain,response,status=result)
            call expect(result,cf_ok)
            do component=1,3
                call close(response(component,1),expected_mean,reference_tolerance)
            end do
            do component=4,nmandel
                call close(response(component,1),0.0_wp,reference_tolerance)
            end do
            ! The stress is unique; classical tangent and parameter derivatives are not.
            call ordinary%eval(strain,response,jacobian,result)
            call expect(result,cf_nondifferentiable)
            call tracked%eval(strain,response,status=result,dsdp=derivative)
            call expect(result,cf_nondifferentiable)
        end do
    end subroutine test_branch_switches

    subroutine test_path_sensitivities()
        integer,parameter :: nparameters=5,npoints=2,nsteps=4
        character(len=18),parameter :: names(nparameters)=[character(len=18) :: "bulk_modulus", &
            "shear_modulus","yield_stress","pressure_slope","hardening_modulus"]
        real(wp),parameter :: parameters(nparameters)=[100.0_wp,60.0_wp,3.0_wp,0.3_wp,15.0_wp]
        real(wp),parameter :: state_modulus=min(2.0_wp*parameters(2),parameters(1)*parameters(4)**2+parameters(5))
        type(material_t) :: tracked
        type(material_t) :: plus(nparameters,2),minus(nparameters,2)
        real(wp) :: path(nmandel,npoints,nsteps),strain(nmandel,npoints),response(nmandel,npoints)
        real(wp) :: derivative(nmandel,nparameters,npoints),state(nstate,npoints)
        real(wp) :: state_derivative(nstate,nparameters,npoints),trial_dev(nmandel,npoints)
        real(wp) :: pplus(nparameters),pminus(nparameters),steps(nparameters,2)
        real(wp) :: yp(nmandel,npoints),ym(nmandel,npoints),hp(nstate,npoints),hm(nstate,npoints)
        real(wp) :: stress_difference(nmandel,npoints,2),state_difference(nstate,npoints,2)
        real(wp) :: stress_tolerance(nmandel,npoints),state_tolerance(nstate,npoints)
        real(wp) :: mean_trial(npoints),solve_scale(npoints),accumulated_scale(npoints),roundoff(npoints)
        integer :: result,load_step,parameter,refinement,point

        ! Accepted branches: [cone,apex], [apex,cone], [elastic,elastic], then [cone,cone].
        ! The last step also changes the deviatoric flow direction.
        path=0.0_wp
        path(1:3,1,1)=[67.0_wp/1500.0_wp,83.0_wp/3000.0_wp,83.0_wp/3000.0_wp]
        path(1:3,2,1)=0.05_wp
        path(1:3,1,2)=0.075_wp
        path(:,2,2)=[0.02_wp,-0.01_wp,-0.01_wp,0.06_wp,0.02_wp,0.0_wp]
        path(:,:,3)=0.95_wp*path(:,:,2)
        path(:,1,4)=[0.02_wp,-0.01_wp,-0.01_wp,0.05_wp,-0.02_wp,0.01_wp]
        path(:,2,4)=[-0.014_wp,0.007_wp,0.007_wp,-0.035_wp,0.025_wp,0.018_wp]
        call tracked%init(drucker_prager_t(parameters(1),parameters(2),parameters(3),parameters(4), &
            parameters(5)),npoints,result,sensitivity_parameters=names)
        call expect(result,cf_ok)
        do parameter=1,nparameters
            ! Balance closed-form arithmetic roundoff against the O(h^2) secant error.
            steps(parameter,1)=(32.0_wp*epsilon(1.0_wp))**(1.0_wp/3.0_wp)*parameters(parameter)
            steps(parameter,2)=0.5_wp*steps(parameter,1)
            do refinement=1,2
                pplus=parameters
                pminus=parameters
                pplus(parameter)=pplus(parameter)+steps(parameter,refinement)
                pminus(parameter)=pminus(parameter)-steps(parameter,refinement)
                call plus(parameter,refinement)%init(drucker_prager_t(pplus(1),pplus(2),pplus(3),pplus(4), &
                    pplus(5)),npoints,result)
                call expect(result,cf_ok)
                call minus(parameter,refinement)%init(drucker_prager_t(pminus(1),pminus(2),pminus(3),pminus(4), &
                    pminus(5)),npoints,result)
                call expect(result,cf_ok)
            end do
        end do
        state=0.0_wp
        accumulated_scale=0.0_wp
        do load_step=1,nsteps
            strain=path(:,:,load_step)
            trial_dev=strain-state(1:nmandel,:)
            mean_trial=parameters(1)*sum(trial_dev(1:3,:),dim=1)
            do point=1,npoints
                trial_dev(1:3,point)=trial_dev(1:3,point)-sum(trial_dev(1:3,point))/3.0_wp
            end do
            trial_dev=2.0_wp*parameters(2)*trial_dev
            solve_scale=max(abs(mean_trial),sqrt(1.5_wp*sum(trial_dev**2,dim=1)), &
                parameters(3)+parameters(5)*state(nstate,:))
            accumulated_scale=accumulated_scale+solve_scale
            call tracked%eval(strain,response,status=result,dsdp=derivative)
            call expect(result,cf_ok)
            call tracked%commit(result)
            call expect(result,cf_ok)
            call tracked%export_committed(state,result,sensitivity_state=state_derivative)
            call expect(result,cf_ok)
            if (.not. all(ieee_is_finite(derivative)).or..not. all(ieee_is_finite(state_derivative))) then
                error stop "Nonfinite Drucker-Prager history derivatives"
            end if
            do parameter=1,nparameters
                do refinement=1,2
                    call plus(parameter,refinement)%eval(strain,yp,status=result)
                    call expect(result,cf_ok)
                    call minus(parameter,refinement)%eval(strain,ym,status=result)
                    call expect(result,cf_ok)
                    call plus(parameter,refinement)%commit(result)
                    call expect(result,cf_ok)
                    call minus(parameter,refinement)%commit(result)
                    call expect(result,cf_ok)
                    call plus(parameter,refinement)%export_committed(hp,result)
                    call expect(result,cf_ok)
                    call minus(parameter,refinement)%export_committed(hm,result)
                    call expect(result,cf_ok)
                    stress_difference(:,:,refinement)=(yp-ym)/(2.0_wp*steps(parameter,refinement))
                    state_difference(:,:,refinement)=(hp-hm)/(2.0_wp*steps(parameter,refinement))
                end do
                if (.not. all(ieee_is_finite(stress_difference)).or..not. all(ieee_is_finite(state_difference))) then
                    error stop "Nonfinite Drucker-Prager history secants"
                end if
                ! Include accepted-path roundoff as well as the measured central-secant truncation.
                roundoff=32.0_wp*epsilon(1.0_wp)*accumulated_scale/steps(parameter,2)
                stress_tolerance=abs(stress_difference(:,:,2)-stress_difference(:,:,1))/3.0_wp+ &
                    spread(roundoff,1,nmandel)
                state_tolerance=abs(state_difference(:,:,2)-state_difference(:,:,1))/3.0_wp+ &
                    spread(roundoff/state_modulus,1,nstate)
                if (any(abs(stress_difference(:,:,2)-derivative(:,parameter,:))>stress_tolerance)) then
                    print *, "Drucker-Prager response history derivative",load_step,parameter
                    error stop "Drucker-Prager parameter derivative mismatch"
                end if
                if (any(abs(state_difference(:,:,2)-state_derivative(:,parameter,:))>state_tolerance)) then
                    print *, "Drucker-Prager internal history derivative",load_step,parameter
                    error stop "Drucker-Prager state derivative mismatch"
                end if
            end do
        end do
    end subroutine test_path_sensitivities

    subroutine expect(actual,expected)
        integer,intent(in) :: actual,expected
        if (actual/=expected) then
            print *, "DP status",actual,expected
            error stop "DP status mismatch"
        end if
    end subroutine expect

    subroutine close(actual,expected,relative)
        real(wp),intent(in) :: actual,expected,relative
        if (.not. ieee_is_finite(actual).or..not. ieee_is_finite(expected)) then
            error stop "Nonfinite Drucker-Prager comparison"
        end if
        if (abs(actual-expected)>relative*max(1.0_wp,abs(expected))) then
            print *, "DP values",actual,expected
            error stop "Drucker-Prager reference mismatch"
        end if
    end subroutine close

end program test_drucker_prager

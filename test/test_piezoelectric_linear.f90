! SPDX-License-Identifier: Apache-2.0

program test_piezoelectric_linear
    !! Test convention: moduli and stress in Pa, electric field in V/m,
    !! electric displacement and piezoelectric coefficients in C/m2,
    !! permittivity in F/m and electric enthalpy per volume in J/m3.
    !!
    !! Coupled stress-charge references, Mandel signs, enthalpy and derivatives.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
    use constituf, only: wp, material_t, cf_ok, cf_invalid_argument, cf_evaluation_failed, cf_no_trial
    use constituf_piezoelectric_linear, only: piezoelectric_linear_t
    implicit none (type, external)

    integer, parameter :: nstrain=6, ninput=9, noutput=10, npoints=2, nparameters=7
    real(wp), parameter :: base(nparameters) = [150000.0_wp,75000.0_wp,-2.3_wp,4.2_wp,3.5_wp,1.2e-8_wp,2.0e-8_wp]
    real(wp), parameter :: reference_values(nparameters) = [100.0_wp,60.0_wp,-2.0_wp,3.0_wp,sqrt(2.0_wp),4.0_wp,6.0_wp]
    real(wp), parameter :: reference_input(ninput) = &
        [0.02_wp,-0.01_wp,0.03_wp,0.04_wp,-0.02_wp,0.01_wp,0.5_wp,-0.25_wp,0.75_wp]
    real(wp), parameter :: reference_response(noutput) = &
        [63.0_wp/10.0_wp,27.0_wp/10.0_wp,15.0_wp/4.0_wp,53.0_wp/10.0_wp,-17.0_wp/5.0_wp,6.0_wp/5.0_wp, &
         49.0_wp/25.0_wp,-23.0_wp/25.0_wp,457.0_wp/100.0_wp,-2067.0_wp/1000.0_wp]
    integer, parameter :: selected_ids(nparameters) = [7,3,2,6,4,5,1]
    character(len=24), parameter :: names(nparameters) = [character(len=24) :: &
        "axial_permittivity", "e31", "shear_modulus", "transverse_permittivity", "e33", "e15", "bulk_modulus"]

    call test_exact_coupling()
    call test_mixed_input_differences()
    call test_uncoupled_and_poling_limits()
    call test_large_shear_modulus()
    call test_deviatoric_accuracy()
    call test_electric_energy_range()
    call test_large_piezoelectric_coupling()
    call test_subnormal_response_products()
    call test_invalid_and_recovery()

contains

    function model_from_values(values) result(model)
        real(wp), intent(in) :: values(nparameters)
        type(piezoelectric_linear_t) :: model

        model=piezoelectric_linear_t(values(1),values(2),values(3),values(4),values(5),values(6),values(7))
    end function model_from_values

    subroutine expect(actual,wanted)
        integer, intent(in) :: actual,wanted

        if (actual/=wanted) then
            print *, "actual status",actual,"wanted",wanted
            error stop "Piezoelectric status mismatch"
        end if
    end subroutine expect

    subroutine close(actual,wanted)
        !! Relative reference comparisons retain small energy derivatives.
        !! All exact zero entries below are structural zeros of the law.
        real(wp), intent(in) :: actual,wanted

        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(wanted)) then
            error stop "Nonfinite piezoelectric comparison"
        end if
        if (abs(actual-wanted)>256.0_wp*epsilon(1.0_wp)*abs(wanted)) then
            print *, "actual",actual,"wanted",wanted
            error stop "Piezoelectric reference mismatch"
        end if
    end subroutine close

    subroutine test_exact_coupling()
        type(material_t) :: material
        real(wp) :: x(ninput,1),y(noutput,1),tangent(noutput,ninput,1)
        real(wp) :: dsdp(noutput,nparameters,1),expected_tangent(noutput,ninput)
        real(wp) :: parameter_reference(noutput,nparameters)
        integer :: status,i,j,k

        call material%init(model_from_values(reference_values),1,status,sensitivity_parameters=names)
        call expect(status,cf_ok)
        if (material%input_size()/=ninput .or. material%output_size()/=noutput) then
            error stop "Piezoelectric packed extents"
        end if
        if (material%output_offset("energy.electric_enthalpy")/=noutput) error stop "Missing electric enthalpy"
        if (material%state_storage_bytes()/=0) error stop "Unexpected piezoelectric history"
        x(:,1)=reference_input
        call material%eval(x,y,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        do i=1,noutput
            call close(y(i,1),reference_response(i))
        end do

        expected_tangent=0.0_wp
        expected_tangent(1:3,1:3)=60.0_wp
        do i=1,3
            expected_tangent(i,i)=180.0_wp
        end do
        do i=4,nstrain
            expected_tangent(i,i)=120.0_wp
        end do
        expected_tangent(1:2,9)=2.0_wp
        expected_tangent(3,9)=-3.0_wp
        expected_tangent(4,8)=-2.0_wp
        expected_tangent(5,7)=-2.0_wp
        expected_tangent(7,5)=2.0_wp
        expected_tangent(8,4)=2.0_wp
        expected_tangent(9,1:2)=-2.0_wp
        expected_tangent(9,3)=3.0_wp
        expected_tangent(7,7)=4.0_wp
        expected_tangent(8,8)=4.0_wp
        expected_tangent(9,9)=6.0_wp
        expected_tangent(noutput,1:nstrain)=reference_response(1:nstrain)
        expected_tangent(noutput,nstrain+1:ninput)=-reference_response(nstrain+1:ninput)
        do j=1,ninput
            do i=1,noutput
                call close(tangent(i,j,1),expected_tangent(i,j))
            end do
        end do
        do j=1,nstrain
            do i=nstrain+1,ninput
                call close(tangent(i,j,1),-tangent(j,i,1))
            end do
        end do
        ! H is an electric enthalpy and may be negative. Its Legendre
        ! transform H+E.D is the positive stored energy for this state.
        call close(y(noutput,1)+dot_product(x(7:9,1),y(7:9,1)),5141.0_wp/2000.0_wp)

        parameter_reference=0.0_wp
        parameter_reference(1:3,1)=1.0_wp/25.0_wp
        parameter_reference(noutput,1)=1.0_wp/1250.0_wp
        parameter_reference(1:nstrain,2)= &
            [1.0_wp/75.0_wp,-7.0_wp/150.0_wp,1.0_wp/30.0_wp,2.0_wp/25.0_wp,-1.0_wp/25.0_wp,1.0_wp/50.0_wp]
        parameter_reference(noutput,2)=89.0_wp/30000.0_wp
        parameter_reference(1:2,3)=-3.0_wp/4.0_wp
        parameter_reference(9,3)=1.0_wp/100.0_wp
        parameter_reference(noutput,3)=-3.0_wp/400.0_wp
        parameter_reference(3,4)=-3.0_wp/4.0_wp
        parameter_reference(9,4)=3.0_wp/100.0_wp
        parameter_reference(noutput,4)=-9.0_wp/400.0_wp
        parameter_reference(4,5)=sqrt(2.0_wp)/4.0_wp
        parameter_reference(5,5)=-sqrt(2.0_wp)/2.0_wp
        parameter_reference(7,5)=-sqrt(2.0_wp)/50.0_wp
        parameter_reference(8,5)=sqrt(2.0_wp)/25.0_wp
        parameter_reference(noutput,5)=sqrt(2.0_wp)/50.0_wp
        parameter_reference(7,6)=0.5_wp
        parameter_reference(8,6)=-0.25_wp
        parameter_reference(noutput,6)=-5.0_wp/32.0_wp
        parameter_reference(9,7)=0.75_wp
        parameter_reference(noutput,7)=-9.0_wp/32.0_wp
        do k=1,nparameters
            do i=1,noutput
                call close(dsdp(i,k,1),parameter_reference(i,selected_ids(k)))
            end do
        end do
    end subroutine test_exact_coupling

    subroutine test_mixed_input_differences()
        type(material_t) :: tracked,plain,perturbed
        real(wp) :: x(ninput,npoints),y(noutput,npoints),tangent(noutput,ninput,npoints)
        real(wp) :: dsdp(noutput,nparameters,npoints),p(nparameters),changed(ninput,npoints)
        real(wp) :: plus(noutput,npoints),minus(noutput,npoints)
        real(wp) :: h,numerical,bound
        integer :: status,i,j,k,ipoint

        call tracked%init(model_from_values(base),npoints,status,sensitivity_parameters=names)
        call expect(status,cf_ok)
        call plain%init(model_from_values(base),npoints,status)
        call expect(status,cf_ok)
        x(:,1)=[0.006_wp,-0.003_wp,0.001_wp,0.0008_wp,-0.001_wp,0.0012_wp,200.0_wp,-150.0_wp,400.0_wp]
        x(:,2)=[-0.004_wp,0.002_wp,0.0005_wp,-0.0006_wp,0.0003_wp,0.001_wp,-70.0_wp,120.0_wp,-90.0_wp]
        call tracked%eval(x,y,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(tangent)) .or. .not. all(ieee_is_finite(dsdp))) then
            error stop "Nonfinite piezoelectric derivatives"
        end if

        ! Stress and displacement are affine; enthalpy is quadratic.
        ! Central differences therefore have no truncation term. Separate
        ! mechanical/electrical input scales avoid mixed-unit step choices.
        do ipoint=1,npoints
            do j=1,ninput
                if (j<=nstrain) then
                    h=0.25_wp*maxval(abs(x(1:nstrain,ipoint)))
                else
                    h=0.25_wp*maxval(abs(x(nstrain+1:ninput,ipoint)))
                end if
                changed=x
                changed(j,ipoint)=x(j,ipoint)+h
                call plain%eval(changed,plus,status=status)
                call expect(status,cf_ok)
                changed(j,ipoint)=x(j,ipoint)-h
                call plain%eval(changed,minus,status=status)
                call expect(status,cf_ok)
                if (.not. all(ieee_is_finite(plus)) .or. .not. all(ieee_is_finite(minus))) then
                    error stop "Nonfinite piezoelectric input secant response"
                end if
                do i=1,noutput
                    numerical=(plus(i,ipoint)-minus(i,ipoint))/(2.0_wp*h)
                    bound=64.0_wp*epsilon(1.0_wp)*(abs(plus(i,ipoint))+abs(minus(i,ipoint)))/(2.0_wp*h)
                    if (.not. all(ieee_is_finite([numerical,bound]))) then
                        error stop "Nonfinite piezoelectric input derivative comparison"
                    end if
                    if (abs(numerical-tangent(i,j,ipoint))>bound) then
                        error stop "Piezoelectric mixed tangent finite difference"
                    end if
                end do
            end do
        end do

        ! Every material coefficient enters all ten outputs affinely.
        do k=1,nparameters
            p=base
            h=0.25_wp*abs(base(selected_ids(k)))
            p(selected_ids(k))=base(selected_ids(k))+h
            call perturbed%init(model_from_values(p),npoints,status)
            call expect(status,cf_ok)
            call perturbed%eval(x,plus,status=status)
            call expect(status,cf_ok)
            p(selected_ids(k))=base(selected_ids(k))-h
            call perturbed%init(model_from_values(p),npoints,status)
            call expect(status,cf_ok)
            call perturbed%eval(x,minus,status=status)
            call expect(status,cf_ok)
            if (.not. all(ieee_is_finite(plus)) .or. .not. all(ieee_is_finite(minus))) then
                error stop "Nonfinite piezoelectric parameter secant response"
            end if
            do ipoint=1,npoints
                do i=1,noutput
                    numerical=(plus(i,ipoint)-minus(i,ipoint))/(2.0_wp*h)
                    bound=64.0_wp*epsilon(1.0_wp)*(abs(plus(i,ipoint))+abs(minus(i,ipoint)))/(2.0_wp*h)
                    if (.not. all(ieee_is_finite([numerical,bound]))) then
                        error stop "Nonfinite piezoelectric parameter derivative comparison"
                    end if
                    if (abs(numerical-dsdp(i,k,ipoint))>bound) then
                        error stop "Piezoelectric parameter finite difference"
                    end if
                end do
            end do
        end do
        call tracked%commit(status)
        call expect(status,cf_ok)
        call tracked%eval(x,y,status=status)
        call expect(status,cf_invalid_argument)
        call tracked%eval(x,y,status=status,dsdp=dsdp)
        call expect(status,cf_ok)
        call tracked%revert()
        call tracked%commit(status)
        call expect(status,cf_no_trial)
    end subroutine test_mixed_input_differences

    subroutine test_uncoupled_and_poling_limits()
        type(material_t) :: uncoupled,reversed
        real(wp), parameter :: uncoupled_response(noutput) = &
            [4.8_wp,1.2_wp,6.0_wp,4.8_wp,-2.4_wp,1.2_wp,2.0_wp,-1.0_wp,4.5_wp,-4109.0_wp/2000.0_wp]
        real(wp) :: p(nparameters),x(ninput,1),y(noutput,1),reverse_y(noutput,1)
        integer :: status,i

        p=reference_values
        p(3:5)=0.0_wp
        call uncoupled%init(model_from_values(p),1,status)
        call expect(status,cf_ok)
        x(:,1)=reference_input
        call uncoupled%eval(x,y,status=status)
        call expect(status,cf_ok)
        p=reference_values
        p(3:5)=-p(3:5)
        call reversed%init(model_from_values(p),1,status)
        call expect(status,cf_ok)
        call reversed%eval(x,reverse_y,status=status)
        call expect(status,cf_ok)
        do i=1,noutput
            call close(y(i,1),uncoupled_response(i))
            call close(0.5_wp*(reference_response(i)+reverse_y(i,1)),uncoupled_response(i))
        end do
    end subroutine test_uncoupled_and_poling_limits

    subroutine test_large_shear_modulus()
        !! Four times G overflows, but the actual normal and shear stiffness
        !! entries and the prescribed shear response are representable.
        type(material_t) :: material
        real(wp) :: x(ninput,1),y(noutput,1),tangent(noutput,ninput,1)
        real(wp) :: shear
        integer :: status

        shear=huge(1.0_wp)/3.5_wp
        call material%init(piezoelectric_linear_t(1.0_wp,shear,0.0_wp,0.0_wp,0.0_wp,1.0_wp,1.0_wp),1,status)
        call expect(status,cf_ok)
        x=0.0_wp
        x(4,1)=0.125_wp
        call material%eval(x,y,tangent,status)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(tangent))) error stop "Nonfinite large-modulus tangent"
        call close(y(4,1)/shear,0.25_wp)
        call close(y(noutput,1)/shear,1.0_wp/64.0_wp)
        call close(tangent(4,4,1)/shear,2.0_wp)
        call close(tangent(1,1,1)/shear,4.0_wp/3.0_wp)
    end subroutine test_large_shear_modulus

    subroutine test_deviatoric_accuracy()
        type(material_t) :: material
        real(wp), parameter :: dilation(2)=[0.125_wp,0.1_wp]
        real(wp) :: x(ninput,1),y(noutput,1),dsdp(noutput,1,1)
        real(wp) :: bulk,h,roundoff
        integer :: status,i,icase

        roundoff=epsilon(1.0_wp)
        bulk=roundoff**2
        call material%init(piezoelectric_linear_t(bulk,1.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp,1.0_wp),1,status, &
            sensitivity_parameters=[character(len=13) :: "shear_modulus"])
        call expect(status,cf_ok)
        do icase=1,size(dilation)
            h=dilation(icase)
            x=0.0_wp
            x(1:3,1)=h
            call material%eval(x,y,status=status,dsdp=dsdp)
            call expect(status,cf_ok)
            do i=1,3
                call close(y(i,1)/(bulk*h),3.0_wp)
            end do
            call close(y(noutput,1)/(bulk*h*h),4.5_wp)
            if (.not. all(ieee_is_finite(dsdp))) error stop "Nonfinite hydrostatic piezoelectric derivative"
            call close(maxval(abs(dsdp)),0.0_wp)
        end do

        ! A tiny Mandel shear on a hydrostatic background has a nonzero
        ! shear-energy derivative that invariant subtraction would erase.
        call material%init(piezoelectric_linear_t(100.0_wp,60.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp,1.0_wp),1,status, &
            sensitivity_parameters=[character(len=13) :: "shear_modulus"])
        call expect(status,cf_ok)
        x(1:3,1)=0.25_wp
        x(4,1)=roundoff
        call material%eval(x,y,status=status,dsdp=dsdp)
        call expect(status,cf_ok)
        call close(dsdp(4,1,1)/roundoff,2.0_wp)
        call close(dsdp(noutput,1,1)/(roundoff**2),1.0_wp)
        if (.not. all(ieee_is_finite(dsdp))) error stop "Nonfinite deviatoric piezoelectric derivative"
        call close(maxval(abs(dsdp(1:3,1,1))),0.0_wp)
    end subroutine test_deviatoric_accuracy

    subroutine test_electric_energy_range()
        type(material_t) :: material
        real(wp) :: x(ninput,1),y(noutput,1),tangent(noutput,ninput,1)
        real(wp) :: dsdp(noutput,nparameters,1)
        real(wp) :: large,field_scale
        integer :: status

        large=huge(1.0_wp)
        call material%init(piezoelectric_linear_t(1.0_wp,1.0_wp,0.0_wp,0.0_wp,0.0_wp,0.75_wp*large,1.0_wp), &
            1,status,sensitivity_parameters=names)
        call expect(status,cf_ok)
        x=0.0_wp
        x(7:8,1)=1.0_wp
        call material%eval(x,y,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(tangent)) .or. .not. all(ieee_is_finite(dsdp))) then
            error stop "Nonfinite electric-energy range derivatives"
        end if
        call close(y(7,1)/large,0.75_wp)
        call close(y(8,1)/large,0.75_wp)
        call close(y(noutput,1)/large,-0.75_wp)
        call close(tangent(7,7,1)/large,0.75_wp)
        call close(dsdp(noutput,4,1),-1.0_wp)

        ! The raw field square overflows; one half of that square, needed
        ! by the selected axial-permittivity energy derivative, is finite.
        field_scale=sqrt(large)
        call material%init(piezoelectric_linear_t(1.0_wp,1.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp,1.0_wp/large), &
            1,status,sensitivity_parameters=names)
        call expect(status,cf_ok)
        x=0.0_wp
        x(9,1)=field_scale*sqrt(1.5_wp)
        call material%eval(x,y,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(tangent)) .or. .not. all(ieee_is_finite(dsdp))) then
            error stop "Nonfinite squared-field derivative"
        end if
        call close(y(9,1)*field_scale,sqrt(1.5_wp))
        call close(y(noutput,1),-0.75_wp)
        call close(dsdp(noutput,1,1)/large,-0.75_wp)
        call close(dsdp(9,1,1)/field_scale,sqrt(1.5_wp))
    end subroutine test_electric_energy_range

    subroutine test_large_piezoelectric_coupling()
        !! The Mandel cross-block coefficient is unrepresentable, but its
        !! products with these fields and the selected derivative are finite.
        type(material_t) :: material
        real(wp) :: x(ninput,1),y(noutput,1),tangent(noutput,ninput,1)
        real(wp) :: dsdp(noutput,1,1)
        real(wp) :: large
        integer :: status

        large=huge(1.0_wp)
        call material%init(piezoelectric_linear_t(100.0_wp,60.0_wp,0.0_wp,0.0_wp,0.75_wp*large,1.0_wp,1.0_wp), &
            1,status,sensitivity_parameters=[character(len=3) :: "e15"])
        call expect(status,cf_ok)
        x=0.0_wp
        call material%eval(x,y,status=status,dsdp=dsdp)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(dsdp))) then
            error stop "Nonfinite zero-field piezoelectric response"
        end if
        call close(maxval(abs(y)),0.0_wp)
        call close(maxval(abs(dsdp)),0.0_wp)
        x(4,1)=0.25_wp
        x(8,1)=0.25_wp
        call material%eval(x,y,status=status,dsdp=dsdp)
        call expect(status,cf_ok)
        call close(y(4,1)/large,-3.0_wp*sqrt(2.0_wp)/16.0_wp)
        call close(y(8,1)/large,3.0_wp*sqrt(2.0_wp)/16.0_wp)
        call close(y(noutput,1)/large,-3.0_wp*sqrt(2.0_wp)/64.0_wp)
        call close(dsdp(4,1,1),-sqrt(2.0_wp)/4.0_wp)
        call close(dsdp(8,1,1),sqrt(2.0_wp)/4.0_wp)
        call close(dsdp(noutput,1,1),-sqrt(2.0_wp)/16.0_wp)
        call material%eval(x,y,tangent,status,dsdp=dsdp)
        call expect(status,cf_evaluation_failed)
        if (material%last_failure_point()/=1) error stop "Piezoelectric coupling tangent overflow point"
        call material%commit(status)
        call expect(status,cf_no_trial)
        call material%eval(x,y,status=status,dsdp=dsdp)
        call expect(status,cf_ok)
    end subroutine test_large_piezoelectric_coupling

    subroutine test_subnormal_response_products()
        !! Keep the Mandel/shear factor in the complete product before
        !! rounding. The inner coefficient/input product alone rounds to zero.
        integer, parameter :: nelectric=ninput-nstrain, ndielectric_parameters=2, ndielectric_cases=3
        real(wp), parameter :: electric_cases(nelectric,ndielectric_cases)=reshape( &
            [1.5_wp,0.0_wp,0.0_wp,1.0_wp,1.0_wp,0.0_wp,0.75_wp,0.75_wp,0.75_wp], &
            [nelectric,ndielectric_cases])
        real(wp), parameter :: rounded_displacement(nelectric,ndielectric_cases)=reshape( &
            [2.0_wp,0.0_wp,0.0_wp,1.0_wp,1.0_wp,0.0_wp,1.0_wp,1.0_wp,1.0_wp], &
            [nelectric,ndielectric_cases])
        real(wp), parameter :: transverse_energy_derivative(ndielectric_cases)= &
            [-9.0_wp/8.0_wp,-1.0_wp,-9.0_wp/16.0_wp]
        real(wp), parameter :: axial_energy_derivative(ndielectric_cases)= &
            [0.0_wp,0.0_wp,-9.0_wp/32.0_wp]
        type(material_t) :: material
        real(wp) :: x(ninput,1),y(noutput,1),tangent(noutput,ninput,1)
        real(wp) :: dsdp(noutput,1,1)
        real(wp) :: dielectric_dsdp(noutput,ndielectric_parameters,1)
        real(wp) :: expected_derivatives(noutput,ndielectric_parameters)
        real(wp) :: smallest
        integer :: status,i,j,icase

        smallest=nearest(0.0_wp,1.0_wp)
        call material%init(piezoelectric_linear_t(1.0_wp,1.0_wp,0.0_wp,0.0_wp,smallest,1.0_wp,1.0_wp), &
            1,status,sensitivity_parameters=[character(len=3) :: "e15"])
        call expect(status,cf_ok)
        x=0.0_wp
        x(5,1)=0.5_wp
        call material%eval(x,y,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        call close(y(7,1)/smallest,1.0_wp)
        call close(dsdp(7,1,1),sqrt(2.0_wp)/2.0_wp)

        x=0.0_wp
        x(7,1)=0.5_wp
        call material%eval(x,y,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        call close(y(5,1)/smallest,-1.0_wp)
        call close(dsdp(5,1,1),-sqrt(2.0_wp)/2.0_wp)

        call material%init(piezoelectric_linear_t(1.0_wp,7.0_wp*smallest,0.0_wp,0.0_wp,0.0_wp,1.0_wp,1.0_wp), &
            1,status,sensitivity_parameters=[character(len=13) :: "shear_modulus"])
        call expect(status,cf_ok)
        x=0.0_wp
        x(6,1)=1.0_wp/16.0_wp
        call material%eval(x,y,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        call close(y(6,1)/smallest,1.0_wp)
        call close(tangent(6,6,1)/smallest,14.0_wp)
        call close(dsdp(6,1,1),1.0_wp/8.0_wp)
        call close(dsdp(noutput,1,1),1.0_wp/256.0_wp)

        ! The shear energy is G*gamma_M**2, without first rounding G*gamma_M.
        call material%init(piezoelectric_linear_t(1.0_wp,smallest,0.0_wp,0.0_wp,0.0_wp,1.0_wp,1.0_wp), &
            1,status,sensitivity_parameters=[character(len=13) :: "shear_modulus"])
        call expect(status,cf_ok)
        x=0.0_wp
        x(6,1)=1.5_wp
        call material%eval(x,y,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        call close(y(6,1)/smallest,3.0_wp)
        call close(y(noutput,1)/smallest,2.0_wp)
        call close(dsdp(noutput,1,1),2.25_wp)

        ! Three weighted shear squares of 0.25d must sum before rounding.
        call material%init(piezoelectric_linear_t(1.0_wp,64.0_wp*smallest,0.0_wp,0.0_wp,0.0_wp,1.0_wp,1.0_wp), &
            1,status,sensitivity_parameters=[character(len=13) :: "shear_modulus"])
        call expect(status,cf_ok)
        x=0.0_wp
        x(4:nstrain,1)=1.0_wp/16.0_wp
        call material%eval(x,y,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        do i=4,nstrain
            call close(y(i,1)/smallest,8.0_wp)
            call close(dsdp(i,1,1),1.0_wp/8.0_wp)
        end do
        call close(y(noutput,1)/smallest,1.0_wp)
        call close(dsdp(noutput,1,1),3.0_wp/256.0_wp)

        ! Hydrostatic energy is K*trace(strain)**2/2 with one final rounding.
        call material%init(piezoelectric_linear_t(smallest,1.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp,1.0_wp), &
            1,status,sensitivity_parameters=[character(len=12) :: "bulk_modulus"])
        call expect(status,cf_ok)
        x=0.0_wp
        x(1:3,1)=0.5_wp
        call material%eval(x,y,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        do i=1,3
            call close(y(i,1)/smallest,2.0_wp)
            call close(dsdp(i,1,1),1.5_wp)
        end do
        call close(y(noutput,1)/smallest,1.0_wp)
        call close(dsdp(noutput,1,1),1.125_wp)

        ! Dielectric enthalpy uses the complete weighted field norm before
        ! rounding. Reusing rounded displacement, or rounding each component
        ! energy separately, loses these representable subnormal energies.
        call material%init(piezoelectric_linear_t(1.0_wp,1.0_wp,0.0_wp,0.0_wp,0.0_wp,smallest,smallest), &
            1,status,sensitivity_parameters=[character(len=24) :: "transverse_permittivity","axial_permittivity"])
        call expect(status,cf_ok)
        do icase=1,ndielectric_cases
            x=0.0_wp
            x(nstrain+1:ninput,1)=electric_cases(:,icase)
            call material%eval(x,y,tangent,status,dsdp=dielectric_dsdp)
            call expect(status,cf_ok)
            if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(tangent)) .or. &
                .not. all(ieee_is_finite(dielectric_dsdp))) then
                error stop "Nonfinite subnormal dielectric enthalpy response"
            end if
            call close(y(noutput,1)/smallest,-1.0_wp)
            do i=1,nelectric
                call close(y(nstrain+i,1)/smallest,rounded_displacement(i,icase))
                do j=1,nelectric
                    call close(tangent(nstrain+i,nstrain+j,1)/smallest,merge(1.0_wp,0.0_wp,i==j))
                end do
            end do
            expected_derivatives=0.0_wp
            expected_derivatives(7:8,1)=electric_cases(1:2,icase)
            expected_derivatives(9,2)=electric_cases(3,icase)
            expected_derivatives(noutput,1)=transverse_energy_derivative(icase)
            expected_derivatives(noutput,2)=axial_energy_derivative(icase)
            do j=1,ndielectric_parameters
                do i=1,noutput
                    call close(dielectric_dsdp(i,j,1),expected_derivatives(i,j))
                end do
            end do
        end do
    end subroutine test_subnormal_response_products

    subroutine test_invalid_and_recovery()
        type(material_t) :: material
        real(wp) :: p(nparameters),x(ninput,npoints),y(noutput,npoints)
        integer :: status

        p=base
        p(6)=0.0_wp
        call material%init(model_from_values(p),npoints,status)
        call expect(status,cf_invalid_argument)
        p=base
        p(3)=ieee_value(0.0_wp,ieee_quiet_nan)
        call material%init(model_from_values(p),npoints,status)
        call expect(status,cf_invalid_argument)
        call material%init(model_from_values(base),npoints,status)
        call expect(status,cf_ok)
        x=0.0_wp
        x(7,2)=ieee_value(0.0_wp,ieee_quiet_nan)
        call material%eval(x,y,status=status)
        call expect(status,cf_evaluation_failed)
        if (material%last_failure_point()/=2) error stop "Piezoelectric failed point"
        call material%commit(status)
        call expect(status,cf_no_trial)
        x(7,2)=huge(1.0_wp)/4.0_wp
        call material%eval(x,y,status=status)
        call expect(status,cf_evaluation_failed)
        if (material%last_failure_point()/=2) error stop "Piezoelectric energy overflow point"
        call material%commit(status)
        call expect(status,cf_no_trial)
        x=0.0_wp
        call material%eval(x,y,status=status)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(y))) error stop "Nonfinite piezoelectric recovery response"
        call close(maxval(abs(y)),0.0_wp)
    end subroutine test_invalid_and_recovery

end program test_piezoelectric_linear

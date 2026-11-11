! SPDX-License-Identifier: Apache-2.0

program test_biot_poroelastic
    !! Test convention: moduli, pore pressure and total stress in Pa.
    !! Strain and fluid content are dimensionless.
    !!
    !! Exact coupled response, mixed potential, limits and arithmetic contracts.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
    use constituf, only: wp,material_t,cf_ok,cf_invalid_argument,cf_evaluation_failed,cf_no_trial
    use constituf_biot_poroelastic, only: biot_poroelastic_t
    implicit none (type, external)

    integer, parameter :: nmandel=6,nfields=7,nparameters=4,npoints=2
    real(wp), parameter :: reference_values(nparameters)=[100.0_wp,60.0_wp,0.75_wp,200.0_wp]
    real(wp), parameter :: reference_input(nfields)=[0.02_wp,-0.01_wp,0.03_wp,0.04_wp,-0.02_wp,0.01_wp,4.0_wp]
    real(wp), parameter :: reference_response(nfields)=[9.0_wp/5.0_wp,-9.0_wp/5.0_wp,3.0_wp, &
        24.0_wp/5.0_wp,-12.0_wp/5.0_wp,6.0_wp/5.0_wp,1.0_wp/20.0_wp]
    integer, parameter :: selected_ids(nparameters)=[4,3,2,1]
    character(len=16), parameter :: parameter_names(nparameters)=[character(len=16) :: &
        "biot_modulus","biot_coefficient","shear_modulus","bulk_modulus"]

    call test_exact_reference_and_potential()
    call test_mixed_input_differences()
    call test_drained_undrained_and_uncoupled()
    call test_large_shear_modulus()
    call test_subnormal_response_products()
    call test_hydrostatic_accuracy()
    call test_storage_scale(tiny(1.0_wp))
    call test_storage_scale(huge(1.0_wp)/4.0_wp)
    call test_storage_derivative_precision()
    call test_requested_output_range()
    call test_invalid_and_recovery()

contains

    function model_from_values(values) result(model)
        real(wp), intent(in) :: values(nparameters)
        type(biot_poroelastic_t) :: model

        model=biot_poroelastic_t(values(1),values(2),values(3),values(4))
    end function model_from_values

    subroutine expect(actual,expected)
        integer, intent(in) :: actual,expected

        if (actual/=expected) then
            print *, "Biot status",actual,expected
            error stop "Biot status mismatch"
        end if
    end subroutine expect

    subroutine close(actual,expected,comparison_scale)
        real(wp), intent(in) :: actual,expected
        real(wp), intent(in), optional :: comparison_scale
        real(wp) :: reference_scale

        reference_scale=abs(expected)
        if (present(comparison_scale)) reference_scale=comparison_scale
        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected) .or. &
            .not. ieee_is_finite(reference_scale) .or. &
            abs(actual-expected)>256.0_wp*epsilon(1.0_wp)*reference_scale) then
            print *, "Biot values",actual,expected
            error stop "Biot reference mismatch"
        end if
    end subroutine close

    pure function mixed_potential(values,x) result(potential)
        !! Independent quadratic mixed potential; not an output of this model.
        real(wp), intent(in) :: values(nparameters),x(nfields)
        real(wp) :: potential,trace_strain,deviatoric_norm_squared

        trace_strain=sum(x(1:3))
        deviatoric_norm_squared=dot_product(x(1:nmandel),x(1:nmandel))-trace_strain**2/3.0_wp
        potential=0.5_wp*values(1)*trace_strain**2+values(2)*deviatoric_norm_squared- &
            values(3)*x(nfields)*trace_strain-0.5_wp*x(nfields)**2/values(4)
    end function mixed_potential

    subroutine test_exact_reference_and_potential()
        type(material_t) :: material
        real(wp) :: x(nfields,1),y(nfields,1),tangent(nfields,nfields,1),dsdp(nfields,nparameters,1)
        real(wp) :: expected_tangent(nfields,nfields),parameter_reference(nfields,nparameters)
        real(wp) :: changed(nfields)
        real(wp) :: potential,positive,negative,h,numerical,wanted,bound
        integer :: status,i,j,k

        call material%init(model_from_values(reference_values),1,status,sensitivity_parameters=parameter_names)
        call expect(status,cf_ok)
        if (material%input_size()/=nfields.or.material%output_size()/=nfields) error stop "Biot shape"
        if (material%state_size()/=0.or.material%state_storage_bytes()/=0) error stop "Biot stateless"
        if (.not. material%field_matches("input","pressure.pore","scalar", &
            "pore_pressure","material",[integer ::])) error stop "Biot pressure schema"
        if (.not. material%field_matches("output","fluid.content","scalar", &
            "fluid_content","material",[integer ::])) error stop "Biot content schema"
        if (.not. material%field_matches("input","strain.small","mandel", &
            "small_strain","spatial",[3,3])) error stop "Biot Mandel strain schema"
        if (.not. material%field_matches("output","stress.total","mandel", &
            "cauchy_stress","spatial",[3,3])) error stop "Biot total-stress schema"

        x(:,1)=reference_input
        call material%eval(x,y,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        do i=1,nfields
            call close(y(i,1),reference_response(i))
        end do
        expected_tangent=0.0_wp
        expected_tangent(1:3,1:3)=60.0_wp
        do i=1,3
            expected_tangent(i,i)=180.0_wp
        end do
        do i=4,nmandel
            expected_tangent(i,i)=120.0_wp
        end do
        expected_tangent(1:3,nfields)=-0.75_wp
        expected_tangent(nfields,1:3)=0.75_wp
        expected_tangent(nfields,nfields)=1.0_wp/200.0_wp
        do j=1,nfields
            do i=1,nfields
                call close(tangent(i,j,1),expected_tangent(i,j))
            end do
        end do

        parameter_reference=0.0_wp
        parameter_reference(1:3,1)=1.0_wp/25.0_wp
        parameter_reference(1:nmandel,2)= &
            [1.0_wp/75.0_wp,-7.0_wp/150.0_wp,1.0_wp/30.0_wp,2.0_wp/25.0_wp,-1.0_wp/25.0_wp,1.0_wp/50.0_wp]
        parameter_reference(1:3,3)=-4.0_wp
        parameter_reference(nfields,3)=1.0_wp/25.0_wp
        parameter_reference(nfields,4)=-1.0_wp/10000.0_wp
        do k=1,nparameters
            do i=1,nfields
                call close(dsdp(i,k,1),parameter_reference(i,selected_ids(k)))
            end do
        end do

        potential=0.5_wp*(dot_product(x(1:nmandel,1),y(1:nmandel,1))-x(nfields,1)*y(nfields,1))
        call close(potential,49.0_wp/500.0_wp)
        call close(potential+x(nfields,1)*y(nfields,1),149.0_wp/500.0_wp)
        call close(mixed_potential(reference_values,x(:,1)),49.0_wp/500.0_wp)
        ! The potential gradient is [sigma,-zeta]. Its quadratic form makes
        ! these finite central perturbations exact apart from roundoff.
        do j=1,nfields
            if (j<=nmandel) then
                h=0.25_wp*maxval(abs(x(1:nmandel,1)))
                wanted=y(j,1)
            else
                h=0.25_wp*abs(x(nfields,1))
                wanted=-y(nfields,1)
            end if
            changed=x(:,1)
            changed(j)=changed(j)+h
            positive=mixed_potential(reference_values,changed)
            changed(j)=x(j,1)-h
            negative=mixed_potential(reference_values,changed)
            if (.not. ieee_is_finite(positive) .or. .not. ieee_is_finite(negative)) then
                error stop "Nonfinite perturbed Biot potential"
            end if
            numerical=(positive-negative)/(2.0_wp*h)
            bound=64.0_wp*epsilon(1.0_wp)*(abs(positive)+abs(negative))/(2.0_wp*h)
            if (.not. ieee_is_finite(numerical) .or. .not. ieee_is_finite(bound) .or. &
                abs(numerical-wanted)>bound) then
                error stop "Biot mixed potential gradient"
            end if
        end do
    end subroutine test_exact_reference_and_potential

    subroutine test_mixed_input_differences()
        type(material_t) :: tracked,plain,upper,lower
        real(wp), parameter :: values(nparameters)=[100.0_wp,40.0_wp,0.7_wp,300.0_wp]
        real(wp) :: x(nfields,npoints),y(nfields,npoints),tangent(nfields,nfields,npoints)
        real(wp) :: dsdp(nfields,nparameters,npoints),positive(nfields,npoints),negative(nfields,npoints)
        real(wp) :: changed(nfields,npoints),parameters(nparameters)
        real(wp) :: h,numerical,bound,relative_step
        integer :: status,i,j,k,id,ipoint

        call tracked%init(model_from_values(values),npoints,status,sensitivity_parameters=parameter_names)
        call expect(status,cf_ok)
        call plain%init(model_from_values(values),npoints,status)
        call expect(status,cf_ok)
        x(:,1)=[0.03_wp,-0.01_wp,0.02_wp,0.04_wp,-0.02_wp,0.01_wp,7.0_wp]
        x(:,2)=[-0.02_wp,0.02_wp,0.005_wp,0.0_wp,0.0_wp,0.0_wp,-5.0_wp]
        call tracked%eval(x,y,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(tangent)) .or. .not. all(ieee_is_finite(dsdp))) then
            error stop "Nonfinite Biot derivatives"
        end if

        ! Both input fields enter affinely; use their separate physical scales.
        do ipoint=1,npoints
            do j=1,nfields
                if (j<=nmandel) then
                    h=0.25_wp*maxval(abs(x(1:nmandel,ipoint)))
                else
                    h=0.25_wp*abs(x(nfields,ipoint))
                end if
                changed=x
                changed(j,ipoint)=x(j,ipoint)+h
                call plain%eval(changed,positive,status=status)
                call expect(status,cf_ok)
                if (.not. all(ieee_is_finite(positive))) error stop "Nonfinite positive Biot input perturbation"
                changed(j,ipoint)=x(j,ipoint)-h
                call plain%eval(changed,negative,status=status)
                call expect(status,cf_ok)
                if (.not. all(ieee_is_finite(negative))) error stop "Nonfinite negative Biot input perturbation"
                do i=1,nfields
                    numerical=(positive(i,ipoint)-negative(i,ipoint))/(2.0_wp*h)
                    bound=64.0_wp*epsilon(1.0_wp)*(abs(positive(i,ipoint))+abs(negative(i,ipoint)))/(2.0_wp*h)
                    if (.not. ieee_is_finite(numerical) .or. .not. ieee_is_finite(bound) .or. &
                        abs(numerical-tangent(i,j,ipoint))>bound) then
                        error stop "Biot input tangent finite difference"
                    end if
                end do
            end do
        end do

        do k=1,nparameters
            id=selected_ids(k)
            if (id==4) then
                relative_step=epsilon(1.0_wp)**(1.0_wp/3.0_wp)
                h=values(id)*relative_step
            else
                relative_step=0.0_wp
                h=0.25_wp*values(id)
            end if
            parameters=values
            parameters(id)=values(id)+h
            call upper%init(model_from_values(parameters),npoints,status)
            call expect(status,cf_ok)
            parameters(id)=values(id)-h
            call lower%init(model_from_values(parameters),npoints,status)
            call expect(status,cf_ok)
            call upper%eval(x,positive,status=status)
            call expect(status,cf_ok)
            if (.not. all(ieee_is_finite(positive))) error stop "Nonfinite positive Biot parameter perturbation"
            call lower%eval(x,negative,status=status)
            call expect(status,cf_ok)
            if (.not. all(ieee_is_finite(negative))) error stop "Nonfinite negative Biot parameter perturbation"
            do ipoint=1,npoints
                do i=1,nfields
                    numerical=(positive(i,ipoint)-negative(i,ipoint))/(2.0_wp*h)
                    ! Only M is nonlinear; use its exact reciprocal central bias.
                    bound=64.0_wp*epsilon(1.0_wp)*(abs(positive(i,ipoint))+abs(negative(i,ipoint)))/(2.0_wp*h)+ &
                        abs(dsdp(i,k,ipoint))*relative_step**2/(1.0_wp-relative_step**2)
                    if (.not. ieee_is_finite(numerical) .or. .not. ieee_is_finite(bound) .or. &
                        abs(numerical-dsdp(i,k,ipoint))>bound) then
                        error stop "Biot parameter finite difference"
                    end if
                end do
            end do
        end do
        call tracked%commit(status)
        call expect(status,cf_ok)
        call tracked%revert()
        call tracked%commit(status)
        call expect(status,cf_no_trial)
    end subroutine test_mixed_input_differences

    subroutine test_drained_undrained_and_uncoupled()
        type(material_t) :: material
        real(wp), parameter :: uncoupled_stress(nmandel)=[4.8_wp,1.2_wp,6.0_wp,4.8_wp,-2.4_wp,1.2_wp]
        real(wp) :: x(nfields,1),y(nfields,1),tangent(nfields,nfields,1),dsdp(nfields,1,1)
        real(wp) :: eliminated,wanted
        integer :: status,i,j

        call material%init(model_from_values(reference_values),1,status)
        call expect(status,cf_ok)
        x=0.0_wp
        x(1:3,1)=-0.01_wp
        call material%eval(x,y,tangent,status)
        call expect(status,cf_ok)
        do i=1,3
            call close(y(i,1),-3.0_wp)
        end do
        call close(y(nfields,1),-9.0_wp/400.0_wp)

        ! Enforcing zeta=0 under compression gives positive pore pressure
        ! and the undrained bulk modulus K+b**2*M=425/2.
        x(nfields,1)=4.5_wp
        call material%eval(x,y,tangent,status)
        call expect(status,cf_ok)
        do i=1,3
            call close(y(i,1),-51.0_wp/8.0_wp)
        end do
        call close(y(nfields,1),0.0_wp,0.03_wp)
        do j=1,nmandel
            do i=1,nmandel
                eliminated=tangent(i,j,1)-tangent(i,nfields,1)*tangent(nfields,j,1)/tangent(nfields,nfields,1)
                wanted=0.0_wp
                if (i<=3.and.j<=3) wanted=172.5_wp
                if (i==j) then
                    if (i<=3) then
                        wanted=292.5_wp
                    else
                        wanted=120.0_wp
                    end if
                end if
                call close(eliminated,wanted)
            end do
        end do

        x=0.0_wp
        x(nfields,1)=4.0_wp
        call material%eval(x,y,status=status)
        call expect(status,cf_ok)
        do i=1,3
            call close(y(i,1),-3.0_wp)
        end do
        call close(y(nfields,1),1.0_wp/50.0_wp)

        call material%init(biot_poroelastic_t(100.0_wp,60.0_wp,0.0_wp,200.0_wp),1,status, &
            sensitivity_parameters=[character(len=16) :: "biot_coefficient"])
        call expect(status,cf_ok)
        x(:,1)=reference_input
        call material%eval(x,y,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        do i=1,nmandel
            call close(y(i,1),uncoupled_stress(i))
            call close(tangent(i,nfields,1),0.0_wp)
            call close(tangent(nfields,i,1),0.0_wp)
        end do
        call close(y(nfields,1),1.0_wp/50.0_wp)
        do i=1,3
            call close(dsdp(i,1,1),-4.0_wp)
        end do
        call close(dsdp(nfields,1,1),1.0_wp/25.0_wp)
    end subroutine test_drained_undrained_and_uncoupled

    subroutine test_large_shear_modulus()
        type(material_t) :: material
        real(wp) :: x(nfields,1),y(nfields,1),tangent(nfields,nfields,1)
        real(wp) :: shear
        integer :: status

        shear=huge(1.0_wp)/3.5_wp
        call material%init(biot_poroelastic_t(1.0_wp,shear,0.0_wp,1.0_wp),1,status)
        call expect(status,cf_ok)
        x=0.0_wp
        x(4,1)=0.125_wp
        call material%eval(x,y,tangent,status)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(tangent))) error stop "Nonfinite large-modulus Biot tangent"
        call close(y(4,1)/shear,0.25_wp)
        call close(tangent(4,4,1)/shear,2.0_wp)
        call close(tangent(1,1,1)/shear,4.0_wp/3.0_wp)
        call close(tangent(1,2,1)/shear,-2.0_wp/3.0_wp)
    end subroutine test_large_shear_modulus

    subroutine test_subnormal_response_products()
        !! Retain original factors until subnormal stress and fluid content round.
        !! The shear products [7/8,-7/8,21/8]*d round to [1,-1,3]*d.
        type(material_t) :: material
        integer, parameter :: nselected=2
        real(wp), parameter :: shear_strain(3)=[1.0_wp,-1.0_wp,3.0_wp]/16.0_wp
        real(wp), parameter :: scaled_stress(3)=[1.0_wp,-1.0_wp,3.0_wp]
        real(wp) :: x(nfields,1),y(nfields,1),tangent(nfields,nfields,1),dsdp(nfields,1,1)
        real(wp) :: fluid_gradient(nfields,nselected,1)
        real(wp) :: smallest,sgn
        integer :: status,i,icase

        smallest=nearest(0.0_wp,1.0_wp)
        call material%init(biot_poroelastic_t(1.0_wp,7.0_wp*smallest,0.0_wp,1.0_wp),1,status, &
            sensitivity_parameters=[character(len=13) :: "shear_modulus"])
        call expect(status,cf_ok)
        x=0.0_wp
        x(4:nmandel,1)=shear_strain
        call material%eval(x,y,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(tangent)) .or. &
            .not. all(ieee_is_finite(dsdp))) error stop "Nonfinite subnormal Biot response"
        do i=1,3
            call close(y(i+3,1)/smallest,scaled_stress(i))
            call close(tangent(i+3,i+3,1)/smallest,14.0_wp)
            call close(dsdp(i+3,1,1),2.0_wp*shear_strain(i))
        end do
        call close(maxval(abs(y(1:3,1))),0.0_wp)
        call close(y(nfields,1),0.0_wp)
        call close(maxval(abs(dsdp(1:3,1,1))),0.0_wp)
        call close(dsdp(nfields,1,1),0.0_wp)

        ! Each fluid-content term is signed 3*d/8, so rounding them separately
        ! loses their signed 3*d/4 sum, which rounds to one subnormal quantum.
        call material%init(biot_poroelastic_t(1.0_wp,1.0_wp,32.0_wp*smallest,8.0_wp),1,status, &
            sensitivity_parameters=parameter_names(1:nselected))
        call expect(status,cf_ok)
        do icase=1,2
            sgn=real(3-2*icase,wp)
            x=0.0_wp
            x(1:3,1)=sgn/256.0_wp
            x(nfields,1)=3.0_wp*sgn*smallest
            call material%eval(x,y,tangent,status,dsdp=fluid_gradient)
            call expect(status,cf_ok)
            if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(tangent)) .or. &
                .not. all(ieee_is_finite(fluid_gradient))) error stop "Nonfinite subnormal Biot fluid content"
            call close(y(nfields,1)/smallest,sgn)
            do i=1,3
                call close(y(i,1),3.0_wp*sgn/256.0_wp)
                call close(fluid_gradient(i,2,1)/smallest,-3.0_wp*sgn)
                call close(tangent(i,nfields,1)/smallest,-32.0_wp)
                call close(tangent(nfields,i,1)/smallest,32.0_wp)
            end do
            call close(maxval(abs(y(4:nmandel,1))),0.0_wp)
            call close(maxval(abs(fluid_gradient(:,1,1))),0.0_wp)
            call close(maxval(abs(fluid_gradient(4:nmandel,2,1))),0.0_wp)
            call close(fluid_gradient(nfields,2,1),3.0_wp*sgn/256.0_wp)
            call close(tangent(nfields,nfields,1),1.0_wp/8.0_wp)
        end do
    end subroutine test_subnormal_response_products

    subroutine test_hydrostatic_accuracy()
        type(material_t) :: material
        real(wp), parameter :: dilations(2)=[0.1_wp,0.125_wp]
        real(wp) :: x(nfields,1),y(nfields,1),dsdp(nfields,1,1)
        real(wp) :: bulk,h
        integer :: status,icase,i

        bulk=epsilon(1.0_wp)**2
        call material%init(biot_poroelastic_t(bulk,1.0_wp,0.0_wp,1.0_wp),1,status, &
            sensitivity_parameters=[character(len=13) :: "shear_modulus"])
        call expect(status,cf_ok)
        do icase=1,size(dilations)
            h=dilations(icase)
            x=0.0_wp
            x(1:3,1)=h
            call material%eval(x,y,status=status,dsdp=dsdp)
            call expect(status,cf_ok)
            if (.not. all(ieee_is_finite(dsdp))) error stop "Nonfinite hydrostatic Biot sensitivity"
            do i=1,3
                call close(y(i,1)/(bulk*h),3.0_wp)
            end do
            call close(maxval(abs(dsdp)),0.0_wp)
        end do
    end subroutine test_hydrostatic_accuracy

    subroutine test_storage_scale(storage_modulus)
        !! Both extremes make M**2 invalid although p/M and its derivative
        !! are finite. Normalization makes the small derivative meaningful.
        real(wp), intent(in) :: storage_modulus
        type(material_t) :: material
        real(wp) :: x(nfields,1),y(nfields,1),tangent(nfields,nfields,1),dsdp(nfields,1,1)
        integer :: status

        call material%init(biot_poroelastic_t(1.0_wp,1.0_wp,0.0_wp,storage_modulus),1,status, &
            sensitivity_parameters=[character(len=12) :: "biot_modulus"])
        call expect(status,cf_ok)
        x=0.0_wp
        x(nfields,1)=storage_modulus
        call material%eval(x,y,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(tangent)) .or. &
            .not. all(ieee_is_finite(dsdp))) error stop "Nonfinite Biot storage response"
        call close(maxval(abs(y(1:nmandel,1))),0.0_wp)
        call close(y(nfields,1),1.0_wp)
        call close(tangent(nfields,nfields,1)*storage_modulus,1.0_wp)
        call close(dsdp(nfields,1,1)*storage_modulus,-1.0_wp)
        call close(maxval(abs(dsdp(1:nmandel,1,1))),0.0_wp)
    end subroutine test_storage_scale

    subroutine test_storage_derivative_precision()
        !! The returned fluid content is subnormal, but its M derivative is
        !! normal. Reusing the rounded p/M intermediate amplifies its error.
        type(material_t) :: material
        real(wp) :: x(nfields,1),y(nfields,1),dsdp(nfields,1,1)
        real(wp) :: storage_modulus
        integer :: status

        storage_modulus=0.75_wp*sqrt(epsilon(1.0_wp))
        call material%init(biot_poroelastic_t(1.0_wp,1.0_wp,0.0_wp,storage_modulus),1,status, &
            sensitivity_parameters=[character(len=12) :: "biot_modulus"])
        call expect(status,cf_ok)
        x=0.0_wp
        x(nfields,1)=nearest(0.0_wp,1.0_wp)
        call material%eval(x,y,status=status,dsdp=dsdp)
        call expect(status,cf_ok)
        call close(dsdp(nfields,1,1)/tiny(1.0_wp),-16.0_wp/9.0_wp)
    end subroutine test_storage_derivative_precision

    subroutine test_requested_output_range()
        type(material_t) :: tracked,plain
        real(wp) :: x(nfields,1),y(nfields,1),tangent(nfields,nfields,1),dsdp(nfields,1,1)
        real(wp) :: storage_modulus
        integer :: status

        ! M remains admissible even if 1/M is unrepresentable: zero response
        ! and zero selected M derivative are valid until a tangent is requested.
        storage_modulus=nearest(0.0_wp,1.0_wp)
        call tracked%init(biot_poroelastic_t(1.0_wp,1.0_wp,0.0_wp,storage_modulus),1,status, &
            sensitivity_parameters=[character(len=12) :: "biot_modulus"])
        call expect(status,cf_ok)
        x=0.0_wp
        call tracked%eval(x,y,status=status,dsdp=dsdp)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(dsdp))) then
            error stop "Nonfinite zero Biot response"
        end if
        call close(maxval(abs(y)),0.0_wp)
        call close(maxval(abs(dsdp)),0.0_wp)
        call tracked%eval(x,y,tangent,status,dsdp=dsdp)
        call expect(status,cf_evaluation_failed)
        if (tracked%last_failure_point()/=1) error stop "Biot tangent overflow point"
        call tracked%commit(status)
        call expect(status,cf_no_trial)
        call tracked%eval(x,y,status=status,dsdp=dsdp)
        call expect(status,cf_ok)

        storage_modulus=tiny(1.0_wp)
        call plain%init(biot_poroelastic_t(1.0_wp,1.0_wp,0.0_wp,storage_modulus),1,status)
        call expect(status,cf_ok)
        call tracked%init(biot_poroelastic_t(1.0_wp,1.0_wp,0.0_wp,storage_modulus),1,status, &
            sensitivity_parameters=[character(len=12) :: "biot_modulus"])
        call expect(status,cf_ok)
        x(nfields,1)=1.0_wp
        call plain%eval(x,y,tangent,status)
        call expect(status,cf_ok)
        call close(y(nfields,1)*storage_modulus,1.0_wp)
        call tracked%eval(x,y,status=status,dsdp=dsdp)
        call expect(status,cf_evaluation_failed)
        call tracked%commit(status)
        call expect(status,cf_no_trial)
        x(nfields,1)=storage_modulus
        call tracked%eval(x,y,status=status,dsdp=dsdp)
        call expect(status,cf_ok)
        call close(dsdp(nfields,1,1)*storage_modulus,-1.0_wp)

        x(nfields,1)=8.0_wp
        call plain%eval(x,y,status=status)
        call expect(status,cf_evaluation_failed)
        call plain%commit(status)
        call expect(status,cf_no_trial)
        x=0.0_wp
        call plain%eval(x,y,status=status)
        call expect(status,cf_ok)
    end subroutine test_requested_output_range

    subroutine test_invalid_and_recovery()
        type(material_t) :: material
        real(wp) :: x(nfields,npoints),y(nfields,npoints)
        integer :: status

        call material%init(biot_poroelastic_t(100.0_wp,40.0_wp,1.1_wp,300.0_wp),npoints,status)
        call expect(status,cf_invalid_argument)
        call material%init(biot_poroelastic_t(100.0_wp,40.0_wp,-0.1_wp,300.0_wp),npoints,status)
        call expect(status,cf_invalid_argument)
        call material%init(biot_poroelastic_t(100.0_wp,40.0_wp,0.7_wp,0.0_wp),npoints,status)
        call expect(status,cf_invalid_argument)
        call material%init(biot_poroelastic_t(100.0_wp,40.0_wp,ieee_value(0.0_wp,ieee_quiet_nan),300.0_wp), &
            npoints,status)
        call expect(status,cf_invalid_argument)
        call material%init(biot_poroelastic_t(100.0_wp,40.0_wp,1.0_wp,300.0_wp),npoints,status)
        call expect(status,cf_ok)
        x=0.0_wp
        x(nfields,2)=ieee_value(0.0_wp,ieee_quiet_nan)
        call material%eval(x,y,status=status)
        call expect(status,cf_evaluation_failed)
        if (material%last_failure_point()/=2) error stop "Biot failed point"
        call material%commit(status)
        call expect(status,cf_no_trial)
        x(nfields,2)=-5.0_wp
        call material%eval(x,y,status=status)
        call expect(status,cf_ok)
        call close(y(1,2),5.0_wp)
        call close(y(nfields,2),-1.0_wp/60.0_wp)
    end subroutine test_invalid_and_recovery

end program test_biot_poroelastic

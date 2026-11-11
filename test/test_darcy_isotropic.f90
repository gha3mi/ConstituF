! SPDX-License-Identifier: Apache-2.0

program test_darcy_isotropic
    !! Test convention: pressure in Pa, length in m and time in s.
    !! Permeability is m2, viscosity Pa*s and Darcy velocity m/s.
    !!
    !! Darcy reference, pressure work, tangent, parameter and API regression.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
    use constituf, only: wp, material_t, cf_ok, cf_invalid_argument, &
                        cf_invalid_shape, cf_evaluation_failed, cf_no_trial, cf_unsupported
    use constituf_darcy_isotropic, only: darcy_isotropic_t
    implicit none (type, external)

    integer, parameter :: ndim = 3

    call test_reference_and_schema()
    call test_tangent_and_covariance()
    call test_parameter_derivatives()
    call test_common_parameter_scale(sqrt(tiny(1.0_wp)))
    call test_common_parameter_scale(sqrt(huge(1.0_wp)))
    call test_sensitivity_range()
    call test_subnormal_mobility()
    call test_invalid_and_recovery()

contains

    subroutine expect(actual, wanted)
        integer, intent(in) :: actual, wanted

        if (actual /= wanted) then
            print *, "Darcy status", actual, "expected", wanted
            error stop "Darcy status mismatch"
        end if
    end subroutine expect

    subroutine close(actual, wanted)
        real(wp), intent(in) :: actual, wanted

        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(wanted)) then
            error stop "Nonfinite Darcy comparison"
        end if
        if (abs(actual-wanted) > 256.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(wanted))) then
            print *, "Darcy actual", actual, "expected", wanted
            error stop "Darcy reference mismatch"
        end if
    end subroutine close

    subroutine test_reference_and_schema()
        type(material_t) :: material
        real(wp) :: gradient(ndim,3), velocity(ndim,3), tangent(ndim,ndim,3)
        integer :: status, ipoint, i, j

        call material%init(darcy_isotropic_t(3.0_wp,2.0_wp),3,status)
        call expect(status,cf_ok)
        if (material%input_size() /= ndim .or. material%output_size() /= ndim) then
            error stop "Darcy field extent mismatch"
        end if
        if (material%state_storage_bytes() /= 0 .or. material%state_size() /= 0) then
            error stop "Darcy is stateless"
        end if
        if (.not. material%field_matches("input", "pressure.gradient", &
            "vector", "pressure_gradient", "spatial", [ndim])) then
            error stop "Darcy input field metadata"
        end if
        if (.not. material%field_matches("output", "fluid.darcy_velocity", &
            "vector", "specific_discharge", "spatial", [ndim])) then
            error stop "Darcy velocity metadata"
        end if
        if (material%field_matches("output", "fluid.darcy_velocity", &
            "vector", "specific_discharge", "reference", [ndim])) then
            error stop "Darcy invalid coordinate frame accepted"
        end if
        if (material%field_matches("input", "pressure.gradient", &
            "vector", "hydraulic_head_gradient", "spatial", [ndim])) then
            error stop "Darcy incompatible gradient measure accepted"
        end if

        gradient(:,1) = [2.0_wp,-4.0_wp,1.0_wp]
        gradient(:,2) = [-1.0_wp,3.0_wp,2.0_wp]
        gradient(:,3) = 0.0_wp
        call material%eval(gradient,velocity,tangent,status)
        call expect(status,cf_ok)
        call close(velocity(1,1),-3.0_wp)
        call close(velocity(2,1),6.0_wp)
        call close(velocity(3,1),-1.5_wp)
        call close(-dot_product(velocity(:,1),gradient(:,1)),63.0_wp/2.0_wp)
        do ipoint = 1, 3
            do i = 1, ndim
                call close(velocity(i,ipoint),-1.5_wp*gradient(i,ipoint))
                do j = 1, ndim
                    if (i == j) then
                        call close(tangent(i,j,ipoint),-1.5_wp)
                    else
                        call close(tangent(i,j,ipoint),0.0_wp)
                    end if
                end do
            end do
            if (-dot_product(velocity(:,ipoint),gradient(:,ipoint)) < 0.0_wp) then
                error stop "Darcy pressure work violates dissipation inequality"
            end if
        end do
        call material%commit(status)
        call expect(status,cf_ok)
        call material%eval(gradient,velocity,status=status)
        call expect(status,cf_ok)
        call material%revert()
        call material%commit(status)
        call expect(status,cf_no_trial)
    end subroutine test_reference_and_schema

    subroutine test_tangent_and_covariance()
        type(material_t) :: batch, single
        real(wp) :: gradient(ndim,2), velocity(ndim,2), tangent(ndim,ndim,2)
        real(wp) :: changed(ndim,2), plus(ndim,2), minus(ndim,2), alone(ndim,1)
        real(wp) :: h, derivative, error_bound
        integer :: status, ipoint, i, j

        call batch%init(darcy_isotropic_t(8.0_wp,5.0_wp),2,status)
        call expect(status,cf_ok)
        call single%init(darcy_isotropic_t(8.0_wp,5.0_wp),1,status)
        call expect(status,cf_ok)
        gradient(:,1) = [1.3_wp,-0.8_wp,2.1_wp]
        gradient(:,2) = [-gradient(2,1),gradient(1,1),gradient(3,1)]
        call batch%eval(gradient,velocity,tangent,status)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(tangent))) error stop "Nonfinite Darcy analytical tangent"
        call single%eval(gradient(:,1:1),alone,status=status)
        call expect(status,cf_ok)
        do i = 1, ndim
            call close(velocity(i,1),alone(i,1))
        end do
        call close(velocity(1,2),-velocity(2,1))
        call close(velocity(2,2),velocity(1,1))
        call close(velocity(3,2),velocity(3,1))

        do ipoint = 1, 2
            do j = 1, ndim
                ! The gradient enters affinely, so a finite central step
                ! avoids cancellation without introducing truncation error.
                h = 0.25_wp*maxval(abs(gradient(:,ipoint)))
                changed = gradient
                changed(j,ipoint) = gradient(j,ipoint)+h
                call batch%eval(changed,plus,status=status)
                call expect(status,cf_ok)
                changed(j,ipoint) = gradient(j,ipoint)-h
                call batch%eval(changed,minus,status=status)
                call expect(status,cf_ok)
                if (.not. all(ieee_is_finite(plus)) .or. .not. all(ieee_is_finite(minus))) then
                    error stop "Nonfinite Darcy tangent secant response"
                end if
                error_bound = 32.0_wp*epsilon(1.0_wp)* &
                    maxval(abs(plus(:,ipoint))+abs(minus(:,ipoint)))/(2.0_wp*h)
                if (.not. ieee_is_finite(error_bound)) error stop "Nonfinite Darcy tangent error bound"
                do i = 1, ndim
                    derivative = (plus(i,ipoint)-minus(i,ipoint))/(2.0_wp*h)
                    if (.not. ieee_is_finite(derivative) .or. &
                        abs(derivative-tangent(i,j,ipoint)) > error_bound) then
                        error stop "Darcy tangent finite difference"
                    end if
                end do
            end do
        end do
    end subroutine test_tangent_and_covariance

    subroutine test_parameter_derivatives()
        type(material_t) :: tracked, upper, lower
        real(wp) :: gradient(ndim,1), velocity(ndim,1), dsdp(ndim,2,1)
        real(wp) :: plus(ndim,1), minus(ndim,1), wrong(ndim,1,1)
        real(wp) :: derivative, h, bound, relative_step
        integer :: status, iparam, i

        call tracked%init(darcy_isotropic_t(3.0_wp,2.0_wp),1,status, &
            sensitivity_parameters=[character(len=17) :: "dynamic_viscosity","permeability"])
        call expect(status,cf_ok)
        if (tracked%sensitivity_size() /= 2) error stop "Darcy parameter count"
        if (tracked%state_storage_bytes() /= 0) error stop "Darcy sensitivity history"
        gradient(:,1) = [1.2_wp,-2.1_wp,0.4_wp]
        call tracked%eval(gradient,velocity,status=status,dsdp=dsdp)
        call expect(status,cf_ok)
        do i = 1, ndim
            call close(dsdp(i,1,1),0.75_wp*gradient(i,1))
            call close(dsdp(i,2,1),-0.5_wp*gradient(i,1))
        end do

        do iparam = 1, 2
            select case (iparam)
            case (1)
                h=2.0_wp*epsilon(1.0_wp)**(1.0_wp/3.0_wp)
                call upper%init(darcy_isotropic_t(3.0_wp,2.0_wp+h),1,status)
                call expect(status,cf_ok)
                call lower%init(darcy_isotropic_t(3.0_wp,2.0_wp-h),1,status)
            case (2)
                ! Permeability enters the velocity affinely.
                h=0.25_wp*3.0_wp
                call upper%init(darcy_isotropic_t(3.0_wp+h,2.0_wp),1,status)
                call expect(status,cf_ok)
                call lower%init(darcy_isotropic_t(3.0_wp-h,2.0_wp),1,status)
            case default
                error stop "Darcy invalid parameter selection"
            end select
            call expect(status,cf_ok)
            call upper%eval(gradient,plus,status=status)
            call expect(status,cf_ok)
            call lower%eval(gradient,minus,status=status)
            call expect(status,cf_ok)
            if (.not. all(ieee_is_finite(plus)) .or. .not. all(ieee_is_finite(minus))) then
                error stop "Nonfinite Darcy parameter secant response"
            end if
            do i = 1, ndim
                derivative = (plus(i,1)-minus(i,1))/(2.0_wp*h)
                bound=32.0_wp*epsilon(1.0_wp)*(abs(plus(i,1))+abs(minus(i,1)))/(2.0_wp*h)
                if (iparam==1) then
                    ! Exact central-difference bias for inverse viscosity.
                    relative_step=h/2.0_wp
                    bound=bound+abs(dsdp(i,iparam,1))*relative_step**2/(1.0_wp-relative_step**2)
                end if
                if (.not. all(ieee_is_finite([derivative,bound]))) then
                    error stop "Nonfinite Darcy parameter comparison"
                end if
                if (abs(derivative-dsdp(i,iparam,1))>bound) then
                    error stop "Darcy parameter derivative finite difference"
                end if
            end do
        end do
        call tracked%commit(status)
        call expect(status,cf_ok)
        call tracked%eval(gradient,velocity,status=status)
        call expect(status,cf_invalid_argument)
        call tracked%eval(gradient,velocity,status=status,dsdp=wrong)
        call expect(status,cf_invalid_shape)
        call tracked%eval(gradient,velocity,status=status,dsdp=dsdp)
        call expect(status,cf_ok)
    end subroutine test_parameter_derivatives

    subroutine test_common_parameter_scale(parameter_scale)
        !! Scaling permeability and viscosity together leaves the mobility
        !! unchanged; the two physical parameter derivatives scale inversely.
        real(wp), intent(in) :: parameter_scale
        type(material_t) :: material
        real(wp) :: gradient(ndim,1),velocity(ndim,1),tangent(ndim,ndim,1)
        real(wp) :: dsdp(ndim,2,1)
        integer :: status,i,j

        call material%init(darcy_isotropic_t(3.0_wp*parameter_scale,2.0_wp*parameter_scale),1,status, &
            sensitivity_parameters=[character(len=17) :: "dynamic_viscosity","permeability"])
        call expect(status,cf_ok)
        gradient(:,1)=[2.0_wp,-4.0_wp,1.0_wp]
        call material%eval(gradient,velocity,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        do i=1,ndim
            call close(velocity(i,1),-1.5_wp*gradient(i,1))
            call close(dsdp(i,1,1)*parameter_scale,0.75_wp*gradient(i,1))
            call close(dsdp(i,2,1)*parameter_scale,-0.5_wp*gradient(i,1))
            do j=1,ndim
                if (i==j) then
                    call close(tangent(i,j,1),-1.5_wp)
                else
                    call close(tangent(i,j,1),0.0_wp)
                end if
            end do
        end do
    end subroutine test_common_parameter_scale

    subroutine test_sensitivity_range()
        type(material_t) :: tracked,plain
        real(wp), parameter :: direction(ndim)=[1.0_wp,-0.5_wp,0.0_wp]
        real(wp) :: gradient(ndim,1),velocity(ndim,1),tangent(ndim,ndim,1)
        real(wp) :: dsdp(ndim,2,1)
        real(wp) :: permeability,viscosity
        integer :: status,i,j

        ! Velocity legitimately rounds to zero, but differentiation of the
        ! original constitutive factors gives a finite nonzero viscosity
        ! derivative. Dividing the rounded velocity loses that derivative.
        permeability=tiny(1.0_wp)
        viscosity=sqrt(permeability)
        call tracked%init(darcy_isotropic_t(permeability,viscosity),1,status, &
            sensitivity_parameters=[character(len=17) :: "dynamic_viscosity","permeability"])
        call expect(status,cf_ok)
        gradient(:,1)=permeability*direction
        call tracked%eval(gradient,velocity,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        do i=1,ndim
            call close(velocity(i,1),0.0_wp)
            call close(dsdp(i,1,1)/permeability,direction(i))
            call close(dsdp(i,2,1)/viscosity,-direction(i))
            do j=1,ndim
                if (i==j) then
                    call close(tangent(i,j,1)/viscosity,-1.0_wp)
                else
                    call close(tangent(i,j,1)/viscosity,0.0_wp)
                end if
            end do
        end do

        ! A finite response does not make an unrepresentable requested
        ! parameter derivative valid; the failed trial must not be committed.
        viscosity=permeability
        call plain%init(darcy_isotropic_t(permeability,viscosity),1,status)
        call expect(status,cf_ok)
        call tracked%init(darcy_isotropic_t(permeability,viscosity),1,status, &
            sensitivity_parameters=[character(len=17) :: "dynamic_viscosity","permeability"])
        call expect(status,cf_ok)
        gradient(:,1)=[8.0_wp,0.0_wp,0.0_wp]
        call plain%eval(gradient,velocity,status=status)
        call expect(status,cf_ok)
        call close(velocity(1,1),-8.0_wp)
        call tracked%eval(gradient,velocity,status=status,dsdp=dsdp)
        call expect(status,cf_evaluation_failed)
        if (tracked%last_failure_point()/=1) error stop "Darcy derivative overflow point"
        call tracked%commit(status)
        call expect(status,cf_no_trial)
        gradient(1,1)=1.0_wp
        call tracked%eval(gradient,velocity,status=status,dsdp=dsdp)
        call expect(status,cf_ok)
        call close(dsdp(1,1,1)*viscosity,1.0_wp)
        call close(dsdp(1,2,1)*viscosity,-1.0_wp)
    end subroutine test_sensitivity_range

    subroutine test_subnormal_mobility()
        !! The tangent rounds the mobility itself. Velocity must retain
        !! the original quotient through multiplication by the gradient.
        type(material_t) :: material
        integer, parameter :: ncases=2
        real(wp), parameter :: permeability_factors(ncases)=[1.0_wp,5.0_wp]
        real(wp), parameter :: viscosities(ncases)=[0.75_wp,4.0_wp]
        real(wp), parameter :: gradients(ncases)=[3.0_wp,4.0_wp]
        real(wp), parameter :: expected_velocity(ncases)=[-4.0_wp,-5.0_wp]
        real(wp), parameter :: viscosity_derivative(ncases)=[5.0_wp,1.0_wp]
        real(wp), parameter :: permeability_derivative(ncases)=[-4.0_wp,-1.0_wp]
        real(wp) :: gradient(ndim,1),velocity(ndim,1),tangent(ndim,ndim,1)
        real(wp) :: dsdp(ndim,2,1)
        real(wp) :: smallest
        integer :: status,icase,i,j

        smallest=nearest(0.0_wp,1.0_wp)
        do icase=1,ncases
            call material%init(darcy_isotropic_t(permeability_factors(icase)*smallest,viscosities(icase)), &
                1,status,sensitivity_parameters=[character(len=17) :: "dynamic_viscosity","permeability"])
            call expect(status,cf_ok)
            gradient(:,1)=[gradients(icase),0.0_wp,0.0_wp]
            call material%eval(gradient,velocity,tangent,status,dsdp=dsdp)
            call expect(status,cf_ok)
            call close(velocity(1,1)/smallest,expected_velocity(icase))
            call close(dsdp(1,1,1)/smallest,viscosity_derivative(icase))
            call close(dsdp(1,2,1),permeability_derivative(icase))
            do i=2,ndim
                call close(velocity(i,1),0.0_wp)
                call close(dsdp(i,1,1),0.0_wp)
                call close(dsdp(i,2,1),0.0_wp)
            end do
            do j=1,ndim
                do i=1,ndim
                    call close(tangent(i,j,1)/smallest,merge(-1.0_wp,0.0_wp,i==j))
                end do
            end do
        end do
    end subroutine test_subnormal_mobility

    subroutine test_invalid_and_recovery()
        type(material_t) :: material
        real(wp) :: gradient(ndim,2), velocity(ndim,2), wrong(ndim-1,2)
        real(wp) :: extra(ndim,1,2)
        integer :: status

        call material%init(darcy_isotropic_t(0.0_wp,1.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(darcy_isotropic_t(1.0_wp,0.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(darcy_isotropic_t(-1.0_wp,1.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(darcy_isotropic_t(1.0_wp,-1.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(darcy_isotropic_t(huge(1.0_wp)/4.0_wp,0.125_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(darcy_isotropic_t(ieee_value(0.0_wp,ieee_quiet_nan), &
                           1.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(darcy_isotropic_t(2.0_wp,1.0_wp),2,status)
        call expect(status,cf_ok)

        gradient = 0.0_wp
        gradient(1,:) = 1.0_wp
        gradient(1,2) = ieee_value(0.0_wp,ieee_positive_inf)
        call material%eval(gradient,velocity,status=status)
        call expect(status,cf_evaluation_failed)
        if (material%last_failure_point() /= 2) error stop "Darcy failing point"
        call material%commit(status)
        call expect(status,cf_no_trial)
        gradient(1,2) = ieee_value(0.0_wp,ieee_quiet_nan)
        call material%eval(gradient,velocity,status=status)
        call expect(status,cf_evaluation_failed)
        gradient(1,2) = huge(1.0_wp)
        call material%eval(gradient,velocity,status=status)
        call expect(status,cf_evaluation_failed)
        gradient(1,2) = 2.0_wp
        call material%eval(gradient,velocity,status=status)
        call expect(status,cf_ok)
        call close(velocity(1,2),-4.0_wp)
        call material%eval(wrong,velocity,status=status)
        call expect(status,cf_invalid_shape)
        call material%eval(gradient,velocity,status=status,dsdp=extra)
        call expect(status,cf_unsupported)
    end subroutine test_invalid_and_recovery

end program test_darcy_isotropic

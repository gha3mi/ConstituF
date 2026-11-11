! SPDX-License-Identifier: Apache-2.0

program test_darcy_exponential
    !! Test convention: pressure in Pa, length in m and time in s.
    !! Permeability is m2, viscosity Pa*s and Darcy velocity m/s.
    !!
    !! Pressure-dependent Darcy flux, 3x4 tangent and four local sensitivities.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
    use constituf, only: wp, material_t, cf_ok, cf_invalid_argument, &
                        cf_invalid_shape, cf_evaluation_failed, cf_no_trial, cf_unsupported
    use constituf_darcy_exponential, only: darcy_exponential_t
    use constituf_darcy_isotropic, only: darcy_isotropic_t
    implicit none (type, external)

    integer, parameter :: ndim = 3, ninput = ndim+1, nparameters = 4
    character(len=23), parameter :: parameter_names(nparameters) = [character(len=23) :: &
        "reference_pressure","pressure_coefficient","dynamic_viscosity","reference_permeability"]

    call test_reference_and_schema()
    call test_tangent_and_rotation()
    call test_sensitivities()
    call test_exact_pressure_reference()
    call test_reference_mobility(tiny(1.0_wp))
    call test_reference_mobility(huge(1.0_wp)/4.0_wp)
    call test_reference_ratio_identity()
    call test_rounded_velocity_derivatives()
    call test_normal_mobility_domain()
    call test_requested_derivative_range()
    call test_invalid_and_recovery()

contains

    function model_from_values(values) result(model)
        real(wp), intent(in) :: values(nparameters)
        type(darcy_exponential_t) :: model

        model = darcy_exponential_t(values(1),values(2),values(3),values(4))
    end function model_from_values

    subroutine expect(actual, wanted)
        integer, intent(in) :: actual, wanted

        if (actual /= wanted) then
            print *, "Status", actual, "expected", wanted
            error stop "Exponential Darcy status mismatch"
        end if
    end subroutine expect

    subroutine close(actual, wanted)
        real(wp), intent(in) :: actual, wanted

        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(wanted) .or. abs(actual-wanted) > &
            512.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(wanted))) then
            print *, "Actual", actual, "expected", wanted
            error stop "Exponential Darcy reference mismatch"
        end if
    end subroutine close

    subroutine test_reference_and_schema()
        type(material_t) :: material, constant
        real(wp) :: x(ninput,3), velocity(ndim,3), tangent(ndim,ninput,3)
        real(wp) :: baseline(ndim,3)
        real(wp) :: mobility
        integer :: status, ipoint, i, j
        integer :: empty_shape(0)

        call material%init(darcy_exponential_t(3.0_wp,2.0_wp,0.1_wp,10.0_wp),3,status)
        call expect(status,cf_ok)
        if (material%input_size() /= ninput .or. material%output_size() /= ndim) then
            error stop "Exponential Darcy field sizes"
        end if
        if (material%state_storage_bytes() /= 0 .or. material%state_size() /= 0) then
            error stop "Exponential Darcy must be stateless"
        end if
        if (.not. material%field_matches("input", "pressure.gradient", &
            "vector", "pressure_gradient", "spatial", [ndim])) then
            error stop "Exponential Darcy gradient metadata"
        end if
        if (.not. material%field_matches("input", "pressure.absolute", &
            "scalar", "absolute_pressure", "material", empty_shape)) then
            error stop "Exponential Darcy absolute pressure metadata"
        end if
        if (.not. material%field_matches("output", "fluid.darcy_velocity", &
            "vector", "specific_discharge", "spatial", [ndim])) then
            error stop "Exponential Darcy velocity metadata"
        end if
        if (material%field_matches("input", "pressure.absolute", &
            "scalar", "absolute_pressure", "spatial", empty_shape)) then
            error stop "Incompatible scalar pressure frame accepted"
        end if

        x(:,1) = [1.0_wp,-2.0_wp,0.5_wp,12.0_wp]
        x(:,2) = [0.0_wp,0.0_wp,0.0_wp,10.0_wp]
        x(:,3) = [-0.25_wp,0.1_wp,1.2_wp,8.0_wp]
        call material%eval(x,velocity,tangent,status)
        call expect(status,cf_ok)
        do ipoint = 1, 3
            mobility = 1.5_wp*exp(0.1_wp*(x(ninput,ipoint)-10.0_wp))
            do i = 1, ndim
                call close(velocity(i,ipoint),-mobility*x(i,ipoint))
                call close(tangent(i,ninput,ipoint),0.1_wp*velocity(i,ipoint))
                do j = 1, ndim
                    if (i == j) then
                        call close(tangent(i,j,ipoint),-mobility)
                    else
                        call close(tangent(i,j,ipoint),0.0_wp)
                    end if
                end do
            end do
            if (-dot_product(x(1:ndim,ipoint),velocity(:,ipoint)) < 0.0_wp) then
                error stop "Exponential Darcy negative dissipation"
            end if
        end do
        call material%commit(status)
        call expect(status,cf_ok)
        call material%eval(x,velocity,status=status)
        call expect(status,cf_ok)
        call material%revert()
        call material%commit(status)
        call expect(status,cf_no_trial)

        call material%init(darcy_exponential_t(3.0_wp,2.0_wp,0.0_wp,10.0_wp),3,status)
        call expect(status,cf_ok)
        call constant%init(darcy_isotropic_t(3.0_wp,2.0_wp),3,status)
        call expect(status,cf_ok)
        call material%eval(x,velocity,tangent,status)
        call expect(status,cf_ok)
        call constant%eval(x(1:ndim,:),baseline,status=status)
        call expect(status,cf_ok)
        do ipoint = 1, 3
            do i = 1, ndim
                call close(velocity(i,ipoint),baseline(i,ipoint))
                call close(tangent(i,ninput,ipoint),0.0_wp)
            end do
        end do
    end subroutine test_reference_and_schema

    subroutine test_tangent_and_rotation()
        type(material_t) :: batch, single
        real(wp) :: x(ninput,2), velocity(ndim,2), tangent(ndim,ninput,2)
        real(wp) :: changed(ninput,2), plus(ndim,2), minus(ndim,2), alone(ndim,1)
        real(wp) :: h, derivative, bound, exponent_step
        integer :: status, ipoint, i, j

        call batch%init(darcy_exponential_t(3.0_wp,2.0_wp,0.12_wp,10.0_wp),2,status)
        call expect(status,cf_ok)
        call single%init(darcy_exponential_t(3.0_wp,2.0_wp,0.12_wp,10.0_wp),1,status)
        call expect(status,cf_ok)
        x(:,1) = [1.3_wp,-0.9_wp,2.1_wp,12.0_wp]
        x(:,2) = [-x(2,1),x(1,1),x(3,1),x(ninput,1)]
        call batch%eval(x,velocity,tangent,status)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(tangent))) error stop "Nonfinite exponential Darcy tangent"
        call single%eval(x(:,1:1),alone,status=status)
        call expect(status,cf_ok)
        do i = 1, ndim
            call close(velocity(i,1),alone(i,1))
        end do
        call close(velocity(1,2),-velocity(2,1))
        call close(velocity(2,2),velocity(1,1))
        call close(velocity(3,2),velocity(3,1))

        do ipoint = 1, 2
            do j = 1, ninput
                if (j==ninput) then
                    h=epsilon(1.0_wp)**(1.0_wp/3.0_wp)/0.12_wp
                else
                    ! The gradient enters affinely, without truncation error.
                    h=0.25_wp*maxval(abs(x(1:ndim,ipoint)))
                end if
                changed = x
                changed(j,ipoint) = x(j,ipoint)+h
                call batch%eval(changed,plus,status=status)
                call expect(status,cf_ok)
                if (.not. all(ieee_is_finite(plus))) error stop "Nonfinite positive input perturbation"
                changed(j,ipoint) = x(j,ipoint)-h
                call batch%eval(changed,minus,status=status)
                call expect(status,cf_ok)
                if (.not. all(ieee_is_finite(minus))) error stop "Nonfinite negative input perturbation"
                do i = 1, ndim
                    derivative = (plus(i,ipoint)-minus(i,ipoint))/(2.0_wp*h)
                    bound=32.0_wp*epsilon(1.0_wp)*(abs(plus(i,ipoint))+abs(minus(i,ipoint)))/(2.0_wp*h)
                    if (j==ninput) then
                        exponent_step=0.12_wp*h
                        bound=bound+abs(tangent(i,j,ipoint))*exponent_step**2*exp(exponent_step)/6.0_wp
                    end if
                    if (.not. ieee_is_finite(derivative) .or. .not. ieee_is_finite(bound) .or. &
                        abs(derivative-tangent(i,j,ipoint))>bound) then
                        error stop "Exponential Darcy 3x4 tangent finite difference"
                    end if
                end do
            end do
        end do
    end subroutine test_tangent_and_rotation

    subroutine test_sensitivities()
        type(material_t) :: tracked, upper, lower
        real(wp) :: values(nparameters), high(nparameters), low(nparameters)
        real(wp) :: x(ninput,1), velocity(ndim,1), dsdp(ndim,nparameters,1)
        real(wp) :: plus(ndim,1), minus(ndim,1), bad(ndim,nparameters-1,1)
        real(wp) :: derivative, h, bound, exponent_step, relative_step
        integer, parameter :: ids(nparameters) = [4,3,2,1]
        integer :: status, iparam, i, id

        values = [3.0_wp,2.0_wp,0.1_wp,10.0_wp]
        call tracked%init(model_from_values(values),1,status, &
            sensitivity_parameters=parameter_names)
        call expect(status,cf_ok)
        if (tracked%sensitivity_size() /= nparameters) error stop "Sensitivity count"
        if (tracked%state_storage_bytes() /= 0) error stop "Derivative history"
        x(:,1) = [1.2_wp,-0.9_wp,2.1_wp,12.0_wp]
        call tracked%eval(x,velocity,status=status,dsdp=dsdp)
        call expect(status,cf_ok)
        do i = 1, ndim
            call close(dsdp(i,1,1),-0.1_wp*velocity(i,1))
            call close(dsdp(i,2,1),2.0_wp*velocity(i,1))
            call close(dsdp(i,3,1),-velocity(i,1)/2.0_wp)
            call close(dsdp(i,4,1),velocity(i,1)/3.0_wp)
        end do

        do iparam = 1, nparameters
            id = ids(iparam)
            exponent_step=0.0_wp
            relative_step=0.0_wp
            select case (id)
            case (1)
                h=0.25_wp*values(1)
            case (2)
                relative_step=epsilon(1.0_wp)**(1.0_wp/3.0_wp)
                h=relative_step*values(2)
            case (3)
                exponent_step=epsilon(1.0_wp)**(1.0_wp/3.0_wp)
                h=exponent_step/abs(x(ninput,1)-values(4))
            case (4)
                exponent_step=epsilon(1.0_wp)**(1.0_wp/3.0_wp)
                h=exponent_step/abs(values(3))
            case default
                error stop "Invalid exponential Darcy parameter index"
            end select
            high = values
            low = values
            high(id) = high(id)+h
            low(id) = low(id)-h
            call upper%init(model_from_values(high),1,status)
            call expect(status,cf_ok)
            call lower%init(model_from_values(low),1,status)
            call expect(status,cf_ok)
            call upper%eval(x,plus,status=status)
            call expect(status,cf_ok)
            if (.not. all(ieee_is_finite(plus))) error stop "Nonfinite positive parameter perturbation"
            call lower%eval(x,minus,status=status)
            call expect(status,cf_ok)
            if (.not. all(ieee_is_finite(minus))) error stop "Nonfinite negative parameter perturbation"
            do i = 1, ndim
                derivative = (plus(i,1)-minus(i,1))/(2.0_wp*h)
                ! Exponential columns have a bounded sinh(z)/z central
                ! bias, viscosity has an exact reciprocal bias, and k0 is affine.
                bound=32.0_wp*epsilon(1.0_wp)*(abs(plus(i,1))+abs(minus(i,1)))/(2.0_wp*h)+ &
                    abs(dsdp(i,iparam,1))*(exponent_step**2*exp(exponent_step)/6.0_wp+ &
                    relative_step**2/(1.0_wp-relative_step**2))
                if (.not. ieee_is_finite(derivative) .or. .not. ieee_is_finite(bound) .or. &
                    abs(derivative-dsdp(i,iparam,1))>bound) then
                    error stop "Exponential Darcy parameter derivative"
                end if
            end do
        end do
        call tracked%commit(status)
        call expect(status,cf_ok)
        call tracked%eval(x,velocity,status=status)
        call expect(status,cf_invalid_argument)
        call tracked%eval(x,velocity,status=status,dsdp=bad)
        call expect(status,cf_invalid_shape)
        call tracked%eval(x,velocity,status=status,dsdp=dsdp)
        call expect(status,cf_ok)
    end subroutine test_sensitivities

    subroutine test_exact_pressure_reference()
        type(material_t) :: material
        real(wp), parameter :: direction(ndim)=[1.0_wp,-2.0_wp,0.5_wp]
        real(wp), parameter :: pressure_offsets(3)=[2.0_wp,0.0_wp,-2.0_wp]
        real(wp) :: x(ninput,3),velocity(ndim,3),tangent(ndim,ninput,3)
        real(wp) :: dsdp(ndim,nparameters,3),mobilities(3),reference(ndim)
        real(wp) :: beta
        integer :: status,icase,ipoint,i,j

        do icase=1,2
            if (icase==1) then
                beta=log(2.0_wp)/2.0_wp
                mobilities=[3.0_wp,1.5_wp,0.75_wp]
            else
                beta=-log(2.0_wp)/2.0_wp
                mobilities=[0.75_wp,1.5_wp,3.0_wp]
            end if
            call material%init(darcy_exponential_t(3.0_wp,2.0_wp,beta,10.0_wp),3,status, &
                sensitivity_parameters=parameter_names)
            call expect(status,cf_ok)
            do ipoint=1,3
                x(1:ndim,ipoint)=direction
                x(ninput,ipoint)=10.0_wp+pressure_offsets(ipoint)
            end do
            call material%eval(x,velocity,tangent,status,dsdp=dsdp)
            call expect(status,cf_ok)
            do ipoint=1,3
                reference=-mobilities(ipoint)*direction
                do i=1,ndim
                    call close(velocity(i,ipoint),reference(i))
                    call close(tangent(i,ninput,ipoint),beta*reference(i))
                    call close(dsdp(i,1,ipoint),-beta*reference(i))
                    call close(dsdp(i,2,ipoint),pressure_offsets(ipoint)*reference(i))
                    call close(dsdp(i,3,ipoint),-0.5_wp*reference(i))
                    call close(dsdp(i,4,ipoint),reference(i)/3.0_wp)
                    do j=1,ndim
                        if (i==j) then
                            call close(tangent(i,j,ipoint),-mobilities(ipoint))
                        else
                            call close(tangent(i,j,ipoint),0.0_wp)
                        end if
                    end do
                end do
                call close(-dot_product(direction,velocity(:,ipoint)),(21.0_wp/4.0_wp)*mobilities(ipoint))
            end do
        end do
    end subroutine test_exact_pressure_reference

    subroutine test_reference_mobility(mobility)
        !! The documented normal-mobility domain includes its actual lower
        !! boundary and large finite values; arbitrary logarithmic margins do not.
        real(wp), intent(in) :: mobility
        type(material_t) :: material
        real(wp), parameter :: direction(ndim)=[1.0_wp,-0.5_wp,0.25_wp]
        real(wp) :: x(ninput,1),velocity(ndim,1),tangent(ndim,ninput,1),dsdp(ndim,nparameters,1)
        integer :: status,i,j

        call material%init(darcy_exponential_t(mobility,1.0_wp,0.0_wp,1.0_wp),1,status, &
            sensitivity_parameters=parameter_names)
        call expect(status,cf_ok)
        x(1:ndim,1)=direction
        x(ninput,1)=1.0_wp
        call material%eval(x,velocity,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        do i=1,ndim
            call close(velocity(i,1)/mobility,-direction(i))
            call close(tangent(i,ninput,1),0.0_wp)
            call close(dsdp(i,1,1),0.0_wp)
            call close(dsdp(i,2,1),0.0_wp)
            call close(dsdp(i,3,1)/mobility,direction(i))
            call close(dsdp(i,4,1),-direction(i))
            do j=1,ndim
                if (i==j) then
                    call close(tangent(i,j,1)/mobility,-1.0_wp)
                else
                    call close(tangent(i,j,1)/mobility,0.0_wp)
                end if
            end do
        end do
    end subroutine test_reference_mobility

    subroutine test_reference_ratio_identity()
        type(material_t) :: material
        real(wp), parameter :: direction(ndim)=[1.0_wp,-0.5_wp,0.25_wp]
        real(wp) :: x(ninput,1),velocity(ndim,1),tangent(ndim,ninput,1),dsdp(ndim,nparameters,1)
        real(wp) :: permeability,viscosity
        integer :: status,i

        permeability=huge(1.0_wp)/8.0_wp
        viscosity=huge(1.0_wp)/16.0_wp
        call material%init(darcy_exponential_t(permeability,viscosity,0.0_wp,1.0_wp),1,status, &
            sensitivity_parameters=parameter_names)
        call expect(status,cf_ok)
        x(1:ndim,1)=direction
        x(ninput,1)=3.0_wp
        call material%eval(x,velocity,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(velocity)) .or. &
            maxval(abs(velocity(:,1)+2.0_wp*direction))>8.0_wp*epsilon(1.0_wp)) then
            error stop "Reference mobility lost through logarithm cancellation"
        end if
        do i=1,ndim
            call close(tangent(i,i,1),-2.0_wp)
            call close(tangent(i,ninput,1),0.0_wp)
            call close(dsdp(i,1,1),0.0_wp)
            ! Beta=0 does not erase its parameter derivative when p differs from p0.
            call close(dsdp(i,2,1),-4.0_wp*direction(i))
            call close(dsdp(i,3,1)*viscosity,2.0_wp*direction(i))
            call close(dsdp(i,4,1)*permeability,-2.0_wp*direction(i))
        end do
    end subroutine test_reference_ratio_identity

    subroutine test_rounded_velocity_derivatives()
        type(material_t) :: material
        real(wp), parameter :: direction(ndim)=[1.0_wp,-0.5_wp,0.0_wp]
        real(wp) :: x(ninput,1),velocity(ndim,1),tangent(ndim,ninput,1),dsdp(ndim,nparameters,1)
        real(wp) :: small,viscosity,beta,reference_pressure,pressure
        integer :: status,i,icase

        small=tiny(1.0_wp)
        viscosity=sqrt(small)
        call material%init(darcy_exponential_t(small,viscosity,0.0_wp,1.0_wp),1,status, &
            sensitivity_parameters=parameter_names)
        call expect(status,cf_ok)
        x(1:ndim,1)=small*direction
        x(ninput,1)=1.0_wp
        call material%eval(x,velocity,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(velocity))) error stop "Nonfinite rounded Darcy velocity"
        call close(maxval(abs(velocity)),0.0_wp)
        do i=1,ndim
            call close(dsdp(i,3,1)/small,direction(i))
            call close(dsdp(i,4,1)/viscosity,-direction(i))
            call close(tangent(i,i,1)/viscosity,-1.0_wp)
        end do

        do icase=1,2
            if (icase==1) then
                beta=1.0_wp/small
                reference_pressure=1.0_wp
                pressure=1.0_wp
            else
                beta=0.0_wp
                reference_pressure=0.0_wp
                pressure=1.0_wp/small
            end if
            call material%init(darcy_exponential_t(8.0_wp*small,1.0_wp,beta,reference_pressure),1,status, &
                sensitivity_parameters=parameter_names)
            call expect(status,cf_ok)
            x(ninput,1)=pressure
            call material%eval(x,velocity,tangent,status,dsdp=dsdp)
            call expect(status,cf_ok)
            if (.not. all(ieee_is_finite(velocity))) error stop "Nonfinite rounded Darcy velocity"
            call close(maxval(abs(velocity)),0.0_wp)
            do i=1,ndim
                call close(tangent(i,i,1)/small,-8.0_wp)
                call close(dsdp(i,3,1),0.0_wp)
                call close(dsdp(i,4,1)/small,-direction(i))
                if (icase==1) then
                    call close(tangent(i,ninput,1)/small,-8.0_wp*direction(i))
                    call close(dsdp(i,1,1)/small,8.0_wp*direction(i))
                    call close(dsdp(i,2,1),0.0_wp)
                else
                    call close(tangent(i,ninput,1),0.0_wp)
                    call close(dsdp(i,1,1),0.0_wp)
                    call close(dsdp(i,2,1)/small,-8.0_wp*direction(i))
                end if
            end do
        end do
    end subroutine test_rounded_velocity_derivatives

    subroutine test_normal_mobility_domain()
        type(material_t) :: material
        real(wp) :: x(ninput,1),velocity(ndim,1)
        integer :: status

        call material%init(darcy_exponential_t(tiny(1.0_wp)/2.0_wp,1.0_wp,0.0_wp,0.0_wp),1,status)
        call expect(status,cf_invalid_argument)
        call material%init(darcy_exponential_t(tiny(1.0_wp),1.0_wp,-log(2.0_wp),0.0_wp),1,status)
        call expect(status,cf_ok)
        x(:,1)=[1.0_wp,0.0_wp,0.0_wp,1.0_wp]
        call material%eval(x,velocity,status=status)
        call expect(status,cf_evaluation_failed)
        if (material%last_failure_point()/=1) error stop "Subnormal mobility failure point"
        call material%commit(status)
        call expect(status,cf_no_trial)
        x(ninput,1)=0.0_wp
        call material%eval(x,velocity,status=status)
        call expect(status,cf_ok)
        call close(velocity(1,1)/tiny(1.0_wp),-1.0_wp)
    end subroutine test_normal_mobility_domain

    subroutine test_requested_derivative_range()
        type(material_t) :: plain,tracked
        real(wp) :: x(ninput,1),velocity(ndim,1),tangent(ndim,ninput,1),dsdp(ndim,nparameters,1)
        real(wp) :: beta
        integer :: status

        beta=huge(1.0_wp)
        call plain%init(darcy_exponential_t(1.0_wp,1.0_wp,beta,1.0_wp),1,status)
        call expect(status,cf_ok)
        call tracked%init(darcy_exponential_t(1.0_wp,1.0_wp,beta,1.0_wp),1,status, &
            sensitivity_parameters=parameter_names)
        call expect(status,cf_ok)
        x(:,1)=[2.0_wp,0.0_wp,0.0_wp,1.0_wp]
        call plain%eval(x,velocity,status=status)
        call expect(status,cf_ok)
        call close(velocity(1,1),-2.0_wp)
        call plain%eval(x,velocity,tangent,status)
        call expect(status,cf_evaluation_failed)
        call plain%commit(status)
        call expect(status,cf_no_trial)
        call tracked%eval(x,velocity,status=status,dsdp=dsdp)
        call expect(status,cf_evaluation_failed)
        if (tracked%last_failure_point()/=1) error stop "Exponential Darcy derivative overflow point"
        call tracked%commit(status)
        call expect(status,cf_no_trial)
        x(1,1)=0.5_wp
        call tracked%eval(x,velocity,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        call close(tangent(1,ninput,1)/beta,-0.5_wp)
        call close(dsdp(1,1,1)/beta,0.5_wp)
    end subroutine test_requested_derivative_range

    subroutine test_invalid_and_recovery()
        type(material_t) :: material
        real(wp) :: x(ninput,2), velocity(ndim,2), invalid(ninput-1,2)
        real(wp) :: extra(ndim,1,2)
        integer :: status

        call material%init(darcy_exponential_t(0.0_wp,2.0_wp,0.1_wp,10.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(darcy_exponential_t(3.0_wp,0.0_wp,0.1_wp,10.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(darcy_exponential_t(3.0_wp,2.0_wp,0.1_wp,-1.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(darcy_exponential_t(3.0_wp,2.0_wp, &
                           ieee_value(0.0_wp,ieee_quiet_nan),10.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(darcy_exponential_t(huge(1.0_wp)/4.0_wp,0.125_wp, &
                           0.0_wp,10.0_wp),2,status)
        call expect(status,cf_invalid_argument)

        call material%init(darcy_exponential_t(3.0_wp,2.0_wp,0.1_wp,10.0_wp),2,status)
        call expect(status,cf_ok)
        x(:,1) = [1.0_wp,-1.0_wp,0.5_wp,12.0_wp]
        x(:,2) = x(:,1)
        x(ninput,2) = -1.0_wp
        call material%eval(x,velocity,status=status)
        call expect(status,cf_evaluation_failed)
        if (material%last_failure_point() /= 2) error stop "Exponential bad point"
        call material%commit(status)
        call expect(status,cf_no_trial)
        x(ninput,2) = 12.0_wp
        x(1,2) = ieee_value(0.0_wp,ieee_quiet_nan)
        call material%eval(x,velocity,status=status)
        call expect(status,cf_evaluation_failed)
        x(1,2) = ieee_value(0.0_wp,ieee_positive_inf)
        call material%eval(x,velocity,status=status)
        call expect(status,cf_evaluation_failed)
        x(1,2) = 1.0_wp
        call material%eval(x,velocity,status=status)
        call expect(status,cf_ok)
        call material%commit(status)
        call expect(status,cf_ok)
        call material%eval(invalid,velocity,status=status)
        call expect(status,cf_invalid_shape)
        call material%eval(x,velocity,status=status,dsdp=extra)
        call expect(status,cf_unsupported)

        x(ninput,2) = 10.0_wp+(log(huge(1.0_wp))+4.0_wp)/0.1_wp
        call material%eval(x,velocity,status=status)
        call expect(status,cf_evaluation_failed)
        x(ninput,2) = 12.0_wp
        x(1,2) = huge(1.0_wp)
        call material%eval(x,velocity,status=status)
        call expect(status,cf_evaluation_failed)
        x(1,2) = 1.0_wp
        call material%eval(x,velocity,status=status)
        call expect(status,cf_ok)

        call material%init(darcy_exponential_t(3.0_wp,2.0_wp,-0.1_wp,10.0_wp),2,status)
        call expect(status,cf_ok)
        x(ninput,2) = 10.0_wp+(-log(tiny(1.0_wp))+4.0_wp)/0.1_wp
        call material%eval(x,velocity,status=status)
        call expect(status,cf_evaluation_failed)
    end subroutine test_invalid_and_recovery

end program test_darcy_exponential

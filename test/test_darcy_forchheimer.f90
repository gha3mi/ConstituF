! SPDX-License-Identifier: Apache-2.0

program test_darcy_forchheimer
    !! Test convention: pressure in Pa, length in m and time in s.
    !! Permeability is m2, viscosity Pa*s and Darcy velocity m/s.
    !!
    !! Independent quadratic-flow, Jacobian, sensitivity and API tests.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
    use constituf, only: wp, material_t, cf_ok, cf_invalid_argument, &
                        cf_invalid_shape, cf_evaluation_failed, cf_no_trial, cf_unsupported
    use constituf_darcy_forchheimer, only: darcy_forchheimer_t
    use constituf_darcy_isotropic, only: darcy_isotropic_t
    implicit none (type, external)

    integer, parameter :: ndim = 3
    character(len=20), parameter :: parameter_names(2) = [character(len=20) :: &
        "quadratic_resistance","linear_resistance"]

    call test_reference_and_schema()
    call test_tangent_and_rotation()
    call test_parameter_derivatives()
    call test_common_resistance_scale(tiny(1.0_wp))
    call test_common_resistance_scale(huge(1.0_wp)/6.0_wp)
    call test_large_zero_field_mobility()
    call test_weak_inertia()
    call test_subnormal_tangent_weight()
    call test_derivative_components()
    call test_large_gradient_norm()
    call test_small_direction_component()
    call test_invalid_and_recovery()

contains

    subroutine expect(actual, wanted)
        integer, intent(in) :: actual, wanted

        if (actual /= wanted) then
            print *, "Forchheimer status", actual, "expected", wanted
            error stop "Forchheimer status mismatch"
        end if
    end subroutine expect

    subroutine close(actual, wanted)
        real(wp), intent(in) :: actual, wanted

        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(wanted)) then
            error stop "Nonfinite Forchheimer comparison"
        end if
        if (abs(actual-wanted) > 512.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(wanted))) then
            print *, "Forchheimer actual", actual, "expected", wanted
            error stop "Forchheimer reference mismatch"
        end if
    end subroutine close

    subroutine test_reference_and_schema()
        type(material_t) :: material, linear
        real(wp) :: gradient(ndim,3), velocity(ndim,3), tangent(ndim,ndim,3)
        real(wp) :: linear_velocity(ndim,3), linear_tangent(ndim,ndim,3)
        real(wp) :: speed, pressure_work
        integer :: status, i, j, ipoint

        call material%init(darcy_forchheimer_t(2.0_wp,3.0_wp),3,status)
        call expect(status,cf_ok)
        if (material%input_size() /= ndim .or. material%output_size() /= ndim) then
            error stop "Forchheimer field sizes"
        end if
        if (material%state_storage_bytes() /= 0 .or. material%state_size() /= 0) then
            error stop "Forchheimer must be stateless"
        end if
        if (.not. material%field_matches("input", "pressure.gradient", &
            "vector", "pressure_gradient", "spatial", [ndim])) then
            error stop "Forchheimer input schema"
        end if
        if (.not. material%field_matches("output", "fluid.darcy_velocity", &
            "vector", "specific_discharge", "spatial", [ndim])) then
            error stop "Forchheimer output schema"
        end if
        if (material%field_matches("input", "pressure.gradient", &
            "vector", "darcy_velocity", "spatial", [ndim])) then
            error stop "Forchheimer incompatible gradient measure accepted"
        end if

        gradient(:,1) = [3.0_wp,4.0_wp,0.0_wp]
        gradient(:,2) = [0.0_wp,0.0_wp,0.0_wp]
        gradient(:,3) = [-2.0_wp,1.0_wp,0.5_wp]
        call material%eval(gradient,velocity,tangent,status)
        call expect(status,cf_ok)
        call close(velocity(1,1),-0.6_wp)
        call close(velocity(2,1),-0.8_wp)
        call close(velocity(3,1),0.0_wp)
        call close(tangent(1,1,1),-0.173_wp)
        call close(tangent(1,2,1),0.036_wp)
        call close(tangent(2,2,1),-0.152_wp)
        call close(tangent(3,3,1),-0.2_wp)
        call close(tangent(1,1,2),-0.5_wp)
        call close(tangent(2,2,2),-0.5_wp)
        call close(tangent(3,3,2),-0.5_wp)
        if (any(abs(velocity(:,2)) > 512.0_wp*epsilon(1.0_wp))) error stop "Forchheimer zero-field velocity"

        do ipoint = 1, 3
            speed = norm2(velocity(:,ipoint))
            pressure_work = -dot_product(gradient(:,ipoint),velocity(:,ipoint))
            if (pressure_work < 0.0_wp) error stop "Forchheimer negative dissipation"
            call close(pressure_work,2.0_wp*speed**2+3.0_wp*speed**3)
            do i = 1, ndim
                call close(gradient(i,ipoint), &
                           -(2.0_wp+3.0_wp*speed)*velocity(i,ipoint))
                do j = 1, ndim
                    call close(tangent(i,j,ipoint),tangent(j,i,ipoint))
                end do
            end do
        end do
        call material%commit(status)
        call expect(status,cf_ok)
        call material%revert()
        call material%commit(status)
        call expect(status,cf_no_trial)

        ! Setting quadratic resistance to zero recovers Darcy k/mu = 1/a.
        call material%init(darcy_forchheimer_t(2.0_wp,0.0_wp),3,status)
        call expect(status,cf_ok)
        call linear%init(darcy_isotropic_t(3.0_wp,6.0_wp),3,status)
        call expect(status,cf_ok)
        call material%eval(gradient,velocity,tangent,status)
        call expect(status,cf_ok)
        call linear%eval(gradient,linear_velocity,linear_tangent,status)
        call expect(status,cf_ok)
        do ipoint = 1, 3
            do j = 1, ndim
                call close(velocity(j,ipoint),linear_velocity(j,ipoint))
                do i = 1, ndim
                    call close(tangent(i,j,ipoint),linear_tangent(i,j,ipoint))
                end do
            end do
        end do
    end subroutine test_reference_and_schema

    subroutine test_tangent_and_rotation()
        type(material_t) :: batch, single
        real(wp) :: gradient(ndim,2), velocity(ndim,2), tangent(ndim,ndim,2)
        real(wp) :: alone(ndim,1), changed(ndim,2), plus(ndim,2), minus(ndim,2)
        real(wp) :: h, derivative, radius, bound
        integer :: status, ipoint, i, j

        call batch%init(darcy_forchheimer_t(2.5_wp,0.7_wp),2,status)
        call expect(status,cf_ok)
        call single%init(darcy_forchheimer_t(2.5_wp,0.7_wp),1,status)
        call expect(status,cf_ok)
        gradient(:,1) = [3.2_wp,-1.7_wp,0.4_wp]
        gradient(:,2) = [-gradient(2,1),gradient(1,1),gradient(3,1)]
        call batch%eval(gradient,velocity,tangent,status)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(tangent))) error stop "Nonfinite Forchheimer analytical tangent"
        call single%eval(gradient(:,1:1),alone,status=status)
        call expect(status,cf_ok)
        do i = 1, ndim
            call close(velocity(i,1),alone(i,1))
        end do
        call close(velocity(1,2),-velocity(2,1))
        call close(velocity(2,2),velocity(1,1))
        call close(velocity(3,2),velocity(3,1))

        do ipoint = 1, 2
            radius=norm2(gradient(:,ipoint))
            do j = 1, ndim
                h = epsilon(1.0_wp)**(1.0_wp/3.0_wp)*radius
                changed = gradient
                changed(j,ipoint) = gradient(j,ipoint)+h
                call batch%eval(changed,plus,status=status)
                call expect(status,cf_ok)
                changed(j,ipoint) = gradient(j,ipoint)-h
                call batch%eval(changed,minus,status=status)
                call expect(status,cf_ok)
                if (.not. all(ieee_is_finite(plus)) .or. .not. all(ieee_is_finite(minus))) then
                    error stop "Nonfinite Forchheimer tangent secant response"
                end if
                do i = 1, ndim
                    derivative = (plus(i,ipoint)-minus(i,ipoint))/(2.0_wp*h)
                    ! For v=-m(r)g, positive a,b imply an entrywise third
                    ! derivative bound 18/(a*r**2). Use the smallest radius
                    ! in the perturbation interval for the central bias.
                    bound=32.0_wp*epsilon(1.0_wp)*(abs(plus(i,ipoint))+abs(minus(i,ipoint)))/(2.0_wp*h)+ &
                        3.0_wp*(h/(radius-h))**2/2.5_wp
                    if (.not. all(ieee_is_finite([derivative,bound]))) then
                        error stop "Nonfinite Forchheimer tangent comparison"
                    end if
                    if (abs(derivative-tangent(i,j,ipoint))>bound) then
                        error stop "Forchheimer exact tangent finite difference"
                    end if
                    call close(tangent(i,j,ipoint),tangent(j,i,ipoint))
                end do
            end do
        end do
    end subroutine test_tangent_and_rotation

    subroutine test_parameter_derivatives()
        type(material_t) :: tracked, upper, lower
        real(wp), parameter :: direction(ndim)=[1.0_wp/3.0_wp,2.0_wp/3.0_wp,-2.0_wp/3.0_wp]
        real(wp), parameter :: exact_tangent(ndim,ndim)=reshape( &
            [-5.0_wp/42.0_wp,1.0_wp/84.0_wp,-1.0_wp/84.0_wp, &
             1.0_wp/84.0_wp,-17.0_wp/168.0_wp,-1.0_wp/42.0_wp, &
             -1.0_wp/84.0_wp,-1.0_wp/42.0_wp,-17.0_wp/168.0_wp],[ndim,ndim])
        real(wp) :: gradient(ndim,1), velocity(ndim,1), dsdp(ndim,2,1)
        real(wp) :: tangent(ndim,ndim,1)
        real(wp) :: plus(ndim,1), minus(ndim,1), wrong(ndim,1,1)
        real(wp) :: derivative,h,radius,minimum_root,maximum_speed,third_bound,bound
        integer :: status, iparam, i, j

        call tracked%init(darcy_forchheimer_t(2.0_wp,3.0_wp),1,status,sensitivity_parameters=parameter_names)
        call expect(status,cf_ok)
        if (tracked%sensitivity_size() /= 2) error stop "Forchheimer sensitivity count"
        if (tracked%state_storage_bytes() /= 0) error stop "Forchheimer sensitivity state"
        ! The independent speed s=2 solves 3*s**2+2*s=16; its two
        ! resistance sensitivities are distinct in every nonzero component.
        radius=16.0_wp
        gradient(:,1)=radius*direction
        call tracked%eval(gradient,velocity,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        do i = 1, ndim
            call close(velocity(i,1),-2.0_wp*direction(i))
            call close(dsdp(i,1,1),(2.0_wp/7.0_wp)*direction(i))
            call close(dsdp(i,2,1),direction(i)/7.0_wp)
            do j=1,ndim
                call close(tangent(i,j,1),exact_tangent(i,j))
            end do
        end do
        call close(-dot_product(gradient(:,1),velocity(:,1)),32.0_wp)

        do iparam = 1, 2
            select case (iparam)
            case (1)
                h=3.0_wp*epsilon(1.0_wp)**(1.0_wp/3.0_wp)
                minimum_root=sqrt(4.0_wp+4.0_wp*(3.0_wp-h)*radius)
                maximum_speed=2.0_wp*radius/(2.0_wp+minimum_root)
                third_bound=6.0_wp*maximum_speed**4/minimum_root**3+ &
                    12.0_wp*radius*maximum_speed**3/minimum_root**4+ &
                    12.0_wp*radius**2*maximum_speed**2/minimum_root**5
                call upper%init(darcy_forchheimer_t(2.0_wp,3.0_wp+h),1,status)
                call expect(status,cf_ok)
                call lower%init(darcy_forchheimer_t(2.0_wp,3.0_wp-h),1,status)
            case (2)
                h=2.0_wp*epsilon(1.0_wp)**(1.0_wp/3.0_wp)
                minimum_root=sqrt((2.0_wp-h)**2+12.0_wp*radius)
                third_bound=6.0_wp*radius*(2.0_wp+h)/minimum_root**5
                call upper%init(darcy_forchheimer_t(2.0_wp+h,3.0_wp),1,status)
                call expect(status,cf_ok)
                call lower%init(darcy_forchheimer_t(2.0_wp-h,3.0_wp),1,status)
            case default
                error stop "Invalid Forchheimer parameter index"
            end select
            call expect(status,cf_ok)
            call upper%eval(gradient,plus,status=status)
            call expect(status,cf_ok)
            call lower%eval(gradient,minus,status=status)
            call expect(status,cf_ok)
            if (.not. all(ieee_is_finite(plus)) .or. .not. all(ieee_is_finite(minus))) then
                error stop "Nonfinite Forchheimer parameter secant response"
            end if
            do i = 1, ndim
                derivative = (plus(i,1)-minus(i,1))/(2.0_wp*h)
                ! Bound the central bias by the analytical third derivative
                ! over the complete parameter perturbation interval.
                bound=32.0_wp*epsilon(1.0_wp)*(abs(plus(i,1))+abs(minus(i,1)))/(2.0_wp*h)+ &
                    abs(direction(i))*h*h*third_bound/6.0_wp
                if (.not. all(ieee_is_finite([derivative,bound]))) then
                    error stop "Nonfinite Forchheimer parameter comparison"
                end if
                if (abs(derivative-dsdp(i,iparam,1))>bound) then
                    error stop "Forchheimer parameter finite difference"
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
        gradient=0.0_wp
        call tracked%eval(gradient,velocity,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(velocity)) .or. .not. all(ieee_is_finite(dsdp))) then
            error stop "Nonfinite zero-field Forchheimer response"
        end if
        call close(maxval(abs(velocity)),0.0_wp)
        call close(maxval(abs(dsdp)),0.0_wp)
        do j=1,ndim
            do i=1,ndim
                if (i==j) then
                    call close(tangent(i,j,1),-0.5_wp)
                else
                    call close(tangent(i,j,1),0.0_wp)
                end if
            end do
        end do
    end subroutine test_parameter_derivatives

    subroutine test_common_resistance_scale(resistance_scale)
        !! The speed remains one under common scaling of a,b and gradient.
        !! For the large scale, d=8*scale overflows although every requested
        !! response and derivative is finite.
        real(wp), intent(in) :: resistance_scale
        type(material_t) :: material
        real(wp), parameter :: exact_tangent(ndim,ndim)=reshape( &
            [-173.0_wp/1000.0_wp,9.0_wp/250.0_wp,0.0_wp, &
             9.0_wp/250.0_wp,-19.0_wp/125.0_wp,0.0_wp, &
             0.0_wp,0.0_wp,-0.2_wp],[ndim,ndim])
        real(wp), parameter :: parameter_reference(ndim)=[3.0_wp/40.0_wp,0.1_wp,0.0_wp]
        real(wp) :: gradient(ndim,1),velocity(ndim,1),tangent(ndim,ndim,1),dsdp(ndim,2,1)
        integer :: status,i,j,k

        call material%init(darcy_forchheimer_t(2.0_wp*resistance_scale,3.0_wp*resistance_scale),1,status, &
            sensitivity_parameters=parameter_names)
        call expect(status,cf_ok)
        gradient(:,1)=resistance_scale*[3.0_wp,4.0_wp,0.0_wp]
        call material%eval(gradient,velocity,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        call close(velocity(1,1),-0.6_wp)
        call close(velocity(2,1),-0.8_wp)
        call close(velocity(3,1),0.0_wp)
        do j=1,ndim
            do i=1,ndim
                call close(tangent(i,j,1)*resistance_scale,exact_tangent(i,j))
            end do
        end do
        do k=1,2
            do i=1,ndim
                call close(dsdp(i,k,1)*resistance_scale,parameter_reference(i))
            end do
        end do
    end subroutine test_common_resistance_scale

    subroutine test_large_zero_field_mobility()
        type(material_t) :: material
        real(wp) :: gradient(ndim,1),velocity(ndim,1),tangent(ndim,ndim,1),dsdp(ndim,2,1)
        real(wp) :: a
        integer :: status,i,j

        ! The admitted mobility 1/a is finite, while forming 2/a first is not.
        a=tiny(1.0_wp)/2.0_wp
        call material%init(darcy_forchheimer_t(a,0.0_wp),1,status,sensitivity_parameters=parameter_names)
        call expect(status,cf_ok)
        gradient=0.0_wp
        call material%eval(gradient,velocity,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(velocity)) .or. .not. all(ieee_is_finite(dsdp))) then
            error stop "Nonfinite large-mobility zero-field response"
        end if
        call close(maxval(abs(velocity)),0.0_wp)
        call close(maxval(abs(dsdp)),0.0_wp)
        do j=1,ndim
            do i=1,ndim
                if (i==j) then
                    call close(tangent(i,j,1)*a,-1.0_wp)
                else
                    call close(tangent(i,j,1)*a,0.0_wp)
                end if
            end do
        end do
    end subroutine test_large_zero_field_mobility

    subroutine test_weak_inertia()
        !! These cases retain a small but representable off-diagonal tangent:
        !! first below subtractive resolution, then below the range of the
        !! dimensionless inertial ratio before its dimensional factors restore it.
        type(material_t) :: material
        real(wp), parameter :: direction(ndim)=[0.6_wp,0.8_wp,0.0_wp]
        real(wp) :: gradient(ndim,1),velocity(ndim,1),tangent(ndim,ndim,1),dsdp(ndim,2,1)
        real(wp) :: a,b,radius,velocity_scale,linear_derivative_scale,quadratic_derivative_scale,coupling_scale
        integer :: status,icase,i

        do icase=1,2
            if (icase==1) then
                a=1.0_wp
                b=epsilon(1.0_wp)**2
                radius=1.0_wp
            else
                a=tiny(1.0_wp)**(1.0_wp/3.0_wp)
                b=tiny(1.0_wp)
                radius=tiny(1.0_wp)
            end if
            velocity_scale=radius/a
            linear_derivative_scale=velocity_scale/a
            quadratic_derivative_scale=linear_derivative_scale*velocity_scale
            coupling_scale=((b/a)/a)*(radius/a)
            call material%init(darcy_forchheimer_t(a,b),1,status,sensitivity_parameters=parameter_names)
            call expect(status,cf_ok)
            gradient(:,1)=radius*direction
            call material%eval(gradient,velocity,tangent,status,dsdp=dsdp)
            call expect(status,cf_ok)
            do i=1,ndim
                call close(velocity(i,1)/velocity_scale,-direction(i))
                call close(dsdp(i,1,1)/quadratic_derivative_scale,direction(i))
                call close(dsdp(i,2,1)/linear_derivative_scale,direction(i))
                call close(tangent(i,i,1)*a,-1.0_wp)
            end do
            call close(tangent(1,2,1)/coupling_scale,0.48_wp)
            call close(tangent(2,1,1)/coupling_scale,0.48_wp)
        end do
    end subroutine test_weak_inertia

    subroutine test_subnormal_tangent_weight()
        !! The unrounded inertial weight is sqrt(2)*d. Multiplication by
        !! the signed direction product +/-1/2 must precede restoration.
        type(material_t) :: material
        integer, parameter :: ncases=2
        real(wp), parameter :: signs(ncases)=[1.0_wp,-1.0_wp]
        real(wp) :: gradient(ndim,1),velocity(ndim,1),tangent(ndim,ndim,1),dsdp(ndim,2,1)
        real(wp) :: smallest
        integer :: status,icase,i

        smallest=nearest(0.0_wp,1.0_wp)
        call material%init(darcy_forchheimer_t(1.0_wp,smallest),1,status,sensitivity_parameters=parameter_names)
        call expect(status,cf_ok)
        do icase=1,ncases
            gradient(:,1)=[signs(icase),1.0_wp,0.0_wp]
            call material%eval(gradient,velocity,tangent,status,dsdp=dsdp)
            call expect(status,cf_ok)
            call close(tangent(1,2,1)/smallest,signs(icase))
            call close(tangent(2,1,1)/smallest,signs(icase))
            do i=1,ndim
                call close(velocity(i,1),-gradient(i,1))
                call close(tangent(i,i,1),-1.0_wp)
                call close(dsdp(i,1,1),sqrt(2.0_wp)*gradient(i,1))
                call close(dsdp(i,2,1),gradient(i,1))
            end do
        end do
    end subroutine test_subnormal_tangent_weight

    subroutine test_derivative_components()
        !! An individual velocity component can round to zero while its
        !! quadratic-resistance derivative remains normal and nonzero.
        type(material_t) :: tracked,linear_only
        real(wp) :: gradient(ndim,1),velocity(ndim,1),tangent(ndim,ndim,1),dsdp(ndim,2,1)
        real(wp) :: linear_dsdp(ndim,1,1)
        real(wp) :: small,field_scale,small_derivative
        integer :: status

        small=nearest(0.0_wp,1.0_wp)
        field_scale=sqrt(huge(1.0_wp))
        small_derivative=small*(field_scale/8.0_wp)
        call tracked%init(darcy_forchheimer_t(2.0_wp,0.0_wp),1,status,sensitivity_parameters=parameter_names)
        call expect(status,cf_ok)
        gradient(:,1)=[small,field_scale,0.0_wp]
        call tracked%eval(gradient,velocity,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        call close(velocity(1,1),0.0_wp)
        call close(velocity(2,1)/field_scale,-0.5_wp)
        call close(dsdp(1,1,1)/small_derivative,1.0_wp)
        call close(dsdp(2,1,1)/huge(1.0_wp),0.125_wp)
        call close(dsdp(2,2,1)/field_scale,0.25_wp)

        ! Increasing the large component makes only the quadratic derivative
        ! unrepresentable. A selected linear derivative is still valid.
        call linear_only%init(darcy_forchheimer_t(2.0_wp,0.0_wp),1,status, &
            sensitivity_parameters=[character(len=17) :: "linear_resistance"])
        call expect(status,cf_ok)
        gradient(2,1)=4.0_wp*field_scale
        call linear_only%eval(gradient,velocity,tangent,status,dsdp=linear_dsdp)
        call expect(status,cf_ok)
        call close(velocity(2,1)/field_scale,-2.0_wp)
        call close(linear_dsdp(2,1,1)/field_scale,1.0_wp)
        call tracked%eval(gradient,velocity,status=status,dsdp=dsdp)
        call expect(status,cf_evaluation_failed)
        if (tracked%last_failure_point()/=1) error stop "Forchheimer derivative overflow point"
        call tracked%commit(status)
        call expect(status,cf_no_trial)
        gradient(2,1)=field_scale
        call tracked%eval(gradient,velocity,status=status,dsdp=dsdp)
        call expect(status,cf_ok)
    end subroutine test_derivative_components

    subroutine test_large_gradient_norm()
        type(material_t) :: material
        real(wp) :: gradient(ndim,1),velocity(ndim,1),tangent(ndim,ndim,1),dsdp(ndim,1,1)
        real(wp) :: large
        integer :: status,i,j

        ! In the Darcy limit the vector norm is unnecessary: its overflow
        ! must not reject finite component responses or selected derivatives.
        large=huge(1.0_wp)
        call material%init(darcy_forchheimer_t(2.0_wp,0.0_wp),1,status, &
            sensitivity_parameters=[character(len=17) :: "linear_resistance"])
        call expect(status,cf_ok)
        gradient=0.75_wp*large
        call material%eval(gradient,velocity,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        do i=1,ndim
            call close(velocity(i,1)/large,-0.375_wp)
            call close(dsdp(i,1,1)/large,0.1875_wp)
            do j=1,ndim
                if (i==j) then
                    call close(tangent(i,j,1),-0.5_wp)
                else
                    call close(tangent(i,j,1),0.0_wp)
                end if
            end do
        end do
    end subroutine test_large_gradient_norm

    subroutine test_small_direction_component()
        type(material_t) :: material
        real(wp) :: gradient(ndim,1),velocity(ndim,1),tangent(ndim,ndim,1)
        real(wp) :: small,resistance,reference
        integer :: status

        ! Dividing the smallest gradient component by its vector norm first
        ! quantizes that component. The final off-diagonal tangent is normal.
        small=nearest(0.0_wp,1.0_wp)
        resistance=tiny(1.0_wp)
        reference=(small/sqrt(resistance))/(4.0_wp*sqrt(sqrt(2.0_wp)))
        call material%init(darcy_forchheimer_t(resistance,resistance),1,status)
        call expect(status,cf_ok)
        gradient(:,1)=[small,1.0_wp,1.0_wp]
        call material%eval(gradient,velocity,tangent,status)
        call expect(status,cf_ok)
        call close(tangent(1,2,1)/reference,1.0_wp)
        call close(tangent(1,3,1)/reference,1.0_wp)
    end subroutine test_small_direction_component

    subroutine test_invalid_and_recovery()
        type(material_t) :: material
        real(wp) :: gradient(ndim,2), velocity(ndim,2), wrong(ndim-1,2)
        real(wp) :: extra(ndim,1,2)
        integer :: status

        call material%init(darcy_forchheimer_t(0.0_wp,1.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(darcy_forchheimer_t(-1.0_wp,1.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(darcy_forchheimer_t(1.0_wp,-1.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(darcy_forchheimer_t(1.0_wp, &
                           ieee_value(0.0_wp,ieee_quiet_nan)),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(darcy_forchheimer_t(2.0_wp,3.0_wp),2,status)
        call expect(status,cf_ok)
        gradient = 0.0_wp
        gradient(1,:) = 1.0_wp
        gradient(2,2) = ieee_value(0.0_wp,ieee_positive_inf)
        call material%eval(gradient,velocity,status=status)
        call expect(status,cf_evaluation_failed)
        if (material%last_failure_point() /= 2) error stop "Forchheimer failing point"
        call material%commit(status)
        call expect(status,cf_no_trial)
        gradient(2,2) = ieee_value(0.0_wp,ieee_quiet_nan)
        call material%eval(gradient,velocity,status=status)
        call expect(status,cf_evaluation_failed)
        gradient(2,2) = 0.0_wp
        call material%eval(gradient,velocity,status=status)
        call expect(status,cf_ok)
        if (.not. all(velocity(1,:) < 0.0_wp)) error stop "Forchheimer recovery"
        call material%eval(wrong,velocity,status=status)
        call expect(status,cf_invalid_shape)
        call material%eval(gradient,velocity,status=status,dsdp=extra)
        call expect(status,cf_unsupported)

        ! A representable zero-field mobility can overflow for a huge
        ! gradient; a failed point may not be committed.
        call material%init(darcy_forchheimer_t(1.0_wp/sqrt(huge(1.0_wp)),0.0_wp),2,status)
        call expect(status,cf_ok)
        gradient(:,1) = [0.0_wp,0.0_wp,0.0_wp]
        gradient(:,2) = [4.0_wp*sqrt(huge(1.0_wp)),0.0_wp,0.0_wp]
        call material%eval(gradient,velocity,status=status)
        call expect(status,cf_evaluation_failed)
        call material%commit(status)
        call expect(status,cf_no_trial)
    end subroutine test_invalid_and_recovery

end program test_darcy_forchheimer

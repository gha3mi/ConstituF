! SPDX-License-Identifier: Apache-2.0

program test_fourier_exponential
    !! Test convention: absolute temperature in K, length in m and time in s.
    !! Conductivity is W/(m*K), temperature gradient K/m and heat flux W/m2.
    !!
    !! Independent Fourier limits, coupled tangents and parameter sensitivity tests.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
    use constituf, only: wp, material_t, cf_ok, cf_invalid_argument, &
                        cf_evaluation_failed, cf_invalid_shape, cf_no_trial
    use constituf_fourier_exponential, only: fourier_exponential_t
    use constituf_fourier_isotropic, only: fourier_isotropic_t
    implicit none (type, external)

    integer, parameter :: ndim = 3, ninput = ndim+1, nparameters = 3

    call test_reference_and_metadata()
    call test_tangent_and_limit()
    call test_parameter_sensitivities()
    call test_reference_range(tiny(1.0_wp))
    call test_reference_range(huge(1.0_wp)/4.0_wp)
    call test_balanced_derivatives(huge(1.0_wp)/16.0_wp,32.0_wp,1.0_wp/32.0_wp)
    call test_balanced_derivatives(8.0_wp*tiny(1.0_wp),1.0_wp/tiny(1.0_wp),tiny(1.0_wp))
    call test_reference_ratio()
    call test_subnormal_conductivity()
    call test_invalid_trials()

contains

    subroutine expect(actual, expected)
        integer, intent(in) :: actual, expected

        if (actual /= expected) then
            print *, "Status", actual, "expected", expected
            error stop "Exponential Fourier status mismatch"
        end if
    end subroutine expect

    subroutine close(actual, expected)
        real(wp), intent(in) :: actual, expected

        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) then
            error stop "Nonfinite exponential Fourier comparison"
        end if
        if (abs(actual-expected) > 512.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(expected))) then
            print *, "Actual", actual, "expected", expected
            error stop "Exponential Fourier reference mismatch"
        end if
    end subroutine close

    subroutine test_reference_and_metadata()
        type(material_t) :: material
        real(wp), parameter :: coefficient=log(2.0_wp)/10.0_wp
        real(wp), parameter :: exact_conductivity(2)=[10.0_wp,5.0_wp]
        real(wp) :: x(ninput,2), flux(ndim,2), tangent(ndim,ninput,2)
        real(wp) :: dsdp(ndim,nparameters,2)
        real(wp) :: conductivity
        integer :: status, ipoint, i, j
        integer :: empty_shape(0)

        call material%init(fourier_exponential_t(5.0_wp,coefficient,300.0_wp),2,status, &
            sensitivity_parameters=[character(len=23) :: "reference_temperature", &
                "temperature_coefficient","reference_conductivity"])
        call expect(status,cf_ok)
        if (material%input_size() /= ninput .or. material%output_size() /= ndim) then
            error stop "Temperature-conduction field sizes"
        end if
        if (material%state_storage_bytes() /= 0) error stop "Unexpected thermal history"
        if (material%input_offset("temperature.gradient") /= 1) then
            error stop "Wrong gradient offset"
        end if
        if (material%input_offset("temperature.absolute") /= ninput) then
            error stop "Wrong temperature offset"
        end if
        if (.not. material%field_matches("input","temperature.gradient", &
            "vector","temperature_gradient","spatial",[ndim])) then
            error stop "Thermal gradient metadata mismatch"
        end if
        if (.not. material%field_matches("input","temperature.absolute", &
            "scalar","absolute_temperature","material",empty_shape)) then
            error stop "Absolute temperature metadata mismatch"
        end if
        if (.not. material%field_matches("output","heat.flux", &
            "vector","heat_flux","spatial",[ndim])) then
            error stop "Heat flux metadata mismatch"
        end if

        x(:,1) = [2.0_wp,-1.0_wp,3.0_wp,310.0_wp]
        x(:,2) = [-0.5_wp,1.5_wp,0.0_wp,300.0_wp]
        call material%eval(x,flux,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        do ipoint = 1, 2
            ! exp(log(2)) doubles the reference conductivity; the second
            ! point is exactly at T0 and must retain k0.
            conductivity = exact_conductivity(ipoint)
            do i = 1, ndim
                call close(flux(i,ipoint),-conductivity*x(i,ipoint))
                call close(tangent(i,ninput,ipoint), &
                           -coefficient*conductivity*x(i,ipoint))
                call close(dsdp(i,1,ipoint),coefficient*conductivity*x(i,ipoint))
                call close(dsdp(i,2,ipoint), &
                    -conductivity*(x(ninput,ipoint)-300.0_wp)*x(i,ipoint))
                call close(dsdp(i,3,ipoint),-(conductivity/5.0_wp)*x(i,ipoint))
                do j = 1, ndim
                    if (i == j) then
                        call close(tangent(i,j,ipoint),-conductivity)
                    else
                        call close(tangent(i,j,ipoint),0.0_wp)
                    end if
                end do
            end do
            if (dot_product(x(1:ndim,ipoint),flux(:,ipoint)) > 0.0_wp) then
                error stop "Fourier thermal dissipation violated"
            end if
        end do
        call close(-dot_product(flux(:,1),x(1:ndim,1)),140.0_wp)
        call material%commit(status)
        call expect(status,cf_ok)
        call material%revert()
        call material%commit(status)
        call expect(status,cf_no_trial)
    end subroutine test_reference_and_metadata

    subroutine test_tangent_and_limit()
        type(material_t) :: material, isotropic
        real(wp) :: x(ninput,2), flux(ndim,2), tangent(ndim,ninput,2)
        real(wp) :: changed(ninput,2), plus(ndim,2), minus(ndim,2)
        real(wp) :: constant_flux(ndim,2)
        real(wp) :: h, derivative, error_bound
        integer :: status, i, j

        call material%init(fourier_exponential_t(7.0_wp,-0.008_wp,295.0_wp),2,status)
        call expect(status,cf_ok)
        x(:,1) = [1.0_wp,-2.0_wp,0.4_wp,305.0_wp]
        x(:,2) = [-x(2,1),x(1,1),x(3,1),305.0_wp]
        call material%eval(x,flux,tangent,status)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(tangent))) error stop "Nonfinite exponential Fourier analytical tangent"
        call close(flux(1,2),-flux(2,1))
        call close(flux(2,2),flux(1,1))
        call close(flux(3,2),flux(3,1))

        do j = 1, ninput
            if (j == ninput) then
                h = epsilon(1.0_wp)**(1.0_wp/3.0_wp)/0.008_wp
            else
                ! Flux is affine in each gradient component.
                h = 0.25_wp*maxval(abs(x(1:ndim,1)))
            end if
            changed = x
            changed(j,1) = changed(j,1)+h
            call material%eval(changed,plus,status=status)
            call expect(status,cf_ok)
            changed(j,1) = changed(j,1)-2.0_wp*h
            call material%eval(changed,minus,status=status)
            call expect(status,cf_ok)
            if (.not. all(ieee_is_finite(plus)) .or. .not. all(ieee_is_finite(minus))) then
                error stop "Nonfinite exponential Fourier tangent response"
            end if
            do i = 1, ndim
                derivative = (plus(i,1)-minus(i,1))/(2.0_wp*h)
                error_bound = 32.0_wp*epsilon(1.0_wp)*maxval(abs(plus)+abs(minus))/(2.0_wp*h)
                if (j == ninput) then
                    ! Central exponential differences have relative truncation
                    ! (beta*h)**2/6, with h set by the exponent scale.
                    error_bound = error_bound+(0.008_wp*h)**2*abs(tangent(i,j,1))/3.0_wp
                end if
                if (.not. all(ieee_is_finite([derivative,error_bound]))) then
                    error stop "Nonfinite exponential Fourier tangent comparison"
                end if
                if (abs(derivative-tangent(i,j,1)) > error_bound) then
                    error stop "Nonlinear Fourier input tangent mismatch"
                end if
            end do
        end do

        call material%init(fourier_exponential_t(7.0_wp,0.0_wp,295.0_wp),2,status)
        call expect(status,cf_ok)
        call isotropic%init(fourier_isotropic_t(7.0_wp),2,status)
        call expect(status,cf_ok)
        call material%eval(x,flux,tangent,status)
        call expect(status,cf_ok)
        call isotropic%eval(x(1:ndim,:),constant_flux,status=status)
        call expect(status,cf_ok)
        do j = 1, 2
            do i = 1, ndim
                call close(flux(i,j),constant_flux(i,j))
                call close(tangent(i,ninput,j),0.0_wp)
            end do
        end do
    end subroutine test_tangent_and_limit

    subroutine test_parameter_sensitivities()
        type(material_t) :: tracked, plus_material, minus_material
        real(wp) :: x(ninput,1), flux(ndim,1), dsdp(ndim,3,1)
        real(wp) :: tangent(ndim,ninput,1), plus(ndim,1), minus(ndim,1)
        real(wp) :: h, derivative, error_bound, exponent_step
        integer :: status, iparam, i

        x(:,1) = [1.0_wp,-0.5_wp,2.5_wp,320.0_wp]
        call tracked%init(fourier_exponential_t(8.0_wp,0.007_wp,300.0_wp),1,status, &
            sensitivity_parameters=[character(len=23) :: "reference_temperature", &
                                    "temperature_coefficient","reference_conductivity"])
        call expect(status,cf_ok)
        if (tracked%sensitivity_size() /= 3) error stop "Missing parameter derivatives"
        if (tracked%state_storage_bytes() /= 0) error stop "Thermal derivative state allocated"
        call tracked%eval(x,flux,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(dsdp))) error stop "Nonfinite exponential Fourier analytical sensitivity"

        do iparam = 1, 3
            select case (iparam)
            case (1)
                h = epsilon(1.0_wp)**(1.0_wp/3.0_wp)/0.007_wp
                exponent_step = 0.007_wp*h
                call plus_material%init( &
                    fourier_exponential_t(8.0_wp,0.007_wp,300.0_wp+h),1,status)
                call expect(status,cf_ok)
                call minus_material%init( &
                    fourier_exponential_t(8.0_wp,0.007_wp,300.0_wp-h),1,status)
            case (2)
                h = epsilon(1.0_wp)**(1.0_wp/3.0_wp)/20.0_wp
                exponent_step = 20.0_wp*h
                call plus_material%init( &
                    fourier_exponential_t(8.0_wp,0.007_wp+h,300.0_wp),1,status)
                call expect(status,cf_ok)
                call minus_material%init( &
                    fourier_exponential_t(8.0_wp,0.007_wp-h,300.0_wp),1,status)
            case (3)
                h = 0.25_wp*8.0_wp
                exponent_step = 0.0_wp
                call plus_material%init( &
                    fourier_exponential_t(8.0_wp+h,0.007_wp,300.0_wp),1,status)
                call expect(status,cf_ok)
                call minus_material%init( &
                    fourier_exponential_t(8.0_wp-h,0.007_wp,300.0_wp),1,status)
            case default
                error stop "Unexpected sensitivity index"
            end select
            call expect(status,cf_ok)
            call plus_material%eval(x,plus,status=status)
            call expect(status,cf_ok)
            call minus_material%eval(x,minus,status=status)
            call expect(status,cf_ok)
            if (.not. all(ieee_is_finite(plus)) .or. .not. all(ieee_is_finite(minus))) then
                error stop "Nonfinite exponential Fourier parameter response"
            end if
            do i = 1, ndim
                derivative = (plus(i,1)-minus(i,1))/(2.0_wp*h)
                error_bound = 32.0_wp*epsilon(1.0_wp)*maxval(abs(plus)+abs(minus))/(2.0_wp*h)+ &
                    exponent_step**2*abs(dsdp(i,iparam,1))/3.0_wp
                if (.not. all(ieee_is_finite([derivative,error_bound]))) then
                    error stop "Nonfinite exponential Fourier parameter comparison"
                end if
                if (abs(derivative-dsdp(i,iparam,1)) > error_bound) then
                    error stop "Thermal parameter sensitivity mismatch"
                end if
            end do
        end do
        call tracked%commit(status)
        call expect(status,cf_ok)
        call tracked%eval(x,flux,status=status)
        call expect(status,cf_invalid_argument)
        call tracked%eval(x,flux,status=status,dsdp=dsdp)
        call expect(status,cf_ok)
    end subroutine test_parameter_sensitivities

    subroutine test_reference_range(reference)
        !! beta=0 and T=T0 preserve finite positive k0 at the range edges.
        real(wp), intent(in) :: reference
        type(material_t) :: material
        real(wp) :: x(ninput,1),flux(ndim,1),tangent(ndim,ninput,1)
        real(wp) :: dsdp(ndim,nparameters,1)
        integer :: status,i,j

        call material%init(fourier_exponential_t(reference,0.0_wp,1.0_wp),1,status, &
            sensitivity_parameters=[character(len=23) :: "reference_temperature", &
                "temperature_coefficient","reference_conductivity"])
        call expect(status,cf_ok)
        x(:,1)=[1.0_wp,-0.5_wp,0.0_wp,1.0_wp]
        call material%eval(x,flux,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        do i=1,ndim
            call close(flux(i,1)/reference,-x(i,1))
            call close(tangent(i,ninput,1),0.0_wp)
            call close(dsdp(i,1,1),0.0_wp)
            call close(dsdp(i,2,1),0.0_wp)
            call close(dsdp(i,3,1),-x(i,1))
            do j=1,ndim
                if (i==j) then
                    call close(tangent(i,j,1)/reference,-1.0_wp)
                else
                    call close(tangent(i,j,1)/reference,0.0_wp)
                end if
            end do
        end do
    end subroutine test_reference_range

    subroutine test_balanced_derivatives(reference,coefficient,gradient_scale)
        !! Product factors must be balanced independently for each derivative:
        !! k*beta can overflow, while beta*q can erase an underflowed flux.
        real(wp), intent(in) :: reference,coefficient,gradient_scale
        type(material_t) :: material
        real(wp) :: x(ninput,1),flux(ndim,1),tangent(ndim,ninput,1)
        real(wp) :: dsdp(ndim,nparameters,1)
        real(wp), parameter :: direction(ndim)=[1.0_wp,-0.5_wp,0.0_wp]
        integer :: status,i

        call material%init(fourier_exponential_t(reference,coefficient,1.0_wp),1,status, &
            sensitivity_parameters=[character(len=23) :: "reference_temperature", &
                "temperature_coefficient","reference_conductivity"])
        call expect(status,cf_ok)
        x(1:ndim,1)=gradient_scale*direction
        x(ninput,1)=1.0_wp
        call material%eval(x,flux,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        do i=1,ndim
            call close(tangent(i,ninput,1)/reference,-direction(i))
            call close(dsdp(i,1,1)/reference,direction(i))
            call close(dsdp(i,2,1),0.0_wp)
            call close(dsdp(i,3,1)/gradient_scale,-direction(i))
            call close(tangent(i,i,1)/reference,-1.0_wp)
        end do
    end subroutine test_balanced_derivatives

    subroutine test_reference_ratio()
        !! k/k0 overflows although multiplying that ratio by the small
        !! gradient gives a finite reference-conductivity derivative.
        type(material_t) :: material
        real(wp) :: x(ninput,1),flux(ndim,1),tangent(ndim,ninput,1),dsdp(ndim,1,1)
        real(wp) :: log_small,error_bound
        integer :: status

        log_small=log(tiny(1.0_wp))
        call material%init(fourier_exponential_t(tiny(1.0_wp),1.0_wp,0.0_wp),1,status, &
            sensitivity_parameters=[character(len=23) :: "reference_conductivity"])
        call expect(status,cf_ok)
        x(:,1)=[tiny(1.0_wp),0.0_wp,0.0_wp,log(16.0_wp)-log_small]
        call material%eval(x,flux,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        ! Forming T from log(tiny) incurs O(eps*abs(log(tiny))) relative error in k.
        error_bound=256.0_wp*epsilon(1.0_wp)*abs(log_small)
        if (.not. all(ieee_is_finite([flux(1,1),tangent(1,1,1),dsdp(1,1,1),error_bound]))) then
            error stop "Nonfinite balanced Fourier ratio"
        end if
        if (abs(flux(1,1)/tiny(1.0_wp)+16.0_wp)>error_bound .or. &
            abs(tangent(1,1,1)+16.0_wp)>error_bound .or. &
            abs(dsdp(1,1,1)+16.0_wp)>error_bound) then
            error stop "Balanced Fourier reference-conductivity derivative"
        end if
    end subroutine test_reference_ratio

    subroutine test_subnormal_conductivity()
        !! Restore each complete product once, including products that amplify
        !! a subnormal conductivity into normal flux or parameter derivatives.
        type(material_t) :: material
        real(wp), parameter :: coefficient=log(1.25_wp)
        real(wp) :: x(ninput,1),flux(ndim,1),tangent(ndim,ninput,1),dsdp(ndim,1,1)
        real(wp) :: all_dsdp(ndim,nparameters,1)
        real(wp) :: smallest
        integer :: status

        smallest=nearest(0.0_wp,1.0_wp)
        call material%init(fourier_exponential_t(smallest,log(2.0_wp),0.0_wp),1,status, &
            sensitivity_parameters=[character(len=23) :: "reference_conductivity"])
        call expect(status,cf_ok)
        x(:,1)=[1.0_wp,0.0_wp,0.0_wp,1.0_wp]
        call material%eval(x,flux,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        call close(flux(1,1)/smallest,-2.0_wp)
        call close(tangent(1,1,1)/smallest,-2.0_wp)
        call close(dsdp(1,1,1),-2.0_wp)

        ! k=1.25*minsub rounds to minsub, but dq/dk0 retains the unrounded
        ! exponential. Multiplication by 1/tiny then gives q=-1.25*epsilon.
        call material%init(fourier_exponential_t(smallest,coefficient,0.0_wp),1,status, &
            sensitivity_parameters=[character(len=23) :: "reference_conductivity", &
                "temperature_coefficient","reference_temperature"])
        call expect(status,cf_ok)
        call material%eval(x,flux,tangent,status,dsdp=all_dsdp)
        call expect(status,cf_ok)
        call close(flux(1,1)/smallest,-1.0_wp)
        call close(tangent(1,1,1)/smallest,-1.0_wp)
        call close(all_dsdp(1,1,1),-1.25_wp)

        x(1,1)=1.0_wp/tiny(1.0_wp)
        call material%eval(x,flux,tangent,status,dsdp=all_dsdp)
        call expect(status,cf_ok)
        call close(flux(1,1)/epsilon(1.0_wp),-1.25_wp)
        call close(tangent(1,1,1)/smallest,-1.0_wp)
        call close(tangent(1,ninput,1)/epsilon(1.0_wp),-1.25_wp*coefficient)
        call close(all_dsdp(1,1,1)*tiny(1.0_wp),-1.25_wp)
        call close(all_dsdp(1,2,1)/epsilon(1.0_wp),-1.25_wp)
        call close(all_dsdp(1,3,1)/epsilon(1.0_wp),1.25_wp*coefficient)
    end subroutine test_subnormal_conductivity

    subroutine test_invalid_trials()
        type(material_t) :: material
        real(wp) :: x(ninput,2), flux(ndim,2), wrong(ninput-1,2)
        integer :: status

        call material%init(fourier_exponential_t(0.0_wp,0.0_wp,300.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(fourier_exponential_t(1.0_wp,0.0_wp,-1.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(fourier_exponential_t(5.0_wp,0.1_wp,300.0_wp),2,status)
        call expect(status,cf_ok)
        x(:,1) = [1.0_wp,-1.0_wp,0.5_wp,300.0_wp]
        x(:,2) = x(:,1)
        x(ninput,2) = -1.0_wp
        call material%eval(x,flux,status=status)
        call expect(status,cf_evaluation_failed)
        if (material%last_failure_point() /= 2) error stop "Bad point not flagged"
        call material%commit(status)
        call expect(status,cf_no_trial)
        x(ninput,2) = 300.0_wp
        x(1,2) = ieee_value(0.0_wp,ieee_quiet_nan)
        call material%eval(x,flux,status=status)
        call expect(status,cf_evaluation_failed)
        x(1,2) = 1.0_wp
        call material%eval(x,flux,status=status)
        call expect(status,cf_ok)
        call material%commit(status)
        call expect(status,cf_ok)
        call material%eval(wrong,flux,status=status)
        call expect(status,cf_invalid_shape)

        ! Exponentials outside the representable k(T) range must fail.
        x(ninput,2) = 300.0_wp+(log(huge(1.0_wp))+4.0_wp)/0.1_wp
        call material%eval(x,flux,status=status)
        call expect(status,cf_evaluation_failed)
        x(ninput,2) = 300.0_wp
        call material%eval(x,flux,status=status)
        call expect(status,cf_ok)
    end subroutine test_invalid_trials

end program test_fourier_exponential

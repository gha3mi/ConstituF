! SPDX-License-Identifier: Apache-2.0

program test_fourier_isotropic
    !! Test convention: absolute temperature in K, length in m and time in s.
    !! Conductivity is W/(m*K), temperature gradient K/m and heat flux W/m2.
    !!
    !! Independent heat-flux, tangent, entropy-production and API regression tests.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use constituf, only: wp, material_t, cf_ok, cf_invalid_argument, &
                        cf_invalid_shape, cf_evaluation_failed, cf_no_trial, cf_unsupported
    use constituf_fourier_isotropic, only: fourier_isotropic_t
    implicit none (type, external)

    integer, parameter :: ndim = 3, npoints = 4

    call test_schema_and_flux()
    call test_tangent_and_rotation()
    call test_sensitivity()
    call test_invalid_trials()

contains

    subroutine expect_status(value, expected)
        integer, intent(in) :: value, expected

        if (value /= expected) then
            print *, "Expected status", expected, "got", value
            error stop "Fourier material status mismatch"
        end if
    end subroutine expect_status

    subroutine expect_close(value, expected)
        real(wp), intent(in) :: value, expected

        if (.not. ieee_is_finite(value) .or. .not. ieee_is_finite(expected)) then
            error stop "Nonfinite Fourier comparison"
        end if
        if (abs(value-expected) > 128.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(expected))) then
            print *, "Expected", expected, "got", value
            error stop "Fourier material numerical mismatch"
        end if
    end subroutine expect_close

    subroutine test_schema_and_flux()
        type(material_t) :: material
        real(wp) :: gradient(ndim,npoints), flux(ndim,npoints)
        real(wp) :: tangent(ndim,ndim,npoints)
        integer :: status, ipoint, i, j

        call material%init(fourier_isotropic_t(12.0_wp), npoints, status)
        call expect_status(status, cf_ok)
        if (material%input_size() /= ndim .or. material%output_size() /= ndim) then
            error stop "Thermal dimension mismatch"
        end if
        if (material%state_size() /= 0 .or. material%state_storage_bytes() /= 0) then
            error stop "Fourier material must be stateless"
        end if
        if (.not. material%field_matches("input", "temperature.gradient", &
            "vector", "temperature_gradient", "spatial", [ndim])) then
            error stop "Temperature-gradient schema mismatch"
        end if
        if (.not. material%field_matches("output", "heat.flux", &
            "vector", "heat_flux", "spatial", [ndim])) then
            error stop "Spatial heat-flux schema mismatch"
        end if
        if (material%field_matches("output", "heat.flux", &
            "vector", "heat_flux", "reference", [ndim])) then
            error stop "Thermal frame mismatch silently accepted"
        end if

        gradient(:,1) = [1.0_wp,-2.0_wp,0.5_wp]
        gradient(:,2) = [-3.0_wp,4.0_wp,1.25_wp]
        gradient(:,3) = [0.0_wp,0.0_wp,0.0_wp]
        gradient(:,4) = [0.75_wp,-0.5_wp,0.125_wp]
        call material%eval(gradient, flux, tangent, status)
        call expect_status(status, cf_ok)
        do ipoint = 1, npoints
            do i = 1, ndim
                call expect_close(flux(i,ipoint), -12.0_wp*gradient(i,ipoint))
                do j = 1, ndim
                    if (i == j) then
                        call expect_close(tangent(i,j,ipoint), -12.0_wp)
                    else
                        call expect_close(tangent(i,j,ipoint), 0.0_wp)
                    end if
                end do
            end do
            if (dot_product(flux(:,ipoint),gradient(:,ipoint)) > 0.0_wp) then
                error stop "Nonphysical negative heat-production inequality"
            end if
        end do
        call material%commit(status)
        call expect_status(status, cf_ok)
        call material%eval(gradient, flux, status=status)
        call expect_status(status, cf_ok)
        call material%revert()
        call material%commit(status)
        call expect_status(status, cf_no_trial)
    end subroutine test_schema_and_flux

    subroutine test_tangent_and_rotation()
        type(material_t) :: batch, single
        real(wp) :: x(ndim,2), y(ndim,2), tangent(ndim,ndim,2)
        real(wp) :: positive(ndim,2), negative(ndim,2), changed(ndim,2)
        real(wp) :: alone(ndim,1)
        real(wp) :: derivative, h, error_bound
        integer :: status, i, j

        call batch%init(fourier_isotropic_t(6.0_wp), 2, status)
        call expect_status(status, cf_ok)
        call single%init(fourier_isotropic_t(6.0_wp), 1, status)
        call expect_status(status, cf_ok)
        x(:,1) = [2.0_wp,-1.0_wp,0.4_wp]
        x(:,2) = [-x(2,1),x(1,1),x(3,1)]
        call batch%eval(x, y, tangent, status)
        call expect_status(status, cf_ok)
        if (.not. all(ieee_is_finite(tangent))) error stop "Nonfinite Fourier analytical tangent"
        call single%eval(x(:,1:1), alone, status=status)
        call expect_status(status, cf_ok)
        do i = 1, ndim
            call expect_close(y(i,1), alone(i,1))
        end do
        call expect_close(y(1,2), -y(2,1))
        call expect_close(y(2,2), y(1,1))
        call expect_close(y(3,2), y(3,1))

        ! Fourier flux is affine in the gradient, so a finite central step
        ! has no truncation error and avoids cancellation in every precision.
        h = 0.25_wp*maxval(abs(x))
        do j = 1, ndim
            changed = x
            changed(j,1) = changed(j,1)+h
            call batch%eval(changed, positive, status=status)
            call expect_status(status, cf_ok)
            changed(j,1) = changed(j,1)-2.0_wp*h
            call batch%eval(changed, negative, status=status)
            call expect_status(status, cf_ok)
            if (.not. all(ieee_is_finite(positive)) .or. .not. all(ieee_is_finite(negative))) then
                error stop "Nonfinite Fourier tangent secant response"
            end if
            error_bound = 32.0_wp*epsilon(1.0_wp)*maxval(abs(positive)+abs(negative))/(2.0_wp*h)
            if (.not. ieee_is_finite(error_bound)) error stop "Nonfinite Fourier tangent error bound"
            do i = 1, ndim
                derivative = (positive(i,1)-negative(i,1))/(2.0_wp*h)
                if (.not. ieee_is_finite(derivative) .or. abs(derivative-tangent(i,j,1)) > error_bound) then
                    error stop "Fourier constitutive tangent mismatch"
                end if
            end do
        end do
    end subroutine test_tangent_and_rotation

    subroutine test_sensitivity()
        type(material_t) :: tracked, positive, negative
        real(wp) :: gradient(ndim,1), flux(ndim,1)
        real(wp) :: dsdp(ndim,1,1), tangent(ndim,ndim,1)
        real(wp) :: plus_flux(ndim,1), minus_flux(ndim,1), bad(ndim-1,1,1)
        real(wp) :: h, numerical, error_bound
        integer :: status, i

        call tracked%init(fourier_isotropic_t(10.0_wp), 1, status, &
            sensitivity_parameters=[character(len=12) :: "conductivity"])
        call expect_status(status, cf_ok)
        if (tracked%sensitivity_size() /= 1) error stop "Missing thermal sensitivity"
        if (tracked%state_storage_bytes() /= 0) error stop "Stateless derivative storage"
        gradient(:,1) = [1.2_wp,-2.3_wp,0.7_wp]
        call tracked%eval(gradient, flux, tangent, status, dsdp=dsdp)
        call expect_status(status, cf_ok)
        ! The conductivity dependence is also affine.
        h = 0.25_wp*10.0_wp
        call positive%init(fourier_isotropic_t(10.0_wp+h), 1, status)
        call expect_status(status, cf_ok)
        call negative%init(fourier_isotropic_t(10.0_wp-h), 1, status)
        call expect_status(status, cf_ok)
        call positive%eval(gradient, plus_flux, status=status)
        call expect_status(status, cf_ok)
        call negative%eval(gradient, minus_flux, status=status)
        call expect_status(status, cf_ok)
        do i = 1, ndim
            call expect_close(dsdp(i,1,1), -gradient(i,1))
            numerical = (plus_flux(i,1)-minus_flux(i,1))/(2.0_wp*h)
            error_bound = 32.0_wp*epsilon(1.0_wp)* &
                (abs(plus_flux(i,1))+abs(minus_flux(i,1)))/(2.0_wp*h)
            if (.not. all(ieee_is_finite([numerical,error_bound]))) then
                error stop "Nonfinite Fourier parameter comparison"
            end if
            if (abs(numerical-dsdp(i,1,1)) > error_bound) then
                error stop "Fourier parameter derivative mismatch"
            end if
        end do
        call tracked%commit(status)
        call expect_status(status, cf_ok)
        call tracked%eval(gradient, flux, status=status)
        call expect_status(status, cf_invalid_argument)
        call tracked%eval(gradient, flux, status=status, dsdp=bad)
        call expect_status(status, cf_invalid_shape)
        call tracked%eval(gradient, flux, status=status, dsdp=dsdp)
        call expect_status(status, cf_ok)
        call tracked%revert()
    end subroutine test_sensitivity

    subroutine test_invalid_trials()
        type(material_t) :: material
        real(wp) :: x(ndim,2), flux(ndim,2), invalid(ndim-1,2)
        real(wp) :: extra_derivative(ndim,1,2)
        integer :: status

        call material%init(fourier_isotropic_t(0.0_wp), 2, status)
        call expect_status(status, cf_invalid_argument)
        call material%init(fourier_isotropic_t(-1.0_wp), 2, status)
        call expect_status(status, cf_invalid_argument)
        call material%init(fourier_isotropic_t(ieee_value(0.0_wp,ieee_quiet_nan)), 2, status)
        call expect_status(status, cf_invalid_argument)
        call material%init(fourier_isotropic_t(4.0_wp), 2, status)
        call expect_status(status, cf_ok)
        x = 0.0_wp
        x(1,:) = 1.0_wp
        x(1,2) = ieee_value(0.0_wp,ieee_quiet_nan)
        call material%eval(x, flux, status=status)
        call expect_status(status, cf_evaluation_failed)
        if (material%last_failure_point() /= 2) error stop "Incorrect failing thermal point"
        call material%commit(status)
        call expect_status(status, cf_no_trial)

        ! Finite inputs can still produce a nonfinite flux; reject the trial
        ! and retain the same failed-point/commit contract as an invalid input.
        x(1,2) = huge(1.0_wp)
        call material%eval(x, flux, status=status)
        call expect_status(status, cf_evaluation_failed)
        if (material%last_failure_point() /= 2) error stop "Incorrect overflowing thermal point"
        call material%commit(status)
        call expect_status(status, cf_no_trial)

        x(1,2) = 1.0_wp
        call material%eval(x, flux, status=status)
        call expect_status(status, cf_ok)
        invalid = 0.0_wp
        call material%eval(invalid, flux, status=status)
        call expect_status(status, cf_invalid_shape)
        call material%init(fourier_isotropic_t(0.0_wp), 2, status)
        call expect_status(status, cf_invalid_argument)
        call material%eval(x, flux, status=status)
        call expect_status(status, cf_ok)
        call expect_close(flux(1,2), -4.0_wp)
        call material%eval(x, flux, status=status, dsdp=extra_derivative)
        call expect_status(status, cf_unsupported)
    end subroutine test_invalid_trials

end program test_fourier_isotropic

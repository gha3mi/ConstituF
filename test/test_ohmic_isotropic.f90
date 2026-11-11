! SPDX-License-Identifier: Apache-2.0

program test_ohmic_isotropic
    !! Test convention: electric field in V/m, current density in A/m2 and
    !! conductivity in S/m. Derived Joule power density is W/m3.
    !!
    !! Verify field metadata, Joule power density, rotation, tangent and failures.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
    use constituf, only: wp, material_t, cf_ok, cf_invalid_argument, &
                        cf_invalid_shape, cf_evaluation_failed, cf_no_trial, cf_unsupported
    use constituf_ohmic_isotropic, only: ohmic_isotropic_t
    implicit none (type, external)

    integer, parameter :: ndim = 3

    call test_reference_and_schema()
    call test_tangent_rotation_batch()
    call test_sensitivity()
    call test_invalid_and_recovery()

contains

    subroutine expect(actual, expected)
        integer, intent(in) :: actual, expected

        if (actual /= expected) then
            print *, "actual status", actual, "expected", expected
            error stop "Ohmic status mismatch"
        end if
    end subroutine expect

    subroutine close(actual, expected)
        real(wp), intent(in) :: actual, expected

        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) then
            error stop "Nonfinite Ohmic comparison"
        end if
        if (abs(actual-expected) > 256.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(expected))) then
            print *, "actual", actual, "expected", expected
            error stop "Ohmic reference mismatch"
        end if
    end subroutine close

    subroutine test_reference_and_schema()
        type(material_t) :: material
        real(wp) :: electric(ndim,2), current_density(ndim,2), tangent(ndim,ndim,2)
        integer :: status, i, j, ipoint

        call material%init(ohmic_isotropic_t(2.5_wp),2,status)
        call expect(status,cf_ok)
        if (material%input_size() /= ndim .or. material%output_size() /= ndim) then
            error stop "Ohmic input/output sizes"
        end if
        if (material%state_storage_bytes() /= 0) error stop "Nonzero ohmic state"
        if (.not. material%field_matches("input", "electric.field", &
            "vector", "electric_field", "spatial", [ndim])) then
            error stop "Ohmic electric-field metadata"
        end if
        if (.not. material%field_matches("output", "electric.current_density", &
            "vector", "electric_current_density", "spatial", [ndim])) then
            error stop "Ohmic current_density metadata"
        end if
        if (material%field_matches("output", "electric.current_density", &
            "vector", "electric_current_density", "reference", [ndim])) then
            error stop "Ohmic incorrect frame accepted"
        end if

        electric(:,1) = [1.2_wp,-2.0_wp,0.5_wp]
        electric(:,2) = 0.0_wp
        call material%eval(electric,current_density,tangent,status)
        call expect(status,cf_ok)
        do ipoint = 1, 2
            do i = 1, ndim
                call close(current_density(i,ipoint),2.5_wp*electric(i,ipoint))
                do j = 1, ndim
                    if (i == j) then
                        call close(tangent(i,j,ipoint),2.5_wp)
                    else
                        call close(tangent(i,j,ipoint),0.0_wp)
                    end if
                end do
            end do
            if (dot_product(electric(:,ipoint),current_density(:,ipoint)) < 0.0_wp) then
                error stop "Negative Joule power density"
            end if
        end do
        call material%commit(status)
        call expect(status,cf_ok)
        call material%eval(electric,current_density,status=status)
        call expect(status,cf_ok)
        call material%revert()
        call material%commit(status)
        call expect(status,cf_no_trial)
    end subroutine test_reference_and_schema

    subroutine test_tangent_rotation_batch()
        type(material_t) :: batch, single
        real(wp) :: x(ndim,2), y(ndim,2), tangent(ndim,ndim,2)
        real(wp) :: alone(ndim,1), trial(ndim,2), plus(ndim,2), minus(ndim,2)
        real(wp) :: h, derivative, error_bound
        integer :: status, i, j

        call batch%init(ohmic_isotropic_t(3.0_wp),2,status)
        call expect(status,cf_ok)
        call single%init(ohmic_isotropic_t(3.0_wp),1,status)
        call expect(status,cf_ok)
        x(:,1) = [1.2_wp,-0.7_wp,2.0_wp]
        x(:,2) = [-x(2,1),x(1,1),x(3,1)]
        call batch%eval(x,y,tangent,status)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(tangent))) error stop "Nonfinite Ohmic analytical tangent"
        call single%eval(x(:,1:1),alone,status=status)
        call expect(status,cf_ok)
        do i = 1, ndim
            call close(y(i,1),alone(i,1))
        end do
        call close(y(1,2),-y(2,1))
        call close(y(2,2),y(1,1))
        call close(y(3,2),y(3,1))

        ! Current is affine in the field. A finite central step has no
        ! truncation error and avoids cancellation at every supported kind.
        h = 0.25_wp*maxval(abs(x))
        do j = 1, ndim
            trial = x
            trial(j,1) = trial(j,1)+h
            call batch%eval(trial,plus,status=status)
            call expect(status,cf_ok)
            trial(j,1) = trial(j,1)-2.0_wp*h
            call batch%eval(trial,minus,status=status)
            call expect(status,cf_ok)
            if (.not. all(ieee_is_finite(plus)) .or. .not. all(ieee_is_finite(minus))) then
                error stop "Nonfinite Ohmic tangent secant response"
            end if
            error_bound = 32.0_wp*epsilon(1.0_wp)*maxval(abs(plus)+abs(minus))/(2.0_wp*h)
            if (.not. ieee_is_finite(error_bound)) error stop "Nonfinite Ohmic tangent error bound"
            do i = 1, ndim
                derivative = (plus(i,1)-minus(i,1))/(2.0_wp*h)
                if (.not. ieee_is_finite(derivative) .or. abs(derivative-tangent(i,j,1)) > error_bound) then
                    error stop "Ohmic tangent mismatch"
                end if
            end do
        end do
    end subroutine test_tangent_rotation_batch

    subroutine test_sensitivity()
        type(material_t) :: tracked, positive, negative
        real(wp) :: electric(ndim,1), current_density(ndim,1), derivative(ndim,1,1)
        real(wp) :: plus(ndim,1), minus(ndim,1)
        real(wp) :: h,numerical,error_bound
        integer :: status, i

        call tracked%init(ohmic_isotropic_t(2.0_wp),1,status, &
            sensitivity_parameters=[character(len=12) :: "conductivity"])
        call expect(status,cf_ok)
        if (tracked%sensitivity_size() /= 1) error stop "Ohmic sensitivity count"
        if (tracked%state_storage_bytes() /= 0) error stop "Ohmic derivative state"
        electric(:,1) = [0.2_wp,-1.0_wp,2.3_wp]
        call tracked%eval(electric,current_density,status=status,dsdp=derivative)
        call expect(status,cf_ok)
        ! Conductivity also enters the constitutive response affinely.
        h = 0.25_wp*2.0_wp
        call positive%init(ohmic_isotropic_t(2.0_wp+h),1,status)
        call expect(status,cf_ok)
        call negative%init(ohmic_isotropic_t(2.0_wp-h),1,status)
        call expect(status,cf_ok)
        call positive%eval(electric,plus,status=status)
        call expect(status,cf_ok)
        call negative%eval(electric,minus,status=status)
        call expect(status,cf_ok)
        if (.not. all(ieee_is_finite(plus)) .or. .not. all(ieee_is_finite(minus))) then
            error stop "Nonfinite Ohmic parameter secant response"
        end if
        do i = 1, ndim
            call close(derivative(i,1,1),electric(i,1))
            numerical=(plus(i,1)-minus(i,1))/(2.0_wp*h)
            error_bound = 32.0_wp*epsilon(1.0_wp)* &
                (abs(plus(i,1))+abs(minus(i,1)))/(2.0_wp*h)
            if (.not. all(ieee_is_finite([numerical,error_bound]))) then
                error stop "Nonfinite Ohmic parameter comparison"
            end if
            if (abs(numerical-derivative(i,1,1)) > error_bound) then
                error stop "Ohmic parameter sensitivity"
            end if
        end do
        call tracked%commit(status)
        call expect(status,cf_ok)
        call tracked%eval(electric,current_density,status=status)
        call expect(status,cf_invalid_argument)
    end subroutine test_sensitivity

    subroutine test_invalid_and_recovery()
        type(material_t) :: material
        real(wp) :: electric(ndim,2), current_density(ndim,2), wrong(ndim-1,2)
        real(wp) :: extra(ndim,1,2)
        integer :: status

        call material%init(ohmic_isotropic_t(0.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(ohmic_isotropic_t(-1.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(ohmic_isotropic_t(ieee_value(0.0_wp,ieee_quiet_nan)),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(ohmic_isotropic_t(ieee_value(0.0_wp,ieee_positive_inf)),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(ohmic_isotropic_t(4.0_wp),2,status)
        call expect(status,cf_ok)
        electric = 0.0_wp
        electric(1,:) = 1.0_wp
        electric(2,2) = ieee_value(0.0_wp,ieee_quiet_nan)
        call material%eval(electric,current_density,status=status)
        call expect(status,cf_evaluation_failed)
        if (material%last_failure_point() /= 2) error stop "Ohmic failed point"
        call material%commit(status)
        call expect(status,cf_no_trial)
        electric(2,2) = 2.0_wp
        call material%eval(electric,current_density,status=status)
        call expect(status,cf_ok)
        call close(current_density(1,2),4.0_wp)
        call material%eval(wrong,current_density,status=status)
        call expect(status,cf_invalid_shape)
        call material%eval(electric,current_density,status=status,dsdp=extra)
        call expect(status,cf_unsupported)
        electric(1,2)=huge(1.0_wp)/2.0_wp
        call material%eval(electric,current_density,status=status)
        call expect(status,cf_evaluation_failed)
        if (material%last_failure_point()/=2) error stop "Ohmic overflowing point"
        call material%commit(status)
        call expect(status,cf_no_trial)
        electric(1,2)=1.0_wp
        call material%eval(electric,current_density,status=status)
        call expect(status,cf_ok)
    end subroutine test_invalid_and_recovery

end program test_ohmic_isotropic

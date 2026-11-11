! SPDX-License-Identifier: Apache-2.0

program test_core_typed_state
    !! Mixed-state ownership, state swaps, failure atomicity and metadata.
    !! Synthetic real fields are dimensionless. Discrete fields are counts and flags.
    use constituf, only: wp, ik, material_t, cf_ok, cf_invalid_argument, &
                        cf_no_trial, cf_unsupported, cf_evaluation_failed
    use constituf_test_coupled, only: coupled_test_model_t
    use constituf_test_typed, only: typed_test_model_t, missing_typed_test_model_t
    implicit none (type, external)

    type(material_t) :: mat
    real(wp) :: x(1,3), y(1,3), tangent(1,1,3)
    integer :: status
    integer :: empty_shape(0)

    call mat%init(typed_test_model_t(), 3, status)
    call check(status, cf_ok)
    if (mat%state_size() /= 1 .or. mat%integer_state_size() /= 1 .or. &
        mat%logical_state_size() /= 1) error stop "Unexpected typed state sizes"
    if (mat%state_storage_bytes() /= 6_ik*int(storage_size(0.0_wp)/8, ik) + &
        6_ik*int(storage_size(0_ik)/8, ik) + &
        6_ik*int(storage_size(.false.)/8, ik)) then
        error stop "Typed state memory accounting"
    end if
    if (.not. mat%field_matches("integer_state", "positive.count", "scalar", "event_count", "none", empty_shape)) then
        error stop "Missing integer metadata"
    end if
    if (.not. mat%field_matches("logical_state", "positive.seen", "scalar", "event_flag", "none", empty_shape)) then
        error stop "Missing logical metadata"
    end if
    if (mat%field_matches("state", "positive.count", &
                          "scalar", "event_count", "none", empty_shape)) then
        error stop "Integer state mistaken for real"
    end if

    x(1,:) = [1.0_wp, -0.5_wp, 0.0_wp]
    call mat%eval(x, y, tangent, status)
    call check(status, cf_ok)
    call check_response(y(1,:), [13.0_wp, 0.5_wp, 1.0_wp])
    if (any(abs(tangent - 1.0_wp) > 1.0e-12_wp)) error stop "Unexpected tangent"
    call mat%revert()
    call mat%commit(status)
    call check(status, cf_no_trial)
    call mat%eval(x, y, status=status)
    call check(status, cf_ok)
    call check_response(y(1,:), [13.0_wp, 0.5_wp, 1.0_wp])
    call mat%commit(status)
    call check(status, cf_ok)

    x(1,:) = [-0.1_wp, 2.0_wp, -1.0_wp]
    call mat%eval(x, y, status=status)
    call check(status, cf_ok)
    call check_response(y(1,:), [12.9_wp, 13.5_wp, 0.0_wp])
    call mat%commit(status)
    call check(status, cf_ok)

    ! A failure in one point invalidates the entire batch trial.
    x(1,:) = [1.0_wp, -10.0_wp, 3.0_wp]
    call mat%eval(x, y, status=status)
    call check(status, cf_evaluation_failed)
    if (mat%last_failure_point() /= 2) error stop "Failed point not reported"
    call mat%commit(status)
    call check(status, cf_no_trial)
    x = 0.0_wp
    call mat%eval(x, y, status=status)
    call check(status, cf_ok)
    call check_response(y(1,:), [12.9_wp, 13.5_wp, 0.0_wp])
    call mat%revert()
    call mat%eval(x, y, status=status)
    call check(status, cf_ok)
    call check_response(y(1,:), [12.9_wp, 13.5_wp, 0.0_wp])

    ! Invalid metadata must preserve a previously initialized handle.
    call mat%init(typed_test_model_t(duplicate_history=.true.), 3, status)
    call check(status, cf_invalid_argument)
    if (mat%integer_state_size() /= 1) error stop "Failed init discarded history"
    call mat%eval(x, y, status=status)
    call check(status, cf_ok)
    call check_response(y(1,:), [12.9_wp, 13.5_wp, 0.0_wp])

    ! Reinitialization to stateless law releases all typed buffers.
    call mat%init(coupled_test_model_t(), 3, status)
    call check(status, cf_ok)
    if (mat%state_storage_bytes() /= 0) error stop "Unused typed state retained"
    if (mat%integer_state_size() /= 0 .or. mat%logical_state_size() /= 0) then
        error stop "Typed metadata retained after reinitialization"
    end if
    block
        real(wp) :: coupled_input(2,3), coupled_response(3,3)

        coupled_input(1,:) = 2.0_wp
        coupled_input(2,:) = 0.0_wp
        call mat%eval(coupled_input, coupled_response, status=status)
        call check(status, cf_ok)
        call check_response(coupled_response(1,:), [4.0_wp, 4.0_wp, 4.0_wp])
        call check_response(coupled_response(2,:), [6.0_wp, 6.0_wp, 6.0_wp])
        call check_response(coupled_response(3,:), [6.0_wp, 6.0_wp, 6.0_wp])
    end block

    ! Missing typed method must report unsupported, not call real-only method.
    call mat%init(missing_typed_test_model_t(), 3, status)
    call check(status, cf_ok)
    if (mat%integer_state_size() /= 1 .or. mat%logical_state_size() /= 0) then
        error stop "Integer-only state layout"
    end if
    call mat%eval(x, y, status=status)
    call check(status, cf_unsupported)
    call mat%commit(status)
    call check(status, cf_no_trial)

contains

    subroutine check(actual, expected)
        integer, intent(in) :: actual, expected

        if (actual /= expected) then
            print *, "Expected", expected, "got", actual
            error stop "Unexpected status"
        end if
    end subroutine check

    subroutine check_response(actual, expected)
        real(wp), intent(in) :: actual(:), expected(:)

        if (any(abs(actual - expected) > 1.0e-12_wp)) then
            print *, "Actual", actual, "expected", expected
            error stop "Incorrect mixed-state response"
        end if
    end subroutine check_response
end program test_core_typed_state

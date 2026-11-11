! SPDX-License-Identifier: Apache-2.0

program test_core_material_api
    !! Core material-handle dispatch, field metadata and transactional safety.
    !! All synthetic model coefficients, fields and real history are dimensionless.
    use constituf, only: wp, material_t, cf_ok, cf_invalid_argument, cf_invalid_shape, &
                        cf_no_trial, cf_unsupported, cf_evaluation_failed, &
                        cf_uninitialized
    use constituf_test_coupled, only: coupled_test_model_t
    use constituf_test_failure, only: failure_test_model_t
    implicit none (type, external)

    call test_coupled_api()
    call test_transactional_failure()

contains

    subroutine expect_status(actual, expected, context)
        integer, intent(in) :: actual, expected
        character(*), intent(in) :: context

        if (actual /= expected) then
            print *, "Expected status ", expected, " got ", actual, ": ", context
            error stop "Unexpected material API status"
        end if
    end subroutine expect_status

    subroutine expect_close(actual, expected, context)
        real(wp), intent(in) :: actual, expected
        character(*), intent(in) :: context
        real(wp) :: tolerance

        tolerance = 128.0_wp * epsilon(1.0_wp) * &
                    max(1.0_wp, abs(actual), abs(expected))
        if (abs(actual - expected) > tolerance) then
            print *, "Expected ", expected, " got ", actual, ": ", context
            error stop "Numerical mismatch"
        end if
    end subroutine expect_close

    subroutine test_coupled_api()
        type(material_t) :: mat
        real(wp) :: x(2,3), y(3,3), tangent(3,2,3), wrong(1,3)
        integer :: ierr, j
        integer :: empty_shape(0)
        integer, allocatable :: shape(:)
        character(:), allocatable :: convention, measure, frame
        logical :: found

        x = 0.0_wp
        call mat%eval(x, y, status=ierr)
        call expect_status(ierr, cf_uninitialized, "uninitialized evaluation")
        call mat%init(coupled_test_model_t(), npoints=3, status=ierr)
        call expect_status(ierr, cf_ok, "coupled initialization")
        if (mat%input_size() /= 2 .or. mat%output_size() /= 3) then
            error stop "Coupled dimensions"
        end if
        if (mat%state_size() /= 0) error stop "Unexpected stateless history"
        if (mat%input_offset("field.x1") /= 1) error stop "First input offset"
        if (mat%input_offset("field.x2") /= 2) error stop "Second input offset"
        if (mat%input_offset("missing") /= 0) error stop "Unknown input offset"
        if (mat%output_offset("response.y2") /= 2) error stop "Second output offset"
        if (mat%output_offset("response.y3") /= 3) error stop "Third output offset"
        if (mat%input_components("field.x1") /= 1) error stop "Input component count"
        if (mat%output_components("response.y2") /= 1) error stop "Output component count"
        if (.not. mat%field_matches("input", "field.x2", "scalar", &
                                    "test_scalar", "none", empty_shape)) then
            error stop "Input binding mismatch"
        end if
        if (mat%field_matches("input", "field.x2", "scalar", &
                              "test_scalar", "spatial", empty_shape)) then
            error stop "Incorrect reference frame accepted"
        end if
        call mat%field_info("input", "field.x2", convention, &
                            measure, frame, shape, found)
        if (.not. found) error stop "Missing field metadata"
        if (convention /= "scalar" .or. &
            measure /= "test_scalar" .or. frame /= "none") then
            error stop "Incorrect field metadata"
        end if
        if (size(shape) /= 0) error stop "Scalar field rank"

        do j = 1, 3
            x(:,j) = [real(j,wp), -2.0_wp]
        end do
        call mat%eval(x, y, tangent, ierr)
        call expect_status(ierr, cf_ok, "coupled evaluation")
        do j = 1, 3
            call expect_close(y(1,j), 2.0_wp*x(1,j) + 3.0_wp*x(2,j), "coupled y1")
            call expect_close(y(2,j), 3.0_wp*x(1,j) + 4.0_wp*x(2,j), "coupled y2")
            call expect_close(y(3,j), 3.0_wp*(x(1,j)-2.0_wp*x(2,j)), "coupled y3")
            call expect_close(tangent(1,1,j), 2.0_wp, "J11")
            call expect_close(tangent(1,2,j), 3.0_wp, "J12")
            call expect_close(tangent(2,1,j), 3.0_wp, "J21")
            call expect_close(tangent(2,2,j), 4.0_wp, "J22")
            call expect_close(tangent(3,1,j), 3.0_wp, "J31")
            call expect_close(tangent(3,2,j), -6.0_wp, "J32")
        end do
        if (mat%last_failure_point() /= 0) error stop "Unexpected failed point"
        call mat%commit(ierr)
        call expect_status(ierr, cf_ok, "accepted trial")
        call mat%commit(ierr)
        call expect_status(ierr, cf_no_trial, "double commit")

        wrong = 0.0_wp
        call mat%eval(wrong, y, status=ierr)
        call expect_status(ierr, cf_invalid_shape, "input shape")
        call mat%commit(ierr)
        call expect_status(ierr, cf_no_trial, "invalid trial cannot commit")
        call mat%init(coupled_test_model_t(), npoints=0, status=ierr)
        call expect_status(ierr, cf_invalid_argument, "invalid point count")
        call mat%init(coupled_test_model_t(), npoints=3, status=ierr, initial_input=wrong)
        call expect_status(ierr, cf_invalid_shape, "initial input shape")
        call mat%init(failure_test_model_t(duplicate_fields=.true.), npoints=3, status=ierr)
        call expect_status(ierr, cf_invalid_argument, "duplicate field names")
        call mat%eval(x, y, status=ierr)
        call expect_status(ierr, cf_ok, "failed initialization preserves model")
        call expect_close(y(1,1), 2.0_wp*x(1,1)+3.0_wp*x(2,1), "previous model")
        call mat%revert()
        call mat%commit(ierr)
        call expect_status(ierr, cf_no_trial, "revert invalidates trial")
    end subroutine test_coupled_api

    subroutine test_transactional_failure()
        type(material_t) :: mat
        real(wp) :: x(1,3), y(1,3), tangent(1,1,3)
        integer :: ierr

        call mat%init(failure_test_model_t(duplicate_fields=.true.), npoints=3, status=ierr)
        call expect_status(ierr, cf_invalid_argument, "duplicate input rejected")
        call mat%init(failure_test_model_t(), npoints=3, status=ierr)
        call expect_status(ierr, cf_ok, "failure test model init")
        x = 1.0_wp
        call mat%eval(x, y, tangent, ierr)
        call expect_status(ierr, cf_unsupported, "unsupported tangent")
        call mat%eval(x, y, status=ierr)
        call expect_status(ierr, cf_ok, "first accepted trial")
        call mat%commit(ierr)
        call expect_status(ierr, cf_ok, "failure test model commit")

        x = 2.0_wp
        x(1,2) = -1.0_wp
        call mat%eval(x, y, status=ierr)
        call expect_status(ierr, cf_evaluation_failed, "per-point failure")
        if (mat%last_failure_point() /= 2) error stop "Wrong failed point"
        call mat%commit(ierr)
        call expect_status(ierr, cf_no_trial, "failed batch cannot commit")

        x = 2.0_wp
        call mat%eval(x, y, status=ierr)
        call expect_status(ierr, cf_ok, "retry after failed batch")
        call expect_close(y(1,1), 3.0_wp, "history not mutated by failed batch")
        call expect_close(y(1,2), 3.0_wp, "failed point history preserved")

        ! A material may derive its initial history from nonzero input fields.
        x = 5.0_wp
        call mat%init(failure_test_model_t(), npoints=3, status=ierr, initial_input=x)
        call expect_status(ierr, cf_ok, "custom state initialization")
        x = 2.0_wp
        call mat%eval(x, y, status=ierr)
        call expect_status(ierr, cf_ok, "custom state evaluation")
        call expect_close(y(1,1), 7.0_wp, "custom initial history applied")
    end subroutine test_transactional_failure
end program test_core_material_api

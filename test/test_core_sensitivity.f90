! SPDX-License-Identifier: Apache-2.0

program test_core_sensitivity
    !! Generic sensitivity selection, buffer validation and nonsquare dispatch.
    !! Physical and path-dependent derivatives are tested per material.
    !! All synthetic model coefficients, fields and their derivatives are dimensionless.
    use constituf, only: wp, material_t, cf_ok, cf_invalid_argument, &
                        cf_invalid_shape, cf_no_trial, cf_unsupported
    use constituf_test_coupled, only: coupled_test_model_t
    use constituf_test_failure, only: failure_test_model_t
    implicit none (type, external)

    call test_parameter_rejection()
    call test_coupled_parameter()

contains

    subroutine expect_status(actual, expected)
        integer, intent(in) :: actual, expected

        if (actual /= expected) then
            print *, "Expected status", expected, "but got", actual
            error stop "Invalid sensitivity status"
        end if
    end subroutine expect_status

    subroutine test_parameter_rejection()
        type(material_t) :: mat
        real(wp) :: x(2,1), y(3,1), dsdp(3,1,1)
        integer :: status
        character(len=16) :: empty_parameters(0)

        x(:,1) = [2.0_wp, -0.5_wp]
        call mat%init(coupled_test_model_t(), 1, status)
        call expect_status(status, cf_ok)
        if (mat%sensitivity_size() /= 0) error stop "Unused sensitivity allocation"
        call mat%eval(x, y, status=status, dsdp=dsdp)
        call expect_status(status, cf_unsupported)
        call mat%init(coupled_test_model_t(), 1, status, &
                      sensitivity_parameters=[character(len=16) :: "unknown"])
        call expect_status(status, cf_invalid_argument)
        call mat%init(coupled_test_model_t(), 1, status, &
                      sensitivity_parameters=[character(len=16) :: "coupling", "coupling"])
        call expect_status(status, cf_invalid_argument)
        call mat%init(coupled_test_model_t(), 1, status, &
                      sensitivity_parameters=empty_parameters)
        call expect_status(status, cf_invalid_argument)
        if (mat%sensitivity_size() /= 0) error stop "Failed init replaced existing model"
        call mat%eval(x, y, status=status)
        call expect_status(status, cf_ok)
        if (any(abs(y(:,1)-[2.5_wp, 4.0_wp, 9.0_wp]) > 1.0e-13_wp)) then
            error stop "Failed sensitivity init changed response"
        end if
        call mat%init(failure_test_model_t(), 1, status, &
                      sensitivity_parameters=[character(len=16) :: "not_supported"])
        call expect_status(status, cf_unsupported)
    end subroutine test_parameter_rejection

    subroutine test_coupled_parameter()
        type(material_t) :: mat
        real(wp) :: x(2,2), y(3,2), tangent(3,2,2), dsdp(3,1,2), wrong(3,2,2)
        integer :: status

        call mat%init(coupled_test_model_t(), 2, status, &
                      sensitivity_parameters=[character(len=16) :: "coupling"])
        call expect_status(status, cf_ok)
        x(1,:) = [0.5_wp, 0.9_wp]
        x(2,:) = [-0.2_wp, 0.4_wp]
        call mat%eval(x, y, tangent, status, dsdp=dsdp)
        call expect_status(status, cf_ok)
        if (any(abs(dsdp(1,1,:)-x(2,:)) > 1.0e-13_wp)) then
            error stop "Coupled first field derivative"
        end if
        if (any(abs(dsdp(2,1,:)-x(1,:)) > 1.0e-13_wp)) then
            error stop "Coupled second field derivative"
        end if
        if (any(abs(dsdp(3,1,:)-(x(1,:)-2.0_wp*x(2,:))) > 1.0e-13_wp)) then
            error stop "Coupled third field derivative"
        end if
        if (any(abs(y(3,:)-3.0_wp*(x(1,:)-2.0_wp*x(2,:))) > 1.0e-13_wp)) then
            error stop "Coupled third field response"
        end if
        if (any(abs(tangent(3,1,:)-3.0_wp) > 1.0e-13_wp) .or. &
            any(abs(tangent(3,2,:)+6.0_wp) > 1.0e-13_wp)) then
            error stop "Coupled non-square tangent"
        end if
        call mat%commit(status)
        call expect_status(status, cf_ok)

        ! All tracked trials must provide correctly shaped derivatives.
        call mat%eval(x, y, status=status)
        call expect_status(status, cf_invalid_argument)
        call mat%commit(status)
        call expect_status(status, cf_no_trial)
        call mat%eval(x, y, status=status, dsdp=wrong)
        call expect_status(status, cf_invalid_shape)
        call mat%eval(x, y, status=status, dsdp=dsdp)
        call expect_status(status, cf_ok)
        call mat%revert()
        call mat%commit(status)
        call expect_status(status, cf_no_trial)
    end subroutine test_coupled_parameter

end program test_core_sensitivity

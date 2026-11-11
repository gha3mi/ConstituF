! SPDX-License-Identifier: Apache-2.0

module constituf_test_failure
    !! Test-only faulty model for verifying transactional evaluation failures.
    !! All synthetic inputs, outputs and real history are dimensionless.
    use constituf_author, only: wp, material_model_t, material_spec_t, &
                               cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    public :: failure_test_model_t

    type, extends(material_model_t) :: failure_test_model_t
        logical :: duplicate_fields = .false.
    contains
        procedure :: describe => describe_failure
        procedure :: eval_batch => evaluate_failure
        procedure :: initialize_state => initialize_failure_state
    end type failure_test_model_t

contains

    subroutine describe_failure(self, spec, status)
        class(failure_test_model_t), intent(in) :: self
        type(material_spec_t), intent(out) :: spec
        integer, intent(out) :: status

        call spec%add_input("x", 1, "scalar", "test_scalar", "none")
        if (self%duplicate_fields) call spec%add_input("x", 1, "scalar", "test_scalar", "none")
        call spec%add_output("y", 1, "scalar", "test_scalar", "none")
        call spec%add_state("history", 1, "scalar", "test_scalar", "none")
        status = cf_ok
    end subroutine describe_failure

    subroutine initialize_failure_state(self, initial_input, state, status)
        class(failure_test_model_t), intent(in) :: self
        real(wp), intent(in) :: initial_input(:,:)
        real(wp), intent(out) :: state(:,:)
        integer, intent(out) :: status

        state(1,:) = initial_input(1,:)
        status = cf_ok
    end subroutine initialize_failure_state

    subroutine evaluate_failure(self, x_old, x_new, state_old, state_trial, &
                                response, dt, tangent, point_status)
        class(failure_test_model_t), intent(in) :: self
        real(wp), intent(in) :: x_old(:,:), x_new(:,:), state_old(:,:)
        real(wp), intent(out) :: state_trial(:,:), response(:,:)
        real(wp), intent(in) :: dt
        real(wp), intent(out), optional :: tangent(:,:,:)
        integer, intent(inout) :: point_status(:)

        state_trial(1,:) = state_old(1,:) + x_new(1,:)
        response(1,:) = state_trial(1,:)
        point_status = cf_ok
        where (x_new(1,:) < 0.0_wp) point_status = cf_invalid_argument
    end subroutine evaluate_failure
end module constituf_test_failure

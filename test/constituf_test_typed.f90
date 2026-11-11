! SPDX-License-Identifier: Apache-2.0

module constituf_test_typed
    !! Test-only mixed real, integer and logical material-point state.
    !! This synthetic switch is not a physical constitutive model.
    !! Synthetic real fields are dimensionless. Discrete fields are counts and flags.
    use constituf_author, only: wp, ik, material_model_t, material_spec_t, &
                               cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    public :: typed_test_model_t, missing_typed_test_model_t

    type, extends(material_model_t) :: typed_test_model_t
        logical :: duplicate_history = .false.
    contains
        procedure :: describe => describe_typed
        procedure :: eval_batch => evaluate_legacy
        procedure :: eval_typed_batch => evaluate_typed
        procedure :: initialize_typed_state => initialize_typed
    end type typed_test_model_t

    type, extends(material_model_t) :: missing_typed_test_model_t
    contains
        procedure :: describe => describe_missing
        procedure :: eval_batch => evaluate_missing
    end type missing_typed_test_model_t

contains

    subroutine describe_typed(self, spec, status)
        class(typed_test_model_t), intent(in) :: self
        type(material_spec_t), intent(out) :: spec
        integer, intent(out) :: status

        call spec%add_input("drive", 1, "scalar", "driving_field", "none")
        call spec%add_output("response", 1, "scalar", "test_response", "none")
        call spec%add_state("integral", 1, "scalar", "history", "none")
        call spec%add_integer_state("positive.count", 1, "scalar", "event_count", "none")
        call spec%add_logical_state("positive.seen", 1, "scalar", "event_flag", "none")
        if (self%duplicate_history) then
            call spec%add_logical_state("positive.count", 1, &
                                        "scalar", "event_flag", "none")
        end if
        call spec%set_has_tangent(.true.)
        status = cf_ok
    end subroutine describe_typed

    subroutine initialize_typed(self, initial_input, integer_state, logical_state, status)
        class(typed_test_model_t), intent(in) :: self
        real(wp), intent(in) :: initial_input(:,:)
        integer(ik), intent(out) :: integer_state(:,:)
        logical, intent(out) :: logical_state(:,:)
        integer, intent(out) :: status

        integer_state = 1_ik
        logical_state = .false.
        status = cf_ok
    end subroutine initialize_typed

    subroutine evaluate_typed(self, x_old, x_new, state_old, state_trial, &
                              integer_old, integer_trial, logical_old, logical_trial, &
                              response, dt, tangent, point_status)
        class(typed_test_model_t), intent(in) :: self
        real(wp), intent(in) :: x_old(:,:), x_new(:,:), state_old(:,:)
        real(wp), intent(out) :: state_trial(:,:), response(:,:)
        integer(ik), intent(in) :: integer_old(:,:)
        integer(ik), intent(out) :: integer_trial(:,:)
        logical, intent(in) :: logical_old(:,:)
        logical, intent(out) :: logical_trial(:,:)
        real(wp), intent(in) :: dt
        real(wp), intent(out), optional :: tangent(:,:,:)
        integer, intent(inout) :: point_status(:)
        integer :: i

        do i = 1, size(x_new, 2)
            state_trial(1,i) = state_old(1,i) + x_new(1,i)
            integer_trial(1,i) = integer_old(1,i)
            if (x_new(1,i) > 0.0_wp) integer_trial(1,i) = integer_trial(1,i) + 1_ik
            logical_trial(1,i) = logical_old(1,i) .or. (x_new(1,i) > 0.0_wp)
            response(1,i) = state_trial(1,i) + real(integer_trial(1,i), wp)
            if (logical_trial(1,i)) response(1,i) = response(1,i) + 10.0_wp
            point_status(i) = cf_ok
            if (x_new(1,i) < -5.0_wp) point_status(i) = cf_invalid_argument
        end do
        if (present(tangent)) tangent(1,1,:) = 1.0_wp
    end subroutine evaluate_typed

    subroutine evaluate_legacy(self, x_old, x_new, state_old, state_trial, &
                               response, dt, tangent, point_status)
        !! This must not be called by the core for models with typed states.
        class(typed_test_model_t), intent(in) :: self
        real(wp), intent(in) :: x_old(:,:), x_new(:,:), state_old(:,:)
        real(wp), intent(out) :: state_trial(:,:), response(:,:)
        real(wp), intent(in) :: dt
        real(wp), intent(out), optional :: tangent(:,:,:)
        integer, intent(inout) :: point_status(:)

        state_trial = 0.0_wp
        response = 0.0_wp
        if (present(tangent)) tangent = 0.0_wp
        point_status = cf_invalid_argument
    end subroutine evaluate_legacy

    subroutine describe_missing(self, spec, status)
        class(missing_typed_test_model_t), intent(in) :: self
        type(material_spec_t), intent(out) :: spec
        integer, intent(out) :: status

        call spec%add_input("x", 1, "scalar", "test_input", "none")
        call spec%add_output("y", 1, "scalar", "test_output", "none")
        call spec%add_integer_state("flag", 1, "scalar", "event_count", "none")
        status = cf_ok
    end subroutine describe_missing

    subroutine evaluate_missing(self, x_old, x_new, state_old, state_trial, &
                                response, dt, tangent, point_status)
        class(missing_typed_test_model_t), intent(in) :: self
        real(wp), intent(in) :: x_old(:,:), x_new(:,:), state_old(:,:)
        real(wp), intent(out) :: state_trial(:,:), response(:,:)
        real(wp), intent(in) :: dt
        real(wp), intent(out), optional :: tangent(:,:,:)
        integer, intent(inout) :: point_status(:)

        state_trial = 0.0_wp
        response = 0.0_wp
        if (present(tangent)) tangent = 0.0_wp
        point_status = cf_invalid_argument
    end subroutine evaluate_missing
end module constituf_test_typed

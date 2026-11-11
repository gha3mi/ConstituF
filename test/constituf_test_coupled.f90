! SPDX-License-Identifier: Apache-2.0

module constituf_test_coupled
    !! Test-only coupled linear operator, not a physical material model.
    !!
    !! This test model exercises multiple field inputs, multiple outputs,
    !! off-diagonal Jacobian blocks and a non-square (3 by 2) Jacobian:
    !! \[y_1=a x_1+b x_2,\quad y_2=b x_1+c x_2,\quad
    !! y_3=b(x_1-2x_2).\]
    !! No scientific constitutive claims are made for the synthetic operator.
    !! All synthetic coefficients, inputs and outputs are dimensionless.
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok
    implicit none (type, external)
    private

    public :: coupled_test_model_t

    type, extends(material_model_t) :: coupled_test_model_t
        real(wp) :: a = 2.0_wp
        real(wp) :: b = 3.0_wp
        real(wp) :: c = 4.0_wp
    contains
        procedure :: describe => describe_coupled
        procedure :: eval_batch => evaluate_coupled
        procedure :: eval_sensitivity_batch => sensitivity_coupled
    end type coupled_test_model_t

contains

    subroutine describe_coupled(self, spec, status)
        class(coupled_test_model_t), intent(in) :: self
        type(material_spec_t), intent(out) :: spec
        integer, intent(out) :: status

        call spec%add_input("field.x1", 1, "scalar", "test_scalar", "none")
        call spec%add_input("field.x2", 1, "scalar", "test_scalar", "none")
        call spec%add_output("response.y1", 1, "scalar", "test_scalar", "none")
        call spec%add_output("response.y2", 1, "scalar", "test_scalar", "none")
        call spec%add_output("response.y3", 1, "scalar", "test_scalar", "none")
        call spec%add_parameter("coupling")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_coupled

    subroutine evaluate_coupled(self, x_old, x_new, state_old, state_trial, &
                                response, dt, tangent, point_status)
        class(coupled_test_model_t), intent(in) :: self
        real(wp), intent(in) :: x_old(:,:), x_new(:,:), state_old(:,:)
        real(wp), intent(out) :: state_trial(:,:), response(:,:)
        real(wp), intent(in) :: dt
        real(wp), intent(out), optional :: tangent(:,:,:)
        integer, intent(inout) :: point_status(:)

        response(1,:) = self%a * x_new(1,:) + self%b * x_new(2,:)
        response(2,:) = self%b * x_new(1,:) + self%c * x_new(2,:)
        response(3,:) = self%b * (x_new(1,:) - 2.0_wp*x_new(2,:))
        state_trial = state_old
        if (present(tangent)) then
            tangent(1,1,:) = self%a
            tangent(1,2,:) = self%b
            tangent(2,1,:) = self%b
            tangent(2,2,:) = self%c
            tangent(3,1,:) = self%b
            tangent(3,2,:) = -2.0_wp*self%b
        end if
        point_status = cf_ok
    end subroutine evaluate_coupled
    subroutine sensitivity_coupled(self, x_old, x_new, state_old, state_trial, &
                                   selected_parameters, sensitivity_old, &
                                   sensitivity_trial, response, dsdp, dt, tangent, &
                                   point_status)
        class(coupled_test_model_t), intent(in) :: self
        real(wp), intent(in) :: x_old(:,:), x_new(:,:), state_old(:,:)
        integer, intent(in) :: selected_parameters(:)
        real(wp), intent(in) :: sensitivity_old(:,:,:)
        real(wp), intent(out) :: state_trial(:,:), sensitivity_trial(:,:,:)
        real(wp), intent(out) :: response(:,:), dsdp(:,:,:)
        real(wp), intent(in) :: dt
        real(wp), intent(out), optional :: tangent(:,:,:)
        integer, intent(inout) :: point_status(:)

        response(1,:) = self%a*x_new(1,:) + self%b*x_new(2,:)
        response(2,:) = self%b*x_new(1,:) + self%c*x_new(2,:)
        response(3,:) = self%b*(x_new(1,:) - 2.0_wp*x_new(2,:))
        state_trial = state_old
        sensitivity_trial = sensitivity_old
        dsdp(1,1,:) = x_new(2,:)
        dsdp(2,1,:) = x_new(1,:)
        dsdp(3,1,:) = x_new(1,:) - 2.0_wp*x_new(2,:)
        if (present(tangent)) then
            tangent(1,1,:) = self%a
            tangent(1,2,:) = self%b
            tangent(2,1,:) = self%b
            tangent(2,2,:) = self%c
            tangent(3,1,:) = self%b
            tangent(3,2,:) = -2.0_wp*self%b
        end if
        point_status = cf_ok
    end subroutine sensitivity_coupled
end module constituf_test_coupled

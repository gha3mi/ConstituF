! SPDX-License-Identifier: Apache-2.0

module constituf_model
    !! Extension interface for independently implemented constitutive models.
    !!
    !! The core never imports concrete material implementations. Model code
    !! receives packed numeric views and writes into caller-owned buffers.
    !! `point_status` is prefilled with failure codes before every batch, so
    !! implementations must set `cf_ok` for each successfully updated point.
    !!
    !! Numerical inputs, parameters, time increments and history use a consistent
    !! unit system chosen by the caller. Model authors document dimensional
    !! relationships and evaluate these values directly.
    !!
    !! `x_new` has shape `(input_components, npoints)`. `x_old` has the same
    !! shape only if the model declares `set_requires_previous_input(.true.)`.
    !! Otherwise it has shape `(0,0)` and must not be indexed. `state_old` and
    !! `state_trial` have shape `(state_components, npoints)`. `response` has
    !! shape `(output_components,npoints)`. If present, `tangent` is the exact
    !! Jacobian of the discrete response with respect to `x_new`, at fixed
    !! previous state and input, shaped `(nout,nin,npoints)`.
    !!
    !! Optional analytical sensitivities differentiate the *discrete* path
    !! update at fixed prescribed input history. The author supplies a fused
    !! batch method. No finite-difference fallback is permitted.
    !! `selected_parameters` contains one-based IDs resolved from the model's
    !! scalar parameter declarations during setup. `sensitivity_old` and
    !! `sensitivity_trial` have shape `(state,nselected,npoints)`. `dsdp` has
    !! shape `(output,nselected,npoints)`. A model must propagate history
    !! derivatives and handle branch nondifferentiability explicitly.
    !! The discrete constitutive map is
    !! \[s_{n+1}=G(x_{n+1},x_n,s_n,p,\Delta t),\quad
    !! y_{n+1}=H(x_{n+1},x_n,s_n,p,\Delta t).\]
    !! For prescribed parameter-independent inputs and time increments:
    !! \[\frac{\mathrm{d}s_{n+1}}{\mathrm{d}p}
    !! =G_p+G_{s_n}\frac{\mathrm{d}s_n}{\mathrm{d}p},\qquad
    !! \frac{\mathrm{d}y_{n+1}}{\mathrm{d}p}
    !! =H_p+H_{s_n}\frac{\mathrm{d}s_n}{\mathrm{d}p}.\]
    !! Partial derivatives above hold the previous converged state fixed.
    !! The total derivatives include all accepted material history steps.
    !! The algorithmic tangent is \(C_{\mathrm{alg}}=H_{x_{n+1}}\) at
    !! fixed \(s_n,x_n,p,\Delta t\). The instantaneous parameter partial
    !! is \(H_p\) at fixed \(s_n,x_{n+1},x_n,\Delta t\). Neither is the
    !! total derivative \(\mathrm{d}y_{n+1}/\mathrm{d}p\) returned in `dsdp`.
    !! The tangent shape is `(noutput,ninput,npoints)`, including non-square
    !! block Jacobians for coupled physical fields.
    !!
    !! When external driving histories depend on the design parameters, the
    !! *general* chain rule additionally includes
    !! \[G_{x_{n+1}}\frac{\mathrm{d}x_{n+1}}{\mathrm{d}p}
    !! +G_{x_n}\frac{\mathrm{d}x_n}{\mathrm{d}p}
    !! +G_{\Delta t}\frac{\mathrm{d}\Delta t}{\mathrm{d}p}\]
    !! in the state update, with corresponding \(H\)-derivative terms in the
    !! response. The current interface intentionally does not accept those
    !! driving-history derivatives. `dsdp` is defined only for prescribed,
    !! parameter-independent inputs and time increments. It supports selected
    !! scalar parameters. Tensor parameters, tangent-parameter derivatives,
    !! adjoint and directional products are not exposed by this interface.
    !! For an implicit local residual \(R(z,s_n,x_{n+1},p)=0\), where
    !! \(R_z\) is nonsingular, the consistent local derivative follows from
    !! \[\frac{\mathrm{d}z}{\mathrm{d}p}=-R_z^{-1}
    !! \left(R_p+R_{s_n}\frac{\mathrm{d}s_n}{\mathrm{d}p}\right).\]
    !! A material author may reuse an existing factorization of \(R_z\),
    !! never replace this derivative with parameter perturbations.
    !! Integer/logical history updates use optional typed batch methods.
    !! The real-only path remains unchanged and does not allocate discrete
    !! history buffers. For typed models, the author overrides
    !! `eval_typed_batch` and, when analytical sensitivities are advertised,
    !! `eval_typed_sensitivity_batch`. Typed initial states default to zero
    !! and `.false.` or may be initialized in `initialize_typed_state`.
    !! All typed history arrays have shape `(component, point)` and the
    !! integer kind is `ik` (int32). An absent typed kind has zero components.
    !! Discrete states are not differentiated as continuous variables.
    !! Branch switching must be reported explicitly where a classical
    !! sensitivity does not exist.
    !! Global FEM solution sensitivities additionally require displacement
    !! and other solver-field derivatives, which are not computed here.
    !! Reference: J. C. Simo and R. L. Taylor, *Consistent Tangent Operators
    !! for Rate Independent Elastoplasticity*, CMAME **48** (1985), 101-118.
    use constituf_kinds, only: wp, ik
    use constituf_spec, only: material_spec_t
    use constituf_status, only: cf_ok, cf_unsupported
    implicit none (type, external)
    private

    public :: material_model_t

    type, abstract :: material_model_t
        !! Abstract material definition, separate from integration-point state.
        !!
        !! A new material implements `describe` and `eval_batch` only.
        !! Analytical parameter derivatives are an optional extension.
        !! The optional initial state defaults to zero and can be overridden.
        !! Model parameters must be treated as immutable after initialization.
    contains
        procedure(describe_interface), deferred :: describe
            !! Describe inputs, outputs, history fields and supported capabilities.
        procedure(evaluate_interface), deferred :: eval_batch
            !! Evaluate the constitutive response and trial history over a batch.
        procedure :: initialize_state => zero_initial_state
            !! Construct the initial real-valued material history.
        procedure :: initialize_typed_state => zero_typed_initial_state
            !! Construct initial integer and logical material histories.
        procedure :: eval_typed_batch => unsupported_typed_batch
            !! Evaluate a batch with integer and logical history components.
        procedure :: eval_typed_sensitivity_batch => unsupported_typed_sensitivity
            !! Evaluate exact sensitivities for a batch with typed history.
        procedure :: initialize_sensitivity => default_initial_sensitivity
            !! Initialize analytical derivatives of the material history.
        procedure :: eval_sensitivity_batch => unsupported_sensitivity
            !! Evaluate the response and total analytical parameter sensitivities.
    end type material_model_t

    abstract interface
        subroutine describe_interface(self, spec, status)
            !! Declare fields, parameter IDs and capabilities in `spec`.
            !! Return cf_ok only for a valid parameterized model. Invalid
            !! descriptions are rejected before material-point allocation.
            import :: material_model_t, material_spec_t
            implicit none (type, external)
            class(material_model_t), intent(in) :: self
                !! Instance of `material_model_t` accessed by this operation.
            type(material_spec_t), intent(out) :: spec
                !! Schema populated with fields, parameter metadata and model capabilities.
            integer, intent(out) :: status
                !! Operation status code; `cf_ok` indicates success.
        end subroutine describe_interface

        subroutine evaluate_interface(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
            !! Batched trial update: `x_new(nin,N)`, `state_old(ns,N)`,
            !! `state_trial(ns,N)` and `response(nout,N)`. When supported,
            !! `tangent(nout,nin,N)` is the exact current-input Jacobian.
            !! `x_old` has shape `(nin,N)` only when previous input is required.
            !! Otherwise it has shape (0,0). `dt` is positive when required.
            !! Set `point_status(i)=cf_ok` only for successfully updated points.
            !! Do not mutate committed arrays or allocate inside point loops.
            import :: material_model_t, wp
            implicit none (type, external)
            class(material_model_t), intent(in) :: self
                !! Instance of `material_model_t` accessed by this operation.
            real(wp), intent(in) :: x_old(:,:)
                !! Previous inputs `(ninput,npoints)` when required; otherwise zero-sized.
            real(wp), intent(in) :: x_new(:,:)
                !! Current prescribed inputs `(ninput,npoints)`.
            real(wp), intent(in) :: state_old(:,:)
                !! Committed real material history `(nstate,npoints)`.
            real(wp), intent(out) :: state_trial(:,:)
                !! Updated trial real history `(nstate,npoints)` before acceptance.
            real(wp), intent(out) :: response(:,:)
                !! Constitutive outputs `(noutput,npoints)` for current prescribed inputs.
            real(wp), intent(in) :: dt
                !! Prescribed time increment; strictly positive when required by the model.
            real(wp), intent(out), optional :: tangent(:,:,:)
                !! Optional exact discrete input Jacobian `(noutput,ninput,npoints)`.
            integer, intent(inout) :: point_status(:)
                !! Per-point status codes; set to `cf_ok` for successful updates.
        end subroutine evaluate_interface
    end interface

contains

    subroutine zero_initial_state(self, initial_input, state, status)
        !! Default state initialization for stateless or zero-history models.
        !! `initial_input` can be used by specialized models overriding this
        !! method to initialize a nonzero history from the reference fields.
        class(material_model_t), intent(in) :: self
            !! Instance of `material_model_t` accessed by this operation.
        real(wp), intent(in) :: initial_input(:,:)
            !! Prescribed initial packed inputs, supplied when required by the model.
        real(wp), intent(out) :: state(:,:)
            !! Real-valued material history to initialize.
        integer, intent(out) :: status
            !! Operation status code; `cf_ok` indicates success.

        state = 0.0_wp
        status = cf_ok
    end subroutine zero_initial_state

    subroutine zero_typed_initial_state(self, initial_input, integer_state, logical_state, status)
        !! Default initial integer and logical history is zero/false.
        !! Override for material-specific nonzero initial discrete states.
        class(material_model_t), intent(in) :: self
            !! Instance of `material_model_t` accessed by this operation.
        real(wp), intent(in) :: initial_input(:,:)
            !! Prescribed initial packed inputs, supplied when required by the model.
        integer(ik), intent(out) :: integer_state(:,:)
            !! Packed accepted or restored integer history `(ninteger,npoints)`.
        logical, intent(out) :: logical_state(:,:)
            !! Packed accepted or restored logical history `(nlogical,npoints)`.
        integer, intent(out) :: status
            !! Operation status code; `cf_ok` indicates success.

        integer_state = 0_ik
        logical_state = .false.
        status = cf_ok
    end subroutine zero_typed_initial_state

    subroutine unsupported_typed_batch(self, x_old, x_new, state_old, state_trial, integer_old, integer_trial, logical_old, logical_trial, response, dt, tangent, point_status)
        !! A model declaring discrete state must explicitly implement its update.
        class(material_model_t), intent(in) :: self
            !! Instance of `material_model_t` accessed by this operation.
        real(wp), intent(in) :: x_old(:,:)
            !! Previous inputs `(ninput,npoints)` when required; otherwise zero-sized.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed inputs `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real material history `(nstate,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Updated trial real history `(nstate,npoints)` before acceptance.
        real(wp), intent(out) :: response(:,:)
            !! Constitutive outputs `(noutput,npoints)` for current prescribed inputs.
        integer(ik), intent(in) :: integer_old(:,:)
            !! Committed integer history `(ninteger,npoints)`.
        integer(ik), intent(out) :: integer_trial(:,:)
            !! Updated trial integer history `(ninteger,npoints)`.
        logical, intent(in) :: logical_old(:,:)
            !! Committed logical history `(nlogical,npoints)`.
        logical, intent(out) :: logical_trial(:,:)
            !! Updated trial logical history `(nlogical,npoints)`.
        real(wp), intent(in) :: dt
            !! Prescribed time increment; strictly positive when required by the model.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact discrete input Jacobian `(noutput,ninput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Per-point status codes; set to `cf_ok` for successful updates.

        state_trial = 0.0_wp
        integer_trial = 0_ik
        logical_trial = .false.
        response = 0.0_wp
        if (present(tangent)) tangent = 0.0_wp
        point_status = cf_unsupported
    end subroutine unsupported_typed_batch

    subroutine unsupported_typed_sensitivity(self, x_old, x_new, state_old, state_trial, integer_old, integer_trial, logical_old, logical_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Discrete active sets have no assumed derivatives. Authors must
        !! implement exact path derivatives on differentiable branches.

        class(material_model_t), intent(in) :: self
            !! Instance of `material_model_t` accessed by this operation.
        real(wp), intent(in) :: x_old(:,:)
            !! Previous inputs `(ninput,npoints)` when required; otherwise zero-sized.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed inputs `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real material history `(nstate,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Updated trial real history `(nstate,npoints)` before acceptance.
        real(wp), intent(out) :: response(:,:)
            !! Constitutive outputs `(noutput,npoints)` for current prescribed inputs.
        integer(ik), intent(in) :: integer_old(:,:)
            !! Committed integer history `(ninteger,npoints)`.
        integer(ik), intent(out) :: integer_trial(:,:)
            !! Updated trial integer history `(ninteger,npoints)`.
        logical, intent(in) :: logical_old(:,:)
            !! Committed logical history `(nlogical,npoints)`.
        logical, intent(out) :: logical_trial(:,:)
            !! Updated trial logical history `(nlogical,npoints)`.
        integer, intent(in) :: selected_parameters(:)
            !! One-based model parameter IDs in the requested derivative-axis order.
        real(wp), intent(in) :: sensitivity_old(:,:,:)
            !! Committed real-state derivatives `(nstate,nselected,npoints)`.
        real(wp), intent(out) :: sensitivity_trial(:,:,:)
            !! Updated trial real-state derivatives `(nstate,nselected,npoints)`.
        real(wp), intent(out) :: dsdp(:,:,:)
            !! Total analytical response derivatives `(noutput,nselected,npoints)` at fixed prescribed driving history.
        real(wp), intent(in) :: dt
            !! Prescribed time increment; strictly positive when required by the model.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact discrete input Jacobian `(noutput,ninput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Per-point status codes; set to `cf_ok` for successful updates.

        state_trial = 0.0_wp
        integer_trial = 0_ik
        logical_trial = .false.
        sensitivity_trial = 0.0_wp
        response = 0.0_wp
        dsdp = 0.0_wp
        if (present(tangent)) tangent = 0.0_wp
        point_status = cf_unsupported
    end subroutine unsupported_typed_sensitivity

    subroutine default_initial_sensitivity(self, initial_input, initial_state, selected_parameters, sensitivity, status)
        !! A stateless material has identically zero initial state derivatives.
        !! State-bearing materials must override this, even if initial state
        !! sensitivities are zero. Silently assuming zero could be incorrect.
        class(material_model_t), intent(in) :: self
            !! Instance of `material_model_t` accessed by this operation.
        real(wp), intent(in) :: initial_input(:,:)
            !! Prescribed initial packed inputs, supplied when required by the model.
        real(wp), intent(in) :: initial_state(:,:)
            !! Initialized real-valued history `(nstate,npoints)`.
        integer, intent(in) :: selected_parameters(:)
            !! One-based model parameter IDs in the requested derivative-axis order.
        real(wp), intent(out) :: sensitivity(:,:,:)
            !! Initial exact derivatives of real history `(nstate,nselected,npoints)`.
        integer, intent(out) :: status
            !! Operation status code; `cf_ok` indicates success.

        sensitivity = 0.0_wp
        status = cf_ok
        if (size(initial_state, 1) > 0) status = cf_unsupported
    end subroutine default_initial_sensitivity

    subroutine unsupported_sensitivity(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Unimplemented sensitivities are never approximated silently.
        class(material_model_t), intent(in) :: self
            !! Instance of `material_model_t` accessed by this operation.
        real(wp), intent(in) :: x_old(:,:)
            !! Previous inputs `(ninput,npoints)` when required; otherwise zero-sized.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed inputs `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real material history `(nstate,npoints)`.
        integer, intent(in) :: selected_parameters(:)
            !! One-based model parameter IDs in the requested derivative-axis order.
        real(wp), intent(in) :: sensitivity_old(:,:,:)
            !! Committed real-state derivatives `(nstate,nselected,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Updated trial real history `(nstate,npoints)` before acceptance.
        real(wp), intent(out) :: sensitivity_trial(:,:,:)
            !! Updated trial real-state derivatives `(nstate,nselected,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Constitutive outputs `(noutput,npoints)` for current prescribed inputs.
        real(wp), intent(out) :: dsdp(:,:,:)
            !! Total analytical response derivatives `(noutput,nselected,npoints)` at fixed prescribed driving history.
        real(wp), intent(in) :: dt
            !! Prescribed time increment; strictly positive when required by the model.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact discrete input Jacobian `(noutput,ninput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Per-point status codes; set to `cf_ok` for successful updates.

        state_trial = 0.0_wp
        sensitivity_trial = 0.0_wp
        response = 0.0_wp
        dsdp = 0.0_wp
        if (present(tangent)) tangent = 0.0_wp
        point_status = cf_unsupported
    end subroutine unsupported_sensitivity
end module constituf_model

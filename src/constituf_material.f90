! SPDX-License-Identifier: Apache-2.0

module constituf_material
    !! User-facing high-performance batch material handle.
    !!
    !! Packed arrays have layout `(component, point)`. Constitutive kernels
    !! are invoked once per batch, with no per-point dynamic dispatch, memory
    !! allocation, or field-name resolution by the framework.
    !!
    !! The caller owns unit consistency. Parameters, prescribed fields, `dt`
    !! and imported history must use the same chosen system throughout the
    !! handle's lifetime. Metadata stores no unit labels. Each tangent entry
    !! has output/input units and each parameter derivative has output/parameter
    !! units, including different units across coupled-field blocks.
    !!
    !! Sensitivities are strictly opt-in at `init(sensitivity_parameters=...)`.
    !! Selected parameter names are resolved once. Every enabled `eval` must
    !! supply `dsdp(noutput,nselected,npoints)` and computes both response and
    !! total path derivatives with respect to parameters at fixed prescribed
    !! input history. `commit/revert` also protect derivative histories.
    !! An unsupported analytical extension returns an error. Never use a
    !! finite-difference fallback.
    !!
    !! For \(N\) points, \(n_s\) real history components per point and
    !! \(m\) selected scalar parameters, retained double-buffered state uses
    !! \[B_{\mathrm{real}}=2Nn_s b_{\mathrm{wp}}\]
    !! bytes and path derivatives add \(2Nn_sm b_{\mathrm{wp}}\) bytes,
    !! where \(b_{\mathrm{wp}}=8\) bytes at `wp=real64`.
    !! Integer and logical histories, when declared, also have committed
    !! and trial buffers. Required previous inputs use two more buffers
    !! of `(ninput,npoints)` but are excluded from `state_storage_bytes()`,
    !! as are metadata, point statuses and caller-owned response arrays.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_kinds, only: wp, ik, lk
    use constituf_model, only: material_model_t
    use constituf_spec, only: material_spec_t
    use constituf_status, only: cf_ok, cf_invalid_argument, cf_invalid_shape, cf_uninitialized, cf_no_trial, cf_invalid_time, cf_unsupported, cf_nondifferentiable, cf_evaluation_failed, cf_allocation_failed
    implicit none (type, external)
    private

    public :: material_t

    type :: material_t
        !! Material model with independently owned, preallocated point states.
        !!
        !! History of each declared kind is owned by the handle. Integer
        !! and logical histories use a distinct optional typed batch update.
        !! A single handle must not be evaluated concurrently from multiple
        !! threads: the trial buffers are mutable. The batch kernel may operate
        !! on disjoint points in parallel. The caller owns response arrays.
        private
        class(material_model_t), allocatable :: model
            !! Owned immutable constitutive model instance.
        type(material_spec_t) :: spec
            !! Validated packed field and parameter schema.
        real(wp), allocatable :: x_committed(:,:)
            !! Packed inputs at the last committed step, when retained.
        real(wp), allocatable :: x_trial(:,:)
            !! Packed inputs associated with the current trial.
        real(wp), allocatable :: state_committed(:,:)
            !! Accepted real-valued internal state for all points.
        real(wp), allocatable :: state_trial(:,:)
            !! Uncommitted real-valued internal state for all points.
        integer(ik), allocatable :: integer_committed(:,:)
            !! Accepted integer-valued material history for all points.
        integer(ik), allocatable :: integer_trial(:,:)
            !! Uncommitted integer-valued material history for all points.
        logical, allocatable :: logical_committed(:,:)
            !! Accepted logical-valued material history for all points.
        logical, allocatable :: logical_trial(:,:)
            !! Uncommitted logical-valued material history for all points.
        integer, allocatable :: point_status(:)
            !! Status code for each point in the last trial evaluation.
        integer, allocatable :: selected_parameters(:)
            !! One-based parameter IDs in sensitivity-axis order.
        real(wp), allocatable :: sensitivity_committed(:,:,:)
            !! Accepted derivatives of real history with respect to selected parameters.
        real(wp), allocatable :: sensitivity_trial(:,:,:)
            !! Trial derivatives of real history with respect to selected parameters.
        integer :: npoints = 0
            !! Number of material points owned by this handle.
        logical :: trial_ready = .false.
            !! True when a successful trial is available for commit.
        logical :: point_failure = .false.
            !! True when the last trial reported at least one point failure.
        logical :: point_results_ready = .false.
            !! True when per-point results are available to query.
    contains
        procedure :: init => initialize_material
            !! Clone the material model and initialize all committed and trial point buffers.
        procedure :: eval => evaluate_material
            !! Evaluate the trial response, optional tangent and selected parameter derivatives.
        procedure :: commit => commit_material
            !! Accept a successful trial as the new committed history.
        procedure :: revert => revert_material
            !! Discard the current trial and retain committed history.
        procedure :: input_size => material_input_size
            !! Return the number of packed input components per material point.
        procedure :: output_size => material_output_size
            !! Return the number of packed response components per point.
        procedure :: state_size => material_state_size
            !! Return the number of real-valued history components per point.
        procedure :: integer_state_size => material_integer_state_size
            !! Return the number of integer history components per point.
        procedure :: logical_state_size => material_logical_state_size
            !! Return the number of logical history components per point.
        procedure :: state_storage_bytes => material_state_storage_bytes
            !! Return the memory footprint of retained history buffers in bytes.
        procedure :: sensitivity_size => material_sensitivity_size
            !! Return the number of selected scalar sensitivity parameters.
        procedure :: parameter_count => material_parameter_count
            !! Return the number of declared scalar material parameters.
        procedure :: parameter_info => material_parameter_info
            !! Get the name of a parameter by one-based index.
        procedure :: point_count => material_point_count
            !! Return the number of points managed by this handle.
        procedure :: previous_input_size => material_previous_input_size
            !! Return the number of previously retained input components per point.
        procedure :: sensitivity_parameter_id => material_sensitivity_parameter_id
            !! Return the model parameter ID at a selected sensitivity index.
        procedure :: input_offset => material_input_offset
            !! Find a named input field's first packed index.
        procedure :: output_offset => material_output_offset
            !! Find a named output field's first packed index.
        procedure :: input_components => material_input_components
            !! Return a named input field's packed component count.
        procedure :: output_components => material_output_components
            !! Return a named output field's packed component count.
        procedure :: field_matches => material_field_matches
            !! Test a field's role and tensor metadata.
        procedure :: field_info => material_field_info
            !! Retrieve metadata for a named field of the specified role.
        procedure :: last_failure_point
            !! Return the index of the first failed point, or zero if none.
        procedure :: point_statuses => material_point_statuses
            !! Copy the last per-point evaluation status codes.
        procedure :: export_committed => material_export_committed
            !! Copy committed real, typed and sensitivity histories into caller buffers.
        procedure :: import_committed => material_import_committed
            !! Replace committed histories after validating a complete state snapshot.
    end type material_t

contains

    subroutine initialize_material(self, model, npoints, status, initial_input, sensitivity_parameters)
        !! Initialize a material and its state without modifying an existing
        !! handle on failure. Setup is permitted to allocate and inspect names.
        !!
        !! `npoints` must be positive. `initial_input`, when supplied, has
        !! shape (ninput,npoints) and must be finite. If omitted, the reference
        !! input is zero. Models unable to use zero must require this argument.
        !! `sensitivity_parameters` is optional. If present it names distinct
        !! scalar parameters in desired derivative-axis order. A state-bearing
        !! model must analytically initialize its parameter sensitivity state.
        !! Names, validation and all allocations occur only in this method.
        !! State-bearing models must provide their exact initial state
        !! sensitivities. By default only stateless models can initialize
        !! derivative tracking without an override.
        !! `status` is cf_ok on success. Invalid names, parameter domains and
        !! nonfinite initial values fail setup without replacing the handle.
        class(material_t), intent(inout) :: self
            !! Instance of `material_t` accessed by this operation.
        class(material_model_t), intent(in) :: model
            !! Concrete constitutive model cloned into the material handle.
        integer, intent(in) :: npoints
            !! Positive number of material points in the batch.
        integer, intent(out) :: status
            !! Operation status code; `cf_ok` indicates success.
        real(wp), intent(in), optional :: initial_input(:,:)
            !! Prescribed initial packed inputs, supplied when required by the model.
        character(*), intent(in), optional :: sensitivity_parameters(:)
            !! Optional distinct scalar parameter names in derivative-axis order.
        type(material_t) :: next
        real(wp), allocatable :: setup_input(:,:)
        integer :: allocation_status
        integer :: i
        integer :: parameter_id
        character(256) :: allocation_message

        status = cf_invalid_argument
        if (npoints < 1) return
        allocate(next%model, source=model, stat=allocation_status, errmsg=allocation_message)
        if (allocation_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        call next%model%describe(next%spec, status)
        if (status /= cf_ok) return
        if (.not. next%spec%valid()) then
            status = cf_invalid_argument
            return
        end if
        if (next%spec%requires_initial_input() .and. .not. present(initial_input)) then
            status = cf_invalid_argument
            return
        end if
        if (present(initial_input)) then
            if (size(initial_input, 1) /= next%spec%input_size() .or. size(initial_input, 2) /= npoints) then
                status = cf_invalid_shape
                return
            end if
            if (.not. all(ieee_is_finite(initial_input))) then
                status = cf_invalid_argument
                return
            end if
        end if

        allocate(setup_input(next%spec%input_size(), npoints), stat=allocation_status, errmsg=allocation_message)
        if (allocation_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        allocate(next%state_committed(next%spec%state_size(), npoints), stat=allocation_status, errmsg=allocation_message)
        if (allocation_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        allocate(next%state_trial(next%spec%state_size(), npoints), stat=allocation_status, errmsg=allocation_message)
        if (allocation_status /= 0) then
            status = cf_allocation_failed
            return
        end if
        if (next%spec%integer_state_size() > 0 .or. next%spec%logical_state_size() > 0) then
            allocate(next%integer_committed(next%spec%integer_state_size(), npoints), stat=allocation_status, errmsg=allocation_message)
            if (allocation_status /= 0) then
                status = cf_allocation_failed
                return
            end if
            allocate(next%integer_trial(next%spec%integer_state_size(), npoints), stat=allocation_status, errmsg=allocation_message)
            if (allocation_status /= 0) then
                status = cf_allocation_failed
                return
            end if
            allocate(next%logical_committed(next%spec%logical_state_size(), npoints), stat=allocation_status, errmsg=allocation_message)
            if (allocation_status /= 0) then
                status = cf_allocation_failed
                return
            end if
            allocate(next%logical_trial(next%spec%logical_state_size(), npoints), stat=allocation_status, errmsg=allocation_message)
            if (allocation_status /= 0) then
                status = cf_allocation_failed
                return
            end if
        end if
        allocate(next%point_status(npoints), stat=allocation_status, errmsg=allocation_message)
        if (allocation_status /= 0) then
            status = cf_allocation_failed
            return
        end if

        setup_input = 0.0_wp
        if (present(initial_input)) setup_input = initial_input
        call next%model%initialize_state(setup_input, next%state_committed, status)
        if (status /= cf_ok) return
        if (.not. all(ieee_is_finite(next%state_committed))) then
            status = cf_invalid_argument
            return
        end if

        if (allocated(next%integer_committed)) then
            call next%model%initialize_typed_state(setup_input, next%integer_committed, next%logical_committed, status)
            if (status /= cf_ok) return
        end if

        if (present(sensitivity_parameters)) then
            status = cf_invalid_argument
            if (size(sensitivity_parameters) == 0) return
            if (.not. next%spec%supports_sensitivity()) then
                status = cf_unsupported
                return
            end if
            allocate(next%selected_parameters(size(sensitivity_parameters)), stat=allocation_status, errmsg=allocation_message)
            if (allocation_status /= 0) then
                status = cf_allocation_failed
                return
            end if
            do i = 1, size(sensitivity_parameters)
                parameter_id = next%spec%parameter_id(sensitivity_parameters(i))
                if (parameter_id == 0) return
                if (any(next%selected_parameters(:i-1) == parameter_id)) return
                next%selected_parameters(i) = parameter_id
            end do
            allocate(next%sensitivity_committed(next%spec%state_size(), size(sensitivity_parameters), npoints), stat=allocation_status, errmsg=allocation_message)
            if (allocation_status /= 0) then
                status = cf_allocation_failed
                return
            end if
            allocate(next%sensitivity_trial(next%spec%state_size(), size(sensitivity_parameters), npoints), stat=allocation_status, errmsg=allocation_message)
            if (allocation_status /= 0) then
                status = cf_allocation_failed
                return
            end if
            call next%model%initialize_sensitivity(setup_input, next%state_committed, next%selected_parameters, next%sensitivity_committed, status)
            if (status /= cf_ok) return
            if (.not. all(ieee_is_finite(next%sensitivity_committed))) then
                status = cf_invalid_argument
                return
            end if
        end if

        if (next%spec%requires_previous_input()) then
            call move_alloc(setup_input, next%x_committed)
            allocate(next%x_trial(next%spec%input_size(), npoints), stat=allocation_status, errmsg=allocation_message)
            if (allocation_status /= 0) then
                status = cf_allocation_failed
                return
            end if
        else
            allocate(next%x_committed(0,0), stat=allocation_status, errmsg=allocation_message)
            if (allocation_status /= 0) then
                status = cf_allocation_failed
                return
            end if
            allocate(next%x_trial(0,0), stat=allocation_status, errmsg=allocation_message)
            if (allocation_status /= 0) then
                status = cf_allocation_failed
                return
            end if
        end if

        self%spec = next%spec
        call move_alloc(next%model, self%model)
        call move_alloc(next%x_committed, self%x_committed)
        call move_alloc(next%x_trial, self%x_trial)
        call move_alloc(next%state_committed, self%state_committed)
        call move_alloc(next%state_trial, self%state_trial)
        if (allocated(next%integer_committed)) then
            call move_alloc(next%integer_committed, self%integer_committed)
            call move_alloc(next%integer_trial, self%integer_trial)
            call move_alloc(next%logical_committed, self%logical_committed)
            call move_alloc(next%logical_trial, self%logical_trial)
        else
            if (allocated(self%integer_committed)) deallocate(self%integer_committed)
            if (allocated(self%integer_trial)) deallocate(self%integer_trial)
            if (allocated(self%logical_committed)) deallocate(self%logical_committed)
            if (allocated(self%logical_trial)) deallocate(self%logical_trial)
        end if
        call move_alloc(next%point_status, self%point_status)
        if (allocated(next%selected_parameters)) then
            call move_alloc(next%selected_parameters, self%selected_parameters)
            call move_alloc(next%sensitivity_committed, self%sensitivity_committed)
            call move_alloc(next%sensitivity_trial, self%sensitivity_trial)
        else
            if (allocated(self%selected_parameters)) deallocate(self%selected_parameters)
            if (allocated(self%sensitivity_committed)) deallocate(self%sensitivity_committed)
            if (allocated(self%sensitivity_trial)) deallocate(self%sensitivity_trial)
        end if
        self%npoints = npoints
        self%trial_ready = .false.
        self%point_failure = .false.
        self%point_status(:) = cf_ok
        status = cf_ok
    end subroutine initialize_material

    subroutine evaluate_material(self, x, response, tangent, status, dt, dsdp)
        !! Compute trial response from committed history. Never commit here.
        !!
        !! `x(input_component,point)` and `response(output_component,point)`
        !! use the field order declared by the material. The optional tangent
        !! has shape `(output_component,input_component,point)` and is the
        !! algorithmic derivative of the discrete update at fixed committed
        !! state and prescribed previous input. `dt` is mandatory for models
        !! whose schema requests a positive time increment.
        !!
        !! When sensitivities were enabled in `init`, `dsdp` is mandatory
        !! and has shape `(noutput,nselected,npoints)`. It is the exact total
        !! derivative of the *discrete* material response with respect to
        !! selected parameters at fixed prescribed input histories. This is
        !! not the derivative of a global equilibrium solution when the
        !! external driving fields themselves depend on parameters. Every
        !! trial starts from committed state and committed state derivatives.
        !!
        !! After an error, output buffers may contain partial values and must
        !! not be consumed. An unsuccessful trial cannot be committed.
        !! Only models declaring typed history take the typed batch path.
        !! Discrete history is never differentiated as real-valued state.
        !! Inputs and outputs must not alias each other. Kernels must report
        !! invalid physical inputs through point statuses. This fast path does
        !! not scan every input component for finiteness.
        class(material_t), intent(inout) :: self
            !! Instance of `material_t` accessed by this operation.
        real(wp), intent(in) :: x(:,:)
            !! Packed current inputs `(ninput,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Constitutive outputs `(noutput,npoints)` for current prescribed inputs.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact discrete input Jacobian `(noutput,ninput,npoints)`.
        integer, intent(out) :: status
            !! Operation status code; `cf_ok` indicates success.
        real(wp), intent(in), optional :: dt
            !! Prescribed time increment; strictly positive when required by the model.
        real(wp), intent(out), optional :: dsdp(:,:,:)
            !! Total analytical response derivatives `(noutput,nselected,npoints)` at fixed prescribed driving history.
        real(wp) :: time_step

        self%trial_ready = .false.
        self%point_failure = .false.
        self%point_results_ready = .false.
        status = cf_uninitialized
        if (.not. allocated(self%model)) return

        status = cf_invalid_shape
        if (size(x, 1) /= self%spec%input_size() .or. size(x, 2) /= self%npoints) return
        if (size(response, 1) /= self%spec%output_size() .or. size(response, 2) /= self%npoints) return
        if (present(tangent)) then
            if (size(tangent, 1) /= self%spec%output_size() .or. size(tangent, 2) /= self%spec%input_size() .or. size(tangent, 3) /= self%npoints) return
            if (.not. self%spec%supports_tangent()) then
                status = cf_unsupported
                return
            end if
        end if

        if (allocated(self%selected_parameters)) then
            status = cf_invalid_argument
            if (.not. present(dsdp)) return
            if (size(dsdp, 1) /= self%spec%output_size() .or. size(dsdp, 2) /= size(self%selected_parameters) .or. size(dsdp, 3) /= self%npoints) then
                status = cf_invalid_shape
                return
            end if
        else
            if (present(dsdp)) then
                status = cf_unsupported
                return
            end if
        end if

        time_step = 0.0_wp
        if (present(dt)) time_step = dt
        if (self%spec%requires_dt()) then
            status = cf_invalid_time
            if (.not. ieee_is_finite(time_step)) return
            if (time_step <= 0.0_wp) return
        end if

        self%point_status(:) = cf_evaluation_failed
        self%point_results_ready = .true.
        if (allocated(self%integer_committed)) then
            if (allocated(self%selected_parameters)) then
                call self%model%eval_typed_sensitivity_batch(self%x_committed, x, self%state_committed, self%state_trial, self%integer_committed, self%integer_trial, self%logical_committed, self%logical_trial, self%selected_parameters, self%sensitivity_committed, self%sensitivity_trial, response, dsdp, time_step, tangent, self%point_status)
            else
                call self%model%eval_typed_batch(self%x_committed, x, self%state_committed, self%state_trial, self%integer_committed, self%integer_trial, self%logical_committed, self%logical_trial, response, time_step, tangent, self%point_status)
            end if
        else if (allocated(self%selected_parameters)) then
            call self%model%eval_sensitivity_batch(self%x_committed, x, self%state_committed, self%state_trial, self%selected_parameters, self%sensitivity_committed, self%sensitivity_trial, response, dsdp, time_step, tangent, self%point_status)
        else
            call self%model%eval_batch(self%x_committed, x, self%state_committed, self%state_trial, response, time_step, tangent, self%point_status)
        end if
        status = cf_evaluation_failed
        if (any(self%point_status /= cf_ok)) then
            self%point_failure = .true.
            if (all(self%point_status == cf_unsupported)) status = cf_unsupported
            if (any(self%point_status == cf_nondifferentiable)) status = cf_nondifferentiable
            return
        end if

        if (self%spec%requires_previous_input()) self%x_trial(:,:) = x
        self%trial_ready = .true.
        status = cf_ok
    end subroutine evaluate_material

    subroutine material_export_committed(self, real_state, status, previous_input, integer_state, logical_state, sensitivity_state)
        !! Copy committed integration-point history to caller-owned arrays.
        !!
        !! real_state has shape (state_size,npoints). Optional buffers are
        !! required when the corresponding history is nonempty:
        !! previous_input(ninput,npoints) for incremental laws,
        !! integer_state(ninteger,npoints), logical_state(nlogical,npoints)
        !! and sensitivity_state(nstate,nselected,npoints) when sensitivities
        !! were selected at initialization. Zero-sized kinds may be omitted.
        !! The snapshot never contains an unaccepted trial. This API does not
        !! perform file I/O or encode a portable on-disk representation.
        !! Reinitialize the same model, parameters, point layout and selected
        !! sensitivity parameters before importing into another handle.
        !! Return cf_invalid_shape for missing or mismatched buffers and
        !! cf_uninitialized when the handle has no material.
        class(material_t), intent(in) :: self
            !! Instance of `material_t` accessed by this operation.
        real(wp), intent(out) :: real_state(:,:)
            !! Accepted real history `(nstate,npoints)`.
        integer, intent(out) :: status
            !! Operation status code; `cf_ok` indicates success.
        real(wp), intent(out), optional :: previous_input(:,:)
            !! Packed input at the last accepted step, when retained by the model.
        integer(ik), intent(out), optional :: integer_state(:,:)
            !! Packed accepted or restored integer history `(ninteger,npoints)`.
        logical, intent(out), optional :: logical_state(:,:)
            !! Packed accepted or restored logical history `(nlogical,npoints)`.
        real(wp), intent(out), optional :: sensitivity_state(:,:,:)
            !! Accepted real-history parameter derivatives `(nstate,nselected,npoints)`.

        status = cf_uninitialized
        if (.not. allocated(self%model)) return
        status = cf_invalid_shape
        if (size(real_state,1) /= size(self%state_committed,1) .or. size(real_state,2) /= self%npoints) return
        if (size(self%x_committed,1) /= 0 .and. .not. present(previous_input)) return
        if (present(previous_input)) then
            if (size(previous_input,1) /= size(self%x_committed,1) .or. size(previous_input,2) /= size(self%x_committed,2)) return
        end if
        if (self%spec%integer_state_size() /= 0 .and. .not. present(integer_state)) return
        if (present(integer_state)) then
            if (size(integer_state,1) /= self%spec%integer_state_size() .or. size(integer_state,2) /= self%npoints) return
        end if
        if (self%spec%logical_state_size() /= 0 .and. .not. present(logical_state)) return
        if (present(logical_state)) then
            if (size(logical_state,1) /= self%spec%logical_state_size() .or. size(logical_state,2) /= self%npoints) return
        end if
        if (allocated(self%selected_parameters) .and. .not. present(sensitivity_state)) return
        if (present(sensitivity_state)) then
            if (size(sensitivity_state,1) /= self%spec%state_size() .or. size(sensitivity_state,2) /= self%sensitivity_size() .or. size(sensitivity_state,3) /= self%npoints) return
        end if

        real_state = self%state_committed
        if (present(previous_input)) previous_input = self%x_committed
        if (present(integer_state)) then
            if (allocated(self%integer_committed)) integer_state = self%integer_committed
        end if
        if (present(logical_state)) then
            if (allocated(self%logical_committed)) logical_state = self%logical_committed
        end if
        if (present(sensitivity_state)) then
            if (allocated(self%sensitivity_committed)) sensitivity_state = self%sensitivity_committed
        end if
        status = cf_ok
    end subroutine material_export_committed

    subroutine material_import_committed(self, real_state, status, previous_input, integer_state, logical_state, sensitivity_state)
        !! Restore an accepted material-point snapshot atomically.
        !!
        !! The shape and finite-value checks occur before any committed state
        !! is modified. A successful import invalidates pending trials. A
        !! rejected import leaves both committed and pending trial data alone.
        !! The caller is responsible for externally recording and verifying
        !! the model, parameter values, precision, field layout, point order
        !! and selected sensitivity parameter names. Equal array shapes alone
        !! cannot certify checkpoint compatibility.
        class(material_t), intent(inout) :: self
            !! Instance of `material_t` accessed by this operation.
        real(wp), intent(in) :: real_state(:,:)
            !! Accepted real history `(nstate,npoints)`.
        integer, intent(out) :: status
            !! Operation status code; `cf_ok` indicates success.
        real(wp), intent(in), optional :: previous_input(:,:)
            !! Packed input at the last accepted step, when retained by the model.
        integer(ik), intent(in), optional :: integer_state(:,:)
            !! Packed accepted or restored integer history `(ninteger,npoints)`.
        logical, intent(in), optional :: logical_state(:,:)
            !! Packed accepted or restored logical history `(nlogical,npoints)`.
        real(wp), intent(in), optional :: sensitivity_state(:,:,:)
            !! Accepted real-history parameter derivatives `(nstate,nselected,npoints)`.

        status = cf_uninitialized
        if (.not. allocated(self%model)) return
        status = cf_invalid_shape
        if (size(real_state,1) /= size(self%state_committed,1) .or. size(real_state,2) /= self%npoints) return
        if (size(self%x_committed,1) /= 0 .and. .not. present(previous_input)) return
        if (present(previous_input)) then
            if (size(previous_input,1) /= size(self%x_committed,1) .or. size(previous_input,2) /= size(self%x_committed,2)) return
        end if
        if (self%spec%integer_state_size() /= 0 .and. .not. present(integer_state)) return
        if (present(integer_state)) then
            if (size(integer_state,1) /= self%spec%integer_state_size() .or. size(integer_state,2) /= self%npoints) return
        end if
        if (self%spec%logical_state_size() /= 0 .and. .not. present(logical_state)) return
        if (present(logical_state)) then
            if (size(logical_state,1) /= self%spec%logical_state_size() .or. size(logical_state,2) /= self%npoints) return
        end if
        if (allocated(self%selected_parameters) .and. .not. present(sensitivity_state)) return
        if (present(sensitivity_state)) then
            if (size(sensitivity_state,1) /= self%spec%state_size() .or. size(sensitivity_state,2) /= self%sensitivity_size() .or. size(sensitivity_state,3) /= self%npoints) return
        end if
        status = cf_invalid_argument
        if (.not. all(ieee_is_finite(real_state))) return
        if (present(previous_input)) then
            if (.not. all(ieee_is_finite(previous_input))) return
        end if
        if (present(sensitivity_state)) then
            if (.not. all(ieee_is_finite(sensitivity_state))) return
        end if

        self%state_committed = real_state
        if (present(previous_input)) self%x_committed = previous_input
        if (present(integer_state)) then
            if (allocated(self%integer_committed)) self%integer_committed = integer_state
        end if
        if (present(logical_state)) then
            if (allocated(self%logical_committed)) self%logical_committed = logical_state
        end if
        if (present(sensitivity_state)) then
            if (allocated(self%sensitivity_committed)) self%sensitivity_committed = sensitivity_state
        end if
        self%trial_ready = .false.
        self%point_failure = .false.
        self%point_results_ready = .false.
        self%point_status = cf_ok
        status = cf_ok
    end subroutine material_import_committed

    subroutine commit_material(self, status)
        !! Accept the last successful trial using allocation-descriptor swaps.
        !! Committed state is unchanged after rejected or failed evaluations.
        class(material_t), intent(inout) :: self
            !! Instance of `material_t` accessed by this operation.
        integer, intent(out) :: status
            !! Operation status code; `cf_ok` indicates success.
        real(wp), allocatable :: scratch(:,:)
        real(wp), allocatable :: sensitivity_scratch(:,:,:)
        integer(ik), allocatable :: integer_scratch(:,:)
        logical, allocatable :: logical_scratch(:,:)

        status = cf_no_trial
        if (.not. self%trial_ready) return

        call move_alloc(self%state_committed, scratch)
        call move_alloc(self%state_trial, self%state_committed)
        call move_alloc(scratch, self%state_trial)

        if (allocated(self%integer_committed)) then
            call move_alloc(self%integer_committed, integer_scratch)
            call move_alloc(self%integer_trial, self%integer_committed)
            call move_alloc(integer_scratch, self%integer_trial)
            call move_alloc(self%logical_committed, logical_scratch)
            call move_alloc(self%logical_trial, self%logical_committed)
            call move_alloc(logical_scratch, self%logical_trial)
        end if

        call move_alloc(self%x_committed, scratch)
        call move_alloc(self%x_trial, self%x_committed)
        call move_alloc(scratch, self%x_trial)

        if (allocated(self%selected_parameters)) then
            call move_alloc(self%sensitivity_committed, sensitivity_scratch)
            call move_alloc(self%sensitivity_trial, self%sensitivity_committed)
            call move_alloc(sensitivity_scratch, self%sensitivity_trial)
        end if

        self%trial_ready = .false.
        status = cf_ok
    end subroutine commit_material

    subroutine revert_material(self)
        !! Discard a trial without copying or modifying committed history.
        class(material_t), intent(inout) :: self
            !! Instance of `material_t` accessed by this operation.

        self%trial_ready = .false.
    end subroutine revert_material

    pure integer function material_integer_state_size(self) result(value)
        !! Number of packed integer history components per point.
        class(material_t), intent(in) :: self
            !! Instance of `material_t` accessed by this operation.

        value = self%spec%integer_state_size()
    end function material_integer_state_size

    pure integer function material_logical_state_size(self) result(value)
        !! Number of packed logical history components per point.
        class(material_t), intent(in) :: self
            !! Instance of `material_t` accessed by this operation.

        value = self%spec%logical_state_size()
    end function material_logical_state_size

    pure integer(lk) function material_state_storage_bytes(self) result(value)
        !! Retained committed/trial history and sensitivity storage in bytes.
        !! Excludes metadata, point statuses and caller-owned buffers.
        !! No memory is allocated by this inquiry. It is safe in setup code.
        class(material_t), intent(in) :: self
            !! Instance of `material_t` accessed by this operation.

        value = 0_lk
        if (allocated(self%state_committed)) then
            value = value + 2_lk * size(self%state_committed, kind=lk) * int(storage_size(0.0_wp)/8, lk)
        end if
        if (allocated(self%integer_committed)) then
            value = value + 2_lk * size(self%integer_committed, kind=lk) * int(storage_size(0_ik)/8, lk)
            value = value + 2_lk * size(self%logical_committed, kind=lk) * int(storage_size(.false.)/8, lk)
        end if
        if (allocated(self%sensitivity_committed)) then
            value = value + 2_lk * size(self%sensitivity_committed, kind=lk) * int(storage_size(0.0_wp)/8, lk)
        end if
    end function material_state_storage_bytes

    pure integer function material_point_count(self) result(value)
        !! Count initialized material points. Zero for an empty handle.
        !! Setup-time inquiry for checkpoint and external-solver adapters.
        class(material_t), intent(in) :: self
            !! Instance of `material_t` accessed by this operation.

        value = self%npoints
    end function material_point_count

    pure integer function material_previous_input_size(self) result(value)
        !! Number of committed previous-input components per point.
        !! Returns zero when this model does not retain prior input.
        class(material_t), intent(in) :: self
            !! Instance of `material_t` accessed by this operation.

        value = 0
        if (self%spec%requires_previous_input()) value = self%spec%input_size()
    end function material_previous_input_size

    pure integer function material_sensitivity_parameter_id(self, index) result(value)
        !! Declared one-based scalar parameter ID at selected derivative index.
        !! Returns zero for an unselected or invalid position. The parameter
        !! ordering is immutable after material initialization. This is an
        !! inquiry only. No allocation or name resolution occurs.
        class(material_t), intent(in) :: self
            !! Instance of `material_t` accessed by this operation.
        integer, intent(in) :: index
            !! One-based index identifying the selected entry.

        value = 0
        if (.not. allocated(self%selected_parameters)) return
        if (index < 1 .or. index > size(self%selected_parameters)) return
        value = self%selected_parameters(index)
    end function material_sensitivity_parameter_id

    pure integer function material_sensitivity_size(self) result(value)
        !! Number of scalar material parameters tracked analytically.
        class(material_t), intent(in) :: self
            !! Instance of `material_t` accessed by this operation.

        value = 0
        if (allocated(self%selected_parameters)) value = size(self%selected_parameters)
    end function material_sensitivity_size

    pure integer function material_parameter_count(self) result(value)
        !! Number of declared scalar model parameters, tracked or not.
        !! Available after initialization. Zero for an empty handle.
        class(material_t), intent(in) :: self
            !! Instance of `material_t` accessed by this operation.

        value=self%spec%parameter_count()
    end function material_parameter_count

    subroutine material_parameter_info(self,index,name,found)
        !! One-based parameter name in author order.
        !! No allocation or lookup takes place in constitutive eval().
        class(material_t), intent(in) :: self
            !! Instance of `material_t` accessed by this operation.
        integer, intent(in) :: index
            !! One-based index identifying the selected entry.
        character(:), allocatable, intent(out) :: name
            !! Name returned for the selected field or parameter.
        logical, intent(out) :: found
            !! True when the requested field or parameter metadata exists.

        call self%spec%parameter_info(index,name,found)
    end subroutine material_parameter_info

    pure integer function material_input_size(self) result(value)
        !! Total packed input components per point. Zero before initialization.
        class(material_t), intent(in) :: self
            !! Instance of `material_t` accessed by this operation.

        value = self%spec%input_size()
    end function material_input_size

    pure integer function material_output_size(self) result(value)
        !! Total packed output components per point. Zero before initialization.
        class(material_t), intent(in) :: self
            !! Instance of `material_t` accessed by this operation.

        value = self%spec%output_size()
    end function material_output_size

    pure integer function material_state_size(self) result(value)
        !! Real-valued history components per point, excluding typed state.
        class(material_t), intent(in) :: self
            !! Instance of `material_t` accessed by this operation.

        value = self%spec%state_size()
    end function material_state_size

    pure integer function material_input_offset(self, name) result(value)
        !! One-based first component of a named input. Zero if not declared.
        !! Resolve during setup, not inside numerical integration loops.
        class(material_t), intent(in) :: self
            !! Instance of `material_t` accessed by this operation.
        character(*), intent(in) :: name
            !! Exact name of the field or scalar parameter.

        value = self%spec%input_offset(name)
    end function material_input_offset

    pure integer function material_output_offset(self, name) result(value)
        !! One-based first component of a named output. Zero if not declared.
        !! Resolve during setup, not inside numerical integration loops.
        class(material_t), intent(in) :: self
            !! Instance of `material_t` accessed by this operation.
        character(*), intent(in) :: name
            !! Exact name of the field or scalar parameter.

        value = self%spec%output_offset(name)
    end function material_output_offset

    pure integer function material_input_components(self, name) result(value)
        !! Number of packed components for an input field, or zero if absent.
        class(material_t), intent(in) :: self
            !! Instance of `material_t` accessed by this operation.
        character(*), intent(in) :: name
            !! Exact name of the field or scalar parameter.

        value = self%spec%input_components(name)
    end function material_input_components

    pure integer function material_output_components(self, name) result(value)
        !! Number of packed components for an output field, or zero if absent.
        class(material_t), intent(in) :: self
            !! Instance of `material_t` accessed by this operation.
        character(*), intent(in) :: name
            !! Exact name of the field or scalar parameter.

        value = self%spec%output_components(name)
    end function material_output_components

    pure logical function material_field_matches(self, role, name, convention, measure, frame, shape) result(matches)
        !! Validate a solver-side field binding exactly once during setup.
        !! The match compares case-sensitive physical measure, frame,
        !! storage convention and physical shape. Unit consistency remains the
        !! caller's responsibility. `role` is "input", "output", "state",
        !! "integer_state" or "logical_state".
        class(material_t), intent(in) :: self
            !! Instance of `material_t` accessed by this operation.
        character(*), intent(in) :: role
            !! Schema role (`input`, `output`, `state`, `integer_state` or `logical_state`).
        character(*), intent(in) :: name
            !! Exact name of the field or scalar parameter.
        character(*), intent(in) :: convention
            !! Physical component packing convention used in exact metadata matching.
        character(*), intent(in) :: measure
            !! Identifier of the physical measure represented by the field.
        character(*), intent(in) :: frame
            !! Physical reference-frame identifier of the field.
        integer, intent(in) :: shape(:)
            !! Physical tensor shape. Use a zero-length array for scalar fields.

        matches = .false.
        if (.not. allocated(self%model)) return
        matches = self%spec%matches_field(role, name, convention, measure, frame, shape)
    end function material_field_matches

    subroutine material_field_info(self, role, name, convention, measure, frame, shape, found)
        !! Return field metadata for solver setup. Allocates returned strings
        !! and shape. It must not be called from the integration-point loop.
        !! Valid roles are "input", "output", "state", "integer_state"
        !! and "logical_state". Unknown names/roles or uninitialized handles
        !! return `found=.false.` with unallocated metadata outputs.
        class(material_t), intent(in) :: self
            !! Instance of `material_t` accessed by this operation.
        character(*), intent(in) :: role
            !! Schema role (`input`, `output`, `state`, `integer_state` or `logical_state`).
        character(*), intent(in) :: name
            !! Exact name of the field or scalar parameter.
        character(:), allocatable, intent(out) :: convention
            !! Component storage convention returned for the selected field.
        character(:), allocatable, intent(out) :: measure
            !! Physical measure identifier returned for the selected field.
        character(:), allocatable, intent(out) :: frame
            !! Frame identifier returned for the selected field.
        integer, allocatable, intent(out) :: shape(:)
            !! Physical tensor dimensions returned for the selected field.
        logical, intent(out) :: found
            !! True when the requested field or parameter metadata exists.

        found = .false.
        if (.not. allocated(self%model)) return
        call self%spec%field_info(role, name, convention, measure, frame, shape, found)
    end subroutine material_field_info

    subroutine material_point_statuses(self, codes, status)
        !! Copy the status of *each point* from the latest invoked batch kernel.
        !!
        !! codes(npoints) has caller ownership. Returns cf_ok when an
        !! evaluated batch has per-point outcomes, even if some points
        !! failed. A setup/shape/time validation failure invalidates the
        !! previous diagnostic record. In that case return cf_no_trial.
        !! Returns cf_uninitialized for an empty material and
        !! cf_invalid_shape for the wrong number of output entries.
        !! This is a setup/debugging query, not an integration hot-path
        !! callback. No allocation or string formatting is performed.
        class(material_t), intent(in) :: self
            !! Instance of `material_t` accessed by this operation.
        integer, intent(out) :: codes(:)
            !! Caller-owned output array containing the most recent point status codes.
        integer, intent(out) :: status
            !! Operation status code; `cf_ok` indicates success.

        status = cf_uninitialized
        if (.not. allocated(self%model)) return
        status = cf_invalid_shape
        if (size(codes) /= self%npoints) return
        status = cf_no_trial
        if (.not. self%point_results_ready) return
        codes = self%point_status
        status = cf_ok
    end subroutine material_point_statuses

    pure integer function last_failure_point(self) result(value)
        !! First failed point index, or zero if all recorded points succeeded.
        class(material_t), intent(in) :: self
            !! Instance of `material_t` accessed by this operation.
        integer :: i

        value = 0
        if (.not. self%point_failure) return
        if (.not. allocated(self%point_status)) return
        do i = 1, size(self%point_status)
            if (self%point_status(i) /= cf_ok) then
                value = i
                return
            end if
        end do
    end function last_failure_point
end module constituf_material

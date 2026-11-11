! SPDX-License-Identifier: Apache-2.0

module constituf_spec
    !! Immutable-after-initialization schema for packed constitutive fields.
    !!
    !! Names and metadata are resolved at setup. The numerical kernel uses
    !! prearranged component offsets. It never searches strings at runtime.
    implicit none (type, external)
    private

    public :: material_spec_t

    type :: field_t
        !! Name, physical metadata and packed component range of one field.
        character(:), allocatable :: name
            !! Exact field identifier.
        character(:), allocatable :: convention
            !! Packing convention, such as `scalar`, `mandel` or `tensor_full`.
        character(:), allocatable :: measure
            !! Identifier of the associated physical tensor or scalar measure.
        character(:), allocatable :: frame
            !! Reference frame associated with the physical field.
        integer, allocatable :: shape(:)
            !! Physical tensor dimensions, independent of packed storage.
        integer :: first = 0
            !! One-based offset of the first packed field component.
        integer :: ncomponents = 0
            !! Number of packed scalar components in this field.
    end type field_t

    type :: field_list_t
        !! Validated ordered field definitions with contiguous packed storage.
        private
        type(field_t), allocatable :: fields(:)
            !! Ordered field definitions for this schema role.
        integer :: total = 0
            !! Sum of the packed component counts in this list.
        logical :: valid = .true.
            !! False if a field declaration or validation has failed.
    contains
        procedure :: add => add_field
            !! Append a validated field to this ordered packed field list.
        procedure :: offset => field_offset
            !! Find a field's first packed component index.
        procedure :: components => field_components
            !! Find a field's number of stored components.
        procedure :: matches => field_matches
            !! Compare a named field with its declared physical metadata.
        procedure :: info => field_info
            !! Return the physical metadata of a named field.
    end type field_list_t

    type :: material_spec_t
        !! Describes one material's packed inputs, outputs and history state.
        !!
        !! Every field declares a storage convention, physical measure,
        !! frame and physical shape. The packed component count may differ
        !! from the number of tensor entries (e.g. 6 for a 3 x 3 Mandel tensor).
        !! Metadata is checked once at setup. Evaluation uses packed arrays.
        !! `tensor_full` and `array` use Fortran column-major flattening.
        !! Three-dimensional symmetric `mandel` tensors use
        !! \[\left(a_{11},a_{22},a_{33},\sqrt{2}a_{23},
        !! \sqrt{2}a_{13},\sqrt{2}a_{12}\right)\]
        !! with physical shape [3,3] and six stored components. In two
        !! dimensions the convention is
        !! \[\left(a_{11},a_{22},\sqrt{2}a_{12}\right)\]
        !! with physical shape [2,2] and three stored components.
        !! Frames may describe two-point tensors (e.g. `spatial_reference`).
        !! The schema stores no units. The caller owns unit consistency.
        !! Integer and logical state fields have distinct packed buffers.
        !! Roles `state`, `integer_state` and `logical_state` identify the
        !! storage kind for setup-time metadata lookup. Names cannot
        !! collide across history kinds.
        private
        type(field_list_t) :: inputs
            !! Ordered input field definitions.
        type(field_list_t) :: outputs
            !! Ordered response field definitions.
        type(field_list_t) :: states
            !! Ordered real history field definitions.
        type(field_list_t) :: integer_states
            !! Ordered integer history field definitions.
        type(field_list_t) :: logical_states
            !! Ordered logical history field definitions.
        type(field_list_t) :: parameters
            !! Ordered scalar parameter definitions.
        logical :: has_sensitivity = .false.
            !! True if the model supplies exact parameter derivatives.
        logical :: needs_dt = .false.
            !! True if evaluation requires a strictly positive time increment.
        logical :: needs_initial_input = .false.
            !! True if the initial state requires prescribed inputs.
        logical :: needs_previous_input = .false.
            !! True if the preceding input must be retained.
        logical :: has_tangent = .false.
            !! True if the model supplies an exact input Jacobian.
    contains
        procedure :: add_input
            !! Append a packed input field to the material schema.
        procedure :: add_output
            !! Append a packed response field to the material schema.
        procedure :: add_state
            !! Append a real-valued material history field.
        procedure :: add_integer_state
            !! Append an integer-valued history field.
        procedure :: add_logical_state
            !! Append a logical history field.
        procedure :: add_parameter
            !! Register a named scalar parameter.
        procedure :: set_has_sensitivity
            !! Declare whether exact parameter sensitivities are supported.
        procedure :: parameter_count
            !! Return the number of declared scalar material parameters.
        procedure :: parameter_id
            !! Return the one-based ID of a named scalar parameter, or zero if it is absent.
        procedure :: parameter_info
            !! Get the name of a parameter by one-based index.
        procedure :: supports_sensitivity
            !! Test support for analytical parameter sensitivities.
        procedure :: set_requires_dt
            !! Declare whether positive time increments are required.
        procedure :: set_requires_initial_input
            !! Declare whether initial inputs must be supplied.
        procedure :: set_requires_previous_input
            !! Declare whether the previous input must be retained.
        procedure :: set_has_tangent
            !! Declare whether an exact input Jacobian is available.
        procedure :: input_size
            !! Return the number of packed input components per material point.
        procedure :: output_size
            !! Return the number of packed response components per point.
        procedure :: state_size
            !! Return the number of real-valued history components per point.
        procedure :: integer_state_size
            !! Return the number of integer history components per point.
        procedure :: logical_state_size
            !! Return the number of logical history components per point.
        procedure :: input_offset
            !! Find a named input field's first packed index.
        procedure :: output_offset
            !! Find a named output field's first packed index.
        procedure :: state_offset
            !! Find a named real history field's first packed index.
        procedure :: integer_state_offset
            !! Find a named integer history field's first packed index.
        procedure :: logical_state_offset
            !! Find a named logical history field's first packed index.
        procedure :: input_components
            !! Return a named input field's packed component count.
        procedure :: output_components
            !! Return a named output field's packed component count.
        procedure :: state_components
            !! Return a named real state field's component count.
        procedure :: matches_field => spec_matches_field
            !! Test whether a named field matches its registered role and physical metadata.
        procedure :: field_info => spec_field_info
            !! Retrieve metadata for a named field of the specified role.
        procedure :: requires_dt
            !! Test whether time increments must be supplied.
        procedure :: requires_initial_input
            !! Test whether initial input data are required.
        procedure :: requires_previous_input
            !! Test whether the preceding input must be available.
        procedure :: supports_tangent
            !! Test whether the model provides an exact algorithmic tangent.
        procedure :: valid => valid_spec
            !! Return whether all schema registrations and metadata are valid.
    end type material_spec_t

contains

    subroutine add_field(self, name, ncomponents, convention, measure, frame, shape)
        !! Register a physical field. All semantic declarations are mandatory.
        !! `shape` is required unless `convention="scalar"`, which has rank zero.
        !! A Mandel 3D symmetric tensor uses `shape=[3,3]`, `ncomponents=6`.
        class(field_list_t), intent(inout) :: self
            !! Instance of `field_list_t` accessed by this operation.
        character(*), intent(in) :: name
            !! Exact name of the field or scalar parameter.
        character(*), intent(in) :: convention
            !! Physical component packing convention used in exact metadata matching.
        character(*), intent(in) :: measure
            !! Identifier of the physical measure represented by the field.
        character(*), intent(in) :: frame
            !! Physical reference-frame identifier of the field.
        integer, intent(in) :: ncomponents
            !! Number of scalar components in the field's packed representation.
        integer, intent(in), optional :: shape(:)
            !! Physical tensor shape. omitted for scalar fields.
        type(field_t), allocatable :: next(:)
        integer :: n
        integer :: i
        integer :: product_size
        integer :: allocation_status
        character(256) :: allocation_message
        logical :: supported

        if (.not. self%valid) return
        self%valid = .false.
        if (len_trim(name) == 0 .or. len_trim(convention) == 0 .or. len_trim(measure) == 0 .or. len_trim(frame) == 0 .or. ncomponents <= 0) return
        if (verify(trim(name), "abcdefghijklmnopqrstuvwxyz0123456789_.") /= 0) return
        if (verify(trim(measure), "abcdefghijklmnopqrstuvwxyz0123456789_.") /= 0) return
        if (verify(trim(frame), "abcdefghijklmnopqrstuvwxyz0123456789_.") /= 0) return
        if (verify(trim(convention), "abcdefghijklmnopqrstuvwxyz0123456789_.") /= 0) return
        if (self%offset(name) /= 0) return
        if (self%total > huge(self%total) - ncomponents) return

        supported = .false.
        select case (convention)
        case ("scalar")
            supported = ncomponents == 1
            if (present(shape)) supported = supported .and. size(shape) == 0
        case ("vector", "tensor_full", "mandel", "array")
            if (.not. present(shape)) return
            if (any(shape <= 0)) return
            if (size(shape) > 6) return
            product_size = 1
            do i = 1, size(shape)
                if (product_size > huge(product_size) / shape(i)) return
                product_size = product_size * shape(i)
            end do
            select case (convention)
            case ("vector")
                supported = size(shape) == 1 .and. ncomponents == product_size
            case ("tensor_full")
                supported = size(shape) == 2 .and. ncomponents == product_size
            case ("array")
                supported = size(shape) > 0 .and. ncomponents == product_size
            case ("mandel")
                if (size(shape) /= 2) return
                supported = (all(shape == [2,2]) .and. ncomponents == 3) .or. (all(shape == [3,3]) .and. ncomponents == 6)
            case default
                supported = .false.
            end select
        case default
            ! Custom packed representations must be explicit and namespaced.
            if (index(convention, "custom.") /= 1) return
            if (len_trim(convention) <= len("custom.")) return
            if (.not. present(shape)) return
            if (size(shape) > 6) return
            if (any(shape <= 0)) return
            supported = .true.
        end select
        if (.not. supported) return

        n = 0
        if (allocated(self%fields)) n = size(self%fields)
        allocate(next(n + 1), stat=allocation_status, errmsg=allocation_message)
        if (allocation_status /= 0) return
        if (n > 0) next(:n) = self%fields
        next(n + 1)%name = trim(name)
        next(n + 1)%convention = trim(convention)
        next(n + 1)%measure = trim(measure)
        next(n + 1)%frame = trim(frame)
        if (present(shape)) then
            next(n + 1)%shape = shape
        else
            allocate(next(n + 1)%shape(0))
        end if
        next(n + 1)%first = self%total + 1
        next(n + 1)%ncomponents = ncomponents
        call move_alloc(next, self%fields)
        self%total = self%total + ncomponents
        self%valid = .true.
    end subroutine add_field

    pure integer function field_offset(self, name) result(offset)
        !! Return a field's one-based first packed component index, or zero if absent.
        class(field_list_t), intent(in) :: self
            !! Instance of `field_list_t` accessed by this operation.
        character(*), intent(in) :: name
            !! Exact name of the field or scalar parameter.
        integer :: i

        offset = 0
        if (.not. allocated(self%fields)) return
        do i = 1, size(self%fields)
            if (self%fields(i)%name == name) then
                offset = self%fields(i)%first
                return
            end if
        end do
    end function field_offset

    pure integer function field_components(self, name) result(ncomponents)
        !! Return a field's packed component count, or zero if absent.
        class(field_list_t), intent(in) :: self
            !! Instance of `field_list_t` accessed by this operation.
        character(*), intent(in) :: name
            !! Exact name of the field or scalar parameter.
        integer :: i

        ncomponents = 0
        if (.not. allocated(self%fields)) return
        do i = 1, size(self%fields)
            if (self%fields(i)%name == name) then
                ncomponents = self%fields(i)%ncomponents
                return
            end if
        end do
    end function field_components

    subroutine add_input(self, name, ncomponents, convention, measure, frame, shape)
        !! Register one input field in packed order. Fields are immutable after init.
        !! Duplicate names, invalid metadata or incompatible shape invalidate the spec.
        class(material_spec_t), intent(inout) :: self
            !! Instance of `material_spec_t` accessed by this operation.
        character(*), intent(in) :: name
            !! Exact name of the field or scalar parameter.
        character(*), intent(in) :: convention
            !! Physical component packing convention used in exact metadata matching.
        character(*), intent(in) :: measure
            !! Identifier of the physical measure represented by the field.
        character(*), intent(in) :: frame
            !! Physical reference-frame identifier of the field.
        integer, intent(in) :: ncomponents
            !! Number of scalar components in the field's packed representation.
        integer, intent(in), optional :: shape(:)
            !! Physical tensor shape. omitted for scalar fields.

        call self%inputs%add(name, ncomponents, convention, measure, frame, shape)
    end subroutine add_input

    subroutine add_output(self, name, ncomponents, convention, measure, frame, shape)
        !! Register one output field in packed order, with the same semantic
        !! metadata and validation rules as `add_input`.
        class(material_spec_t), intent(inout) :: self
            !! Instance of `material_spec_t` accessed by this operation.
        character(*), intent(in) :: name
            !! Exact name of the field or scalar parameter.
        character(*), intent(in) :: convention
            !! Physical component packing convention used in exact metadata matching.
        character(*), intent(in) :: measure
            !! Identifier of the physical measure represented by the field.
        character(*), intent(in) :: frame
            !! Physical reference-frame identifier of the field.
        integer, intent(in) :: ncomponents
            !! Number of scalar components in the field's packed representation.
        integer, intent(in), optional :: shape(:)
            !! Physical tensor shape. omitted for scalar fields.

        call self%outputs%add(name, ncomponents, convention, measure, frame, shape)
    end subroutine add_output

    subroutine add_state(self, name, ncomponents, convention, measure, frame, shape)
        !! Register real-valued committed/trial history components per point.
        !! Names must not collide with integer or logical history fields.
        class(material_spec_t), intent(inout) :: self
            !! Instance of `material_spec_t` accessed by this operation.
        character(*), intent(in) :: name
            !! Exact name of the field or scalar parameter.
        character(*), intent(in) :: convention
            !! Physical component packing convention used in exact metadata matching.
        character(*), intent(in) :: measure
            !! Identifier of the physical measure represented by the field.
        character(*), intent(in) :: frame
            !! Physical reference-frame identifier of the field.
        integer, intent(in) :: ncomponents
            !! Number of scalar components in the field's packed representation.
        integer, intent(in), optional :: shape(:)
            !! Physical tensor shape. omitted for scalar fields.

        if (self%integer_states%offset(name) /= 0 .or. self%logical_states%offset(name) /= 0) then
            self%states%valid = .false.
            return
        end if
        call self%states%add(name, ncomponents, convention, measure, frame, shape)
    end subroutine add_state

    subroutine add_integer_state(self, name, ncomponents, convention, measure, frame, shape)
        !! Declare a named `integer(ik)` history field, packed independently of reals.
        !! Initial values are set by `initialize_typed_state` (zero by default).
        class(material_spec_t), intent(inout) :: self
            !! Instance of `material_spec_t` accessed by this operation.
        character(*), intent(in) :: name
            !! Exact name of the field or scalar parameter.
        character(*), intent(in) :: convention
            !! Physical component packing convention used in exact metadata matching.
        character(*), intent(in) :: measure
            !! Identifier of the physical measure represented by the field.
        character(*), intent(in) :: frame
            !! Physical reference-frame identifier of the field.
        integer, intent(in) :: ncomponents
            !! Number of scalar components in the field's packed representation.
        integer, intent(in), optional :: shape(:)
            !! Physical tensor shape. omitted for scalar fields.

        if (self%states%offset(name) /= 0 .or. self%logical_states%offset(name) /= 0) then
            self%integer_states%valid = .false.
            return
        end if
        call self%integer_states%add(name, ncomponents, convention, measure, frame, shape)
    end subroutine add_integer_state

    subroutine add_logical_state(self, name, ncomponents, convention, measure, frame, shape)
        !! Declare a named logical history field, separately packed per point.
        !! No classical derivative is defined with respect to a discrete flag.
        class(material_spec_t), intent(inout) :: self
            !! Instance of `material_spec_t` accessed by this operation.
        character(*), intent(in) :: name
            !! Exact name of the field or scalar parameter.
        character(*), intent(in) :: convention
            !! Physical component packing convention used in exact metadata matching.
        character(*), intent(in) :: measure
            !! Identifier of the physical measure represented by the field.
        character(*), intent(in) :: frame
            !! Physical reference-frame identifier of the field.
        integer, intent(in) :: ncomponents
            !! Number of scalar components in the field's packed representation.
        integer, intent(in), optional :: shape(:)
            !! Physical tensor shape. omitted for scalar fields.

        if (self%states%offset(name) /= 0 .or. self%integer_states%offset(name) /= 0) then
            self%logical_states%valid = .false.
            return
        end if
        call self%logical_states%add(name, ncomponents, convention, measure, frame, shape)
    end subroutine add_logical_state

    pure logical function field_matches(self, name, convention, measure, frame, shape) result(matches)
        !! Test whether a field's stored metadata matches the requested description.
        class(field_list_t), intent(in) :: self
            !! Instance of `field_list_t` accessed by this operation.
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
        integer :: i

        matches = .false.
        if (.not. allocated(self%fields)) return
        do i = 1, size(self%fields)
            if (self%fields(i)%name == name) then
                if (self%fields(i)%convention /= convention) return
                if (self%fields(i)%measure /= measure) return
                if (self%fields(i)%frame /= frame) return
                if (size(self%fields(i)%shape) /= size(shape)) return
                if (any(self%fields(i)%shape /= shape)) return
                matches = .true.
                return
            end if
        end do
    end function field_matches

    subroutine field_info(self, name, convention, measure, frame, shape, found)
        !! Retrieve a field's physical metadata and report whether it exists.
        class(field_list_t), intent(in) :: self
            !! Instance of `field_list_t` accessed by this operation.
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
        integer :: i

        found = .false.
        if (.not. allocated(self%fields)) return
        do i = 1, size(self%fields)
            if (self%fields(i)%name == name) then
                convention = self%fields(i)%convention
                measure = self%fields(i)%measure
                frame = self%fields(i)%frame
                shape = self%fields(i)%shape
                found = .true.
                return
            end if
        end do
    end subroutine field_info

    pure logical function spec_matches_field(self, role, name, convention, measure, frame, shape) result(matches)
        !! Compare the declared field role, measure, frame, packing and shape.
        !! Use at solver setup. No string matching occurs in `eval_batch`.
        class(material_spec_t), intent(in) :: self
            !! Instance of `material_spec_t` accessed by this operation.
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
        select case (role)
        case ("input")
            matches = self%inputs%matches(name, convention, measure, frame, shape)
        case ("output")
            matches = self%outputs%matches(name, convention, measure, frame, shape)
        case ("state")
            matches = self%states%matches(name, convention, measure, frame, shape)
        case ("integer_state")
            matches = self%integer_states%matches(name, convention, measure, frame, shape)
        case ("logical_state")
            matches = self%logical_states%matches(name, convention, measure, frame, shape)
        case default
            matches = .false.
        end select
    end function spec_matches_field

    subroutine spec_field_info(self, role, name, convention, measure, frame, shape, found)
        !! Query metadata during setup, without exposing extra derived types.
        !! For unknown fields, `found=.false.` and outputs stay unallocated.
        class(material_spec_t), intent(in) :: self
            !! Instance of `material_spec_t` accessed by this operation.
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
        select case (role)
        case ("input")
            call self%inputs%info(name, convention, measure, frame, shape, found)
        case ("output")
            call self%outputs%info(name, convention, measure, frame, shape, found)
        case ("state")
            call self%states%info(name, convention, measure, frame, shape, found)
        case ("integer_state")
            call self%integer_states%info(name, convention, measure, frame, shape, found)
        case ("logical_state")
            call self%logical_states%info(name, convention, measure, frame, shape, found)
        case default
            found = .false.
        end select
    end subroutine spec_field_info

    subroutine add_parameter(self, name)
        !! Declare a scalar differentiable parameter with a stable one-based
        !! index. Only scalar parameters are supported in this first release.
        class(material_spec_t), intent(inout) :: self
            !! Instance of `material_spec_t` accessed by this operation.
        character(*), intent(in) :: name
            !! Exact name of the field or scalar parameter.

        call self%parameters%add(name, 1, "scalar", "model_parameter", "none")
    end subroutine add_parameter

    subroutine set_has_sensitivity(self, value)
        !! Opt-in: a model must implement exact discrete path sensitivities.
        class(material_spec_t), intent(inout) :: self
            !! Instance of `material_spec_t` accessed by this operation.
        logical, intent(in) :: value
            !! Value to register or report in the field schema.

        self%has_sensitivity = value
    end subroutine set_has_sensitivity

    pure integer function parameter_count(self) result(value)
        !! Number of declared scalar differentiable parameters.
        class(material_spec_t), intent(in) :: self
            !! Instance of `material_spec_t` accessed by this operation.

        value = self%parameters%total
    end function parameter_count

    pure integer function parameter_id(self, name) result(value)
        !! Resolve names during setup, never in evaluation loops.
        class(material_spec_t), intent(in) :: self
            !! Instance of `material_spec_t` accessed by this operation.
        character(*), intent(in) :: name
            !! Exact name of the field or scalar parameter.

        value = self%parameters%offset(name)
    end function parameter_id

    subroutine parameter_info(self, index, name, found)
        !! Return the name of one scalar parameter in declared order.
        !! One-based index. Absent indices return an empty name and false.
        !! This is a setup-time metadata lookup, never a hot-path operation.
        class(material_spec_t), intent(in) :: self
            !! Instance of `material_spec_t` accessed by this operation.
        integer, intent(in) :: index
            !! One-based index identifying the selected entry.
        character(:), allocatable, intent(out) :: name
            !! Name returned for the selected field or parameter.
        logical, intent(out) :: found
            !! True when the requested field or parameter metadata exists.

        name=""
        found=.false.
        if (.not. allocated(self%parameters%fields)) return
        if (index<1.or.index>size(self%parameters%fields)) return
        name=self%parameters%fields(index)%name
        found=.true.
    end subroutine parameter_info

    pure logical function supports_sensitivity(self) result(value)
        !! True only when the author declared analytical path sensitivity support.
        class(material_spec_t), intent(in) :: self
            !! Instance of `material_spec_t` accessed by this operation.

        value = self%has_sensitivity
    end function supports_sensitivity

    subroutine set_requires_dt(self, value)
        !! Require a strictly positive finite `dt` in each trial evaluation.
        class(material_spec_t), intent(inout) :: self
            !! Instance of `material_spec_t` accessed by this operation.
        logical, intent(in) :: value
            !! Value to register or report in the field schema.

        self%needs_dt = value
    end subroutine set_requires_dt

    subroutine set_requires_initial_input(self, value)
        !! Require `initial_input(ninput,npoints)` when initializing the handle.
        class(material_spec_t), intent(inout) :: self
            !! Instance of `material_spec_t` accessed by this operation.
        logical, intent(in) :: value
            !! Value to register or report in the field schema.

        self%needs_initial_input = value
    end subroutine set_requires_initial_input

    subroutine set_requires_previous_input(self, value)
        !! Keep the last committed input for use in incremental constitutive updates.
        class(material_spec_t), intent(inout) :: self
            !! Instance of `material_spec_t` accessed by this operation.
        logical, intent(in) :: value
            !! Value to register or report in the field schema.

        self%needs_previous_input = value
    end subroutine set_requires_previous_input

    subroutine set_has_tangent(self, value)
        !! Advertise an exact current-input response Jacobian.
        class(material_spec_t), intent(inout) :: self
            !! Instance of `material_spec_t` accessed by this operation.
        logical, intent(in) :: value
            !! Value to register or report in the field schema.

        self%has_tangent = value
    end subroutine set_has_tangent

    pure integer function input_size(self) result(value)
        !! Total packed input components per material point.
        class(material_spec_t), intent(in) :: self
            !! Instance of `material_spec_t` accessed by this operation.

        value = self%inputs%total
    end function input_size

    pure integer function output_size(self) result(value)
        !! Total packed output components per material point.
        class(material_spec_t), intent(in) :: self
            !! Instance of `material_spec_t` accessed by this operation.

        value = self%outputs%total
    end function output_size

    pure integer function state_size(self) result(value)
        !! Total real history components per material point.
        class(material_spec_t), intent(in) :: self
            !! Instance of `material_spec_t` accessed by this operation.

        value = self%states%total
    end function state_size

    pure integer function integer_state_size(self) result(value)
        !! Integer history components per material point, not including real state.
        class(material_spec_t), intent(in) :: self
            !! Instance of `material_spec_t` accessed by this operation.

        value = self%integer_states%total
    end function integer_state_size

    pure integer function logical_state_size(self) result(value)
        !! Logical history components per material point, not including real state.
        class(material_spec_t), intent(in) :: self
            !! Instance of `material_spec_t` accessed by this operation.

        value = self%logical_states%total
    end function logical_state_size

    pure integer function input_offset(self, name) result(value)
        !! One-based first component of a named input, or zero if absent.
        class(material_spec_t), intent(in) :: self
            !! Instance of `material_spec_t` accessed by this operation.
        character(*), intent(in) :: name
            !! Exact name of the field or scalar parameter.

        value = self%inputs%offset(name)
    end function input_offset

    pure integer function output_offset(self, name) result(value)
        !! One-based first component of a named output, or zero if absent.
        class(material_spec_t), intent(in) :: self
            !! Instance of `material_spec_t` accessed by this operation.
        character(*), intent(in) :: name
            !! Exact name of the field or scalar parameter.

        value = self%outputs%offset(name)
    end function output_offset

    pure integer function state_offset(self, name) result(value)
        !! One-based first real-state component, or zero if absent.
        class(material_spec_t), intent(in) :: self
            !! Instance of `material_spec_t` accessed by this operation.
        character(*), intent(in) :: name
            !! Exact name of the field or scalar parameter.

        value = self%states%offset(name)
    end function state_offset

    pure integer function integer_state_offset(self, name) result(value)
        !! One-based first integer-state component, or zero if absent.
        class(material_spec_t), intent(in) :: self
            !! Instance of `material_spec_t` accessed by this operation.
        character(*), intent(in) :: name
            !! Exact name of the field or scalar parameter.

        value = self%integer_states%offset(name)
    end function integer_state_offset

    pure integer function logical_state_offset(self, name) result(value)
        !! One-based first logical-state component, or zero if absent.
        class(material_spec_t), intent(in) :: self
            !! Instance of `material_spec_t` accessed by this operation.
        character(*), intent(in) :: name
            !! Exact name of the field or scalar parameter.

        value = self%logical_states%offset(name)
    end function logical_state_offset

    pure integer function input_components(self, name) result(value)
        !! Number of packed input components, or zero if the name is unknown.
        class(material_spec_t), intent(in) :: self
            !! Instance of `material_spec_t` accessed by this operation.
        character(*), intent(in) :: name
            !! Exact name of the field or scalar parameter.

        value = self%inputs%components(name)
    end function input_components

    pure integer function output_components(self, name) result(value)
        !! Number of packed output components, or zero if the name is unknown.
        class(material_spec_t), intent(in) :: self
            !! Instance of `material_spec_t` accessed by this operation.
        character(*), intent(in) :: name
            !! Exact name of the field or scalar parameter.

        value = self%outputs%components(name)
    end function output_components

    pure integer function state_components(self, name) result(value)
        !! Number of real history components, or zero if the name is unknown.
        class(material_spec_t), intent(in) :: self
            !! Instance of `material_spec_t` accessed by this operation.
        character(*), intent(in) :: name
            !! Exact name of the field or scalar parameter.

        value = self%states%components(name)
    end function state_components

    pure logical function requires_dt(self) result(value)
        !! Whether this material requires a positive time increment.
        class(material_spec_t), intent(in) :: self
            !! Instance of `material_spec_t` accessed by this operation.

        value = self%needs_dt
    end function requires_dt

    pure logical function requires_initial_input(self) result(value)
        !! Whether initialization requires a prescribed reference input.
        class(material_spec_t), intent(in) :: self
            !! Instance of `material_spec_t` accessed by this operation.

        value = self%needs_initial_input
    end function requires_initial_input

    pure logical function requires_previous_input(self) result(value)
        !! Whether the committed prior input is retained between steps.
        class(material_spec_t), intent(in) :: self
            !! Instance of `material_spec_t` accessed by this operation.

        value = self%needs_previous_input
    end function requires_previous_input

    pure logical function supports_tangent(self) result(value)
        !! Whether the model supports the exact current-input Jacobian.
        class(material_spec_t), intent(in) :: self
            !! Instance of `material_spec_t` accessed by this operation.

        value = self%has_tangent
    end function supports_tangent

    pure logical function valid_spec(self) result(value)
        !! True if field/parameter declarations are consistent and both inputs
        !! and outputs exist. Sensitivity-enabled models also require parameters.
        class(material_spec_t), intent(in) :: self
            !! Instance of `material_spec_t` accessed by this operation.

        value = self%inputs%valid .and. self%outputs%valid .and. self%states%valid .and. self%integer_states%valid .and. self%logical_states%valid .and. self%parameters%valid .and. self%inputs%total > 0 .and. self%outputs%total > 0
        if (self%has_sensitivity) then
            value = value .and. self%parameters%total > 0
        end if
    end function valid_spec
end module constituf_spec

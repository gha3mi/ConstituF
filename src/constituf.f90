! SPDX-License-Identifier: Apache-2.0

module constituf
    !! Public entry point for constitutive material-point evaluation.
    !!
    !! Import `material_t`, `wp` and status constants from this module.
    !! Import a concrete material constructor from its own
    !! `constituf_<model>` module. The core has no dependency on concrete
    !! material modules.
    !!
    !! The numerical API carries no units. The caller chooses one consistent
    !! unit system for material parameters, inputs, time increments and history.
    !! Responses and derivatives use the corresponding derived units.
    !!
    !! `material_t` owns accepted and trial history for a batch of points.
    !! Call `init` to configure a model, `eval` to compute a trial response,
    !! `commit` to accept it, or `revert` to discard it. Packed arrays
    !! use `(component,point)` layout. Optional tangents and analytical
    !! parameter sensitivities depend on the selected material.
    !!
    !! The module also exports optional verification, calibration and
    !! checkpoint utilities. For complete array and lifecycle contracts,
    !! see [[constituf_material(module):material_t(type)]].
    use constituf_calibration, only: calibrate_parameters
    use constituf_checkpoint, only: save_checkpoint, load_checkpoint
    use constituf_kinds, only: wp, ik, lk
    use constituf_material, only: material_t
    use constituf_potential_verify, only: verify_potential
    use constituf_status, only: cf_ok, cf_invalid_argument, cf_invalid_shape, cf_uninitialized, cf_no_trial, cf_invalid_time, cf_unsupported, cf_nondifferentiable, cf_evaluation_failed, cf_allocation_failed, cf_checkpoint_mismatch, cf_checkpoint_corrupt, cf_io_error
    use constituf_verify, only: verify_tangent
    implicit none (type, external)
    private

    public :: wp, ik, lk, material_t, verify_tangent, verify_potential
    public :: save_checkpoint, load_checkpoint, calibrate_parameters
    public :: cf_ok, cf_invalid_argument, cf_invalid_shape
    public :: cf_uninitialized, cf_no_trial, cf_invalid_time
    public :: cf_unsupported, cf_nondifferentiable, cf_evaluation_failed, cf_allocation_failed
    public :: cf_checkpoint_mismatch, cf_checkpoint_corrupt, cf_io_error
end module constituf

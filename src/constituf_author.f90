! SPDX-License-Identifier: Apache-2.0

module constituf_author
    !! Single import for authors of independent constitutive models.
    !!
    !! Implement `material_model_t%describe` and `%eval_batch`. All user-facing
    !! state management, buffer ownership and commit/revert are provided by the
    !! core. This module has no dependency on concrete materials.
    !! Field declarations describe physical measures, frames and packing.
    !! Document dimensional relationships in each model. Numerical values use
    !! the caller's consistent unit system without unit metadata or conversion.
    use constituf_kinds, only: wp, ik
    use constituf_model, only: material_model_t
    use constituf_spec, only: material_spec_t
    use constituf_status, only: cf_ok, cf_invalid_argument, cf_invalid_shape, cf_invalid_time, cf_unsupported, cf_nondifferentiable, cf_evaluation_failed, cf_allocation_failed
    implicit none (type, external)
    private

    public :: wp, ik, material_model_t, material_spec_t
    public :: cf_ok, cf_invalid_argument, cf_invalid_shape
    public :: cf_invalid_time, cf_unsupported, cf_nondifferentiable, cf_evaluation_failed
    public :: cf_allocation_failed
end module constituf_author

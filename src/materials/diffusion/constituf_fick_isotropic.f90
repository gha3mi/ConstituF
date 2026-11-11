! SPDX-License-Identifier: Apache-2.0

module constituf_fick_isotropic
    !! Isotropic single-species Fickian molar diffusion in three dimensions.
    !!
    !! A prescribed spatial molar concentration \(c\), with dimensions of
    !! amount per volume, defines the gradient \(\boldsymbol{g}=\nabla c\).
    !! The gradient has concentration-per-length units. Diffusive molar
    !! flux \(\boldsymbol{j}\), in amount per area per time, obeys Fick's first law:
    !! \[\boldsymbol{j}=-D\boldsymbol{g},\qquad D>0.\]
    !! Here \(D\) is a finite, constant, strictly positive scalar
    !! diffusivity with dimensions of length squared per time.
    !! Exact local derivatives at fixed prescribed gradient:
    !! \[\frac{\partial j_i}{\partial g_j}=-D\delta_{ij},\qquad
    !! \frac{\partial\boldsymbol{j}}{\partial D}=-\boldsymbol{g}.\]
    !! The down-gradient transport inequality is
    !! \[-\boldsymbol{g}\cdot\boldsymbol{j}
    !! =D\|\boldsymbol{g}\|^2\geq0.\]
    !! This is not an entropy production rate unless a chemical potential
    !! law and associated thermodynamic assumptions are supplied separately.
    !!
    !! Input `concentration.gradient` (concentration per length) and output
    !! `species.molar_flux` (amount per area per time) are spatial 3-vectors in
    !! Cartesian component order \([1,2,3]\). One material handle represents
    !! one prescribed species. The caller supplies the gradient, valid
    !! concentrations and all initial/boundary conditions. There is no
    !! implicit unit/frame conversion or mass-conservation equation solve.
    !! The returned diffusive flux is relative to the carrier or medium.
    !! It excludes the advective contribution \(c\boldsymbol{v}\) from carrier velocity \(\boldsymbol{v}\).
    !! Positive components denote diffusion toward positive spatial axes.
    !!
    !! The model is stateless, isothermal, dilute and nonreactive, without
    !! time-step dependence, convection, concentration dependence,
    !! electromigration, cross-diffusion or finite-rate memory. Local
    !! derivatives do not include a parameter-dependent global solution.
    !! Nonfinite input/response makes the trial invalid. Point loops are
    !! allocation-free and I/O-free.
    !!
    !! Scientific law: A. Fick, *Ueber Diffusion*, Annalen der Physik
    !! 170(1), 59-86 (1855), DOI 10.1002/andp.18551700105.
    !! Modern molar-flux convention: IUPAC Gold Book, "diffusion coefficient",
    !! D01719, DOI 10.1351/goldbook.D01719:
    !! https://goldbook.iupac.org/terms/view/D01719 .
    !! Code provenance: this 3D molar specialization and its analytical
    !! derivatives are independently derived from the stated law. No
    !! third-party software code is reused in this revision.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: ndim = 3

    public :: fick_isotropic_t

    type, extends(material_model_t) :: fick_isotropic_t
        !! Spatial single-species isotropic diffusion with constant diffusivity.
        !!
        !! Construct as `fick_isotropic_t(diffusivity)` and pass to
        !! `material_t%init`. The material handle clones the coefficient.
        !! No conservation-equation solver, reaction or finite-rate memory.
        real(wp) :: diffusivity = 0.0_wp
            !! Finite, strictly positive diffusivity \(D\), with dimensions of length squared per time.
    contains
        procedure :: describe => describe_fick
            !! Declare the constitutive fields, parameter metadata and supported capabilities.
        procedure :: eval_batch => evaluate_fick
            !! Evaluate the batched constitutive response, optional tangent and trial history.
        procedure :: eval_sensitivity_batch => sensitivity_fick
            !! Evaluate exact parameter-path derivatives and corresponding trial history.
    end type fick_isotropic_t

contains

    subroutine describe_fick(self, spec, status)
        !! Declare a 3D concentration-gradient input and molar-flux response.
        !!
        !! Parameter `diffusivity` is scalar, length squared per time. Stateless, without
        !! required `dt`. Exact tangent and parameter sensitivity supported.
        class(fick_isotropic_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema receiving field layouts, parameters and supported operations.
        integer, intent(out) :: status
            !! Operation status. `cf_ok` indicates a valid model or successful initialization.

        status = cf_invalid_argument
        if (.not. ieee_is_finite(self%diffusivity)) return
        if (self%diffusivity <= 0.0_wp) return

        call spec%add_input("concentration.gradient", ndim, "vector", "concentration_gradient", "spatial", [ndim])
        call spec%add_output("species.molar_flux", ndim, "vector", "molar_flux", "spatial", [ndim])
        call spec%add_parameter("diffusivity")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_fick

    subroutine evaluate_fick(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Evaluate \(\boldsymbol{j}=-D\nabla c\) without persistent state.
        !!
        !! The caller owns `x_new(ndim,npoints)`, `response(ndim,npoints)`
        !! and optional `tangent(ndim,ndim,npoints)`. A failed point makes the
        !! entire trial invalid. Response buffers are unspecified on failure.
        class(fick_isotropic_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed input `(ninput,npoints)` if required, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `concentration.gradient` input `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued state `(nstate,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Trial real-valued state `(nstate,npoints)` returned by this evaluation.
        real(wp), intent(out) :: response(:,:)
            !! Packed `species.molar_flux` response `(noutput,npoints)` for the evaluated material points.
        real(wp), intent(in) :: dt
            !! Prescribed time increment. It is positive when the material requires time integration.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact discrete response Jacobian `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(inout) :: point_status(:)
            !! Per-point status values; successful evaluations must set each entry to `cf_ok`.

        state_trial = state_old
        call evaluate_points(self, x_new, response, point_status, tangent)
    end subroutine evaluate_fick

    subroutine sensitivity_fick(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Exact derivative with respect to the selected diffusivity parameter.
        !!
        !! `dsdp(ndim,1,npoints)` has concentration-gradient units and the
        !! zero-component sensitivity state remains unchanged.
        class(fick_isotropic_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed input `(ninput,npoints)` if required, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `concentration.gradient` input `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued state `(nstate,npoints)`.
        integer, intent(in) :: selected_parameters(:)
            !! One-based IDs of selected model parameters, in sensitivity-axis order.
        real(wp), intent(in) :: sensitivity_old(:,:,:)
            !! Committed real-state derivatives `(nstate,nselected,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Trial real-valued state `(nstate,npoints)` returned by this evaluation.
        real(wp), intent(out) :: sensitivity_trial(:,:,:)
            !! Trial real-state derivatives `(nstate,nselected,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `species.molar_flux` response `(noutput,npoints)` for the evaluated material points.
        real(wp), intent(out) :: dsdp(:,:,:)
            !! Exact total response derivatives `(noutput,nselected,npoints)` at fixed prescribed input history.
        real(wp), intent(in) :: dt
            !! Prescribed time increment. It is positive when the material requires time integration.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact discrete response Jacobian `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(inout) :: point_status(:)
            !! Per-point status values; successful evaluations must set each entry to `cf_ok`.

        state_trial = state_old
        sensitivity_trial = sensitivity_old
        call evaluate_points(self, x_new, response, point_status, tangent, selected_parameters, dsdp)
    end subroutine sensitivity_fick

    subroutine evaluate_points(self, x_new, response, point_status, tangent, selected_parameters, dsdp)
        !! Evaluate Fick flux and optional exact derivatives without allocation.
        !! Failed points retain non-success status. Caller discards the trial.
        class(fick_isotropic_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `concentration.gradient` input `(ninput,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `species.molar_flux` response `(noutput,npoints)` for the evaluated material points.
        integer, intent(inout) :: point_status(:)
            !! Per-point status values; successful evaluations must set each entry to `cf_ok`.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact discrete response Jacobian `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(in), optional :: selected_parameters(:)
            !! One-based IDs of selected model parameters, in sensitivity-axis order.
        real(wp), intent(out), optional :: dsdp(:,:,:)
            !! Exact total response derivatives `(noutput,nselected,npoints)` at fixed prescribed input history.
        integer :: i
        integer :: j
        integer :: ipoint

        point_loop: do ipoint = 1, size(x_new,2)
            point_status(ipoint) = cf_invalid_argument
            if (.not. all(ieee_is_finite(x_new(:,ipoint)))) cycle point_loop
            response(:,ipoint) = -self%diffusivity*x_new(:,ipoint)
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop

            if (present(tangent)) then
                tangent(:,:,ipoint) = 0.0_wp
                do i = 1, ndim
                    tangent(i,i,ipoint) = -self%diffusivity
                end do
            end if

            if (present(selected_parameters)) then
                parameter_loop: do j = 1, size(selected_parameters)
                    select case (selected_parameters(j))
                    case (1)
                        dsdp(:,j,ipoint) = -x_new(:,ipoint)
                    case default
                        exit parameter_loop
                    end select
                end do parameter_loop
                if (j <= size(selected_parameters)) cycle point_loop
            end if
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_points

end module constituf_fick_isotropic

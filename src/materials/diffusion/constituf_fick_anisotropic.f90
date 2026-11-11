! SPDX-License-Identifier: Apache-2.0

module constituf_fick_anisotropic
    !! Symmetric anisotropic single-species Fickian diffusion in 3D.
    !!
    !! Spatial molar concentration \(c\) has dimensions of amount per volume.
    !! Its prescribed gradient \(\boldsymbol{g}=\nabla c\), in concentration
    !! per length, drives the diffusive molar flux \(\boldsymbol{j}\), in
    !! amount per area per time:
    !! \[\boldsymbol{j}=-\boldsymbol{D}\boldsymbol{g},\qquad
    !! \boldsymbol{D}=\begin{bmatrix}
    !! D_{11}&D_{12}&D_{13}\\[2pt]
    !! D_{12}&D_{22}&D_{23}\\[2pt]
    !! D_{13}&D_{23}&D_{33}
    !! \end{bmatrix}.\]
    !! The six independent tensor entries \((11,22,33,12,13,23)\) all have
    !! units of length squared per time. The tensor must be finite, symmetric and strictly
    !! positive definite in the spatial frame of the supplied fields.
    !! The down-gradient transport property holds for every input:
    !! \[-\boldsymbol{g}\cdot\boldsymbol{j}
    !! =\boldsymbol{g}^{\mathsf T}\boldsymbol{D}\boldsymbol{g}\geq0.\]
    !! This algebraic property is not equated with entropy production
    !! without an appropriate chemical-potential constitutive relation.
    !! Exact local derivatives at fixed prescribed \(\boldsymbol{g}\):
    !! \[\frac{\partial\boldsymbol{j}}{\partial\boldsymbol{g}}
    !! =-\boldsymbol{D},\qquad
    !! \frac{\partial j_i}{\partial D_{aa}}=-\delta_{ia}g_a,\qquad
    !! \frac{\partial j_i}{\partial D_{ab}}=
    !! -\delta_{ia}g_b-\delta_{ib}g_a\quad(a\ne b).\]
    !! Off-diagonal parameter perturbations affect two symmetric entries.
    !! Selected sensitivities are reported in caller selection order.
    !!
    !! Input `concentration.gradient` (concentration per length) and output
    !! `species.molar_flux` (amount per area per time) are spatial 3-vectors with
    !! Cartesian order \([1,2,3]\). When axes are rotated, the caller rotates
    !! both gradient and tensor. No unit or frame conversions occur.
    !! The diffusive flux is relative to the carrier or medium and excludes
    !! the advective contribution \(c\boldsymbol{v}\), where \(\boldsymbol{v}\) is the carrier velocity. Positive components denote diffusion
    !! toward positive spatial axes.
    !! Scaled Cholesky validation tests positive definiteness, rejecting
    !! tensors too near singularity for the chosen floating-point precision.
    !! This stateless, isothermal, dilute, nonreactive law has no time-step
    !! dependence, convection, nonlinear/cross diffusion, electromigration,
    !! multicomponent constraints or mass-balance equation solver.
    !! Parameter derivatives are local, not global solution sensitivities.
    !! Point loops have no heap allocations, strings or I/O. Nonfinite
    !! input and results invalidate the trial.
    !!
    !! Original law: A. Fick, *Ueber Diffusion*, Annalen der Physik
    !! 170(1), 59-86 (1855), DOI 10.1002/andp.18551700105.
    !! Modern molar-flux convention: IUPAC Gold Book, "diffusion coefficient",
    !! D01719, DOI 10.1351/goldbook.D01719:
    !! https://goldbook.iupac.org/terms/view/D01719 .
    !! Continuum context: J. Crank, *The Mathematics of Diffusion*,
    !! 2nd ed., Clarendon Press, Oxford (1975), ISBN 0-19-853344-6.
    !! Section 1.4, pp. 5--7, equation (1.9) gives the coupled anisotropic
    !! flux components. Section 1.2, p. 2, equation (1.1) discusses mass
    !! and molar concentration conventions with the same diffusivity units.
    !! The 3D symmetric SPD tensor restriction is a ConstituF adaptation,
    !! not attributed to the original scalar Fick law. Code provenance:
    !! this tensor specialization and its analytical derivatives directly
    !! implement the stated equations. No third-party software code is
    !! reused in this revision.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: ndim = 3
    integer, parameter :: nparameters = 6

    public :: fick_anisotropic_t

    type, extends(material_model_t) :: fick_anisotropic_t
        !! Spatial diffusivity tensor with six independent coefficients in consistent units.
        !!
        !! Construct as `fick_anisotropic_t(D11,D22,D33,D12,D13,D23)`.
        !! `material_t%init` clones this definition. Invalid or non-SPD
        !! tensors are rejected at initialization, not silently regularized.
        real(wp) :: diffusivity_11 = 0.0_wp
            !! Diagonal coefficient \(D_{11}\), length squared per time.
        real(wp) :: diffusivity_22 = 0.0_wp
            !! Diagonal coefficient \(D_{22}\), length squared per time.
        real(wp) :: diffusivity_33 = 0.0_wp
            !! Diagonal coefficient \(D_{33}\), length squared per time.
        real(wp) :: diffusivity_12 = 0.0_wp
            !! Symmetric coefficient \(D_{12}=D_{21}\), length squared per time.
        real(wp) :: diffusivity_13 = 0.0_wp
            !! Symmetric coefficient \(D_{13}=D_{31}\), length squared per time.
        real(wp) :: diffusivity_23 = 0.0_wp
            !! Symmetric coefficient \(D_{23}=D_{32}\), length squared per time.
    contains
        procedure :: describe => describe_fick
            !! Declare the constitutive fields, parameter metadata and supported capabilities.
        procedure :: eval_batch => evaluate_fick
            !! Evaluate the batched constitutive response, optional tangent and trial history.
        procedure :: eval_sensitivity_batch => sensitivity_fick
            !! Evaluate exact parameter-path derivatives and corresponding trial history.
    end type fick_anisotropic_t

contains

    subroutine describe_fick(self, spec, status)
        !! Declare the spatial gradient, flux and six scalar design parameters.
        !!
        !! An admissible tensor has a strictly positive Cholesky factorization.
        !! Supports exact \(3\times3\) tangent and optional local parameter derivatives.
        class(fick_anisotropic_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema receiving field layouts, parameters and supported operations.
        integer, intent(out) :: status
            !! Operation status. `cf_ok` indicates a valid model or successful initialization.
        real(wp) :: diffusivity(ndim,ndim)

        call build_tensor(self, diffusivity, status)
        if (status /= cf_ok) return

        call spec%add_input("concentration.gradient", ndim, "vector", "concentration_gradient", "spatial", [ndim])
        call spec%add_output("species.molar_flux", ndim, "vector", "molar_flux", "spatial", [ndim])
        call spec%add_parameter("diffusivity_11")
        call spec%add_parameter("diffusivity_22")
        call spec%add_parameter("diffusivity_33")
        call spec%add_parameter("diffusivity_12")
        call spec%add_parameter("diffusivity_13")
        call spec%add_parameter("diffusivity_23")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
    end subroutine describe_fick

    subroutine evaluate_fick(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Evaluate the flux and optional input Jacobian for a point batch.
        !!
        !! `x_new` and response have shape `(3,npoints)`. Optional tangent has
        !! shape `(3,3,npoints)`. Empty state remains empty. Failed trials
        !! cannot be committed and output contents on failure are unspecified.
        class(fick_anisotropic_t), intent(in) :: self
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
        !! Evaluate selected local analytical flux sensitivities.
        !!
        !! `dsdp(3,nselected,npoints)` follows the exact setup selection order.
        !! The gradient is held fixed. No parameter-dependent prescribed
        !! histories or global solution derivatives are included. Neither
        !! ordinary nor sensitivity evaluation retains history or allocates.
        class(fick_anisotropic_t), intent(in) :: self
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

    subroutine evaluate_points(self, gradient, response, point_status, tangent, selected_parameters, dsdp)
        !! Build diffusivity once per batch, with no per-point allocations.
        class(fick_anisotropic_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        real(wp), intent(in) :: gradient(:,:)
            !! Prescribed `concentration.gradient` gradient components at each point.
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
        real(wp) :: diffusivity(ndim,ndim)
        integer :: ipoint
        integer :: iparam
        integer :: id
        integer :: status
        integer :: i

        point_status = cf_invalid_argument
        call build_tensor(self, diffusivity, status)
        if (status /= cf_ok) return
        if (present(selected_parameters)) then
            if (any(selected_parameters < 1) .or. any(selected_parameters > nparameters)) return
        end if

        point_loop: do ipoint = 1, size(gradient,2)
            if (.not. all(ieee_is_finite(gradient(:,ipoint)))) cycle point_loop
            do i = 1, ndim
                response(i,ipoint) = -dot_product(diffusivity(i,:),gradient(:,ipoint))
                if (.not. ieee_is_finite(response(i,ipoint))) then
                    response(i,ipoint) = scaled_flux_component(diffusivity(i,:),gradient(:,ipoint))
                end if
            end do
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop

            if (present(selected_parameters)) then
                do iparam = 1, size(selected_parameters)
                    id = selected_parameters(iparam)
                    dsdp(:,iparam,ipoint) = 0.0_wp
                    select case (id)
                    case (1:3)
                        dsdp(id,iparam,ipoint) = -gradient(id,ipoint)
                    case (4)
                        dsdp(1,iparam,ipoint) = -gradient(2,ipoint)
                        dsdp(2,iparam,ipoint) = -gradient(1,ipoint)
                    case (5)
                        dsdp(1,iparam,ipoint) = -gradient(3,ipoint)
                        dsdp(3,iparam,ipoint) = -gradient(1,ipoint)
                    case (6)
                        dsdp(2,iparam,ipoint) = -gradient(3,ipoint)
                        dsdp(3,iparam,ipoint) = -gradient(2,ipoint)
                    case default
                        cycle point_loop
                    end select
                end do
                if (.not. all(ieee_is_finite(dsdp(:,:,ipoint)))) cycle point_loop
            end if
            if (present(tangent)) tangent(:,:,ipoint) = -diffusivity
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_points

    pure function scaled_flux_component(diffusivity_row, gradient) result(flux)
        !! Retry an overflowing dot product using bounded product fractions.
        !! Sum larger terms first so cancellation preserves small finite terms.
        real(wp), intent(in) :: diffusivity_row(:)
            !! Finite coefficients of one spatial diffusivity row.
        real(wp), intent(in) :: gradient(ndim)
            !! Finite prescribed spatial molar concentration gradient.
        real(wp) :: flux
        real(wp) :: term_fraction(ndim)
        real(wp) :: sum_fraction
        integer :: term_exponent(ndim)
        integer :: sum_exponent
        integer :: common_exponent
        integer :: largest
        integer :: i
        integer :: j

        term_fraction = 0.0_wp
        term_exponent = 0
        product_loop: do i = 1, ndim
            if (abs(diffusivity_row(i)) <= 0.0_wp .or. abs(gradient(i)) <= 0.0_wp) cycle product_loop
            term_fraction(i) = fraction(diffusivity_row(i))*fraction(gradient(i))
            term_exponent(i) = exponent(diffusivity_row(i))+exponent(gradient(i))+exponent(term_fraction(i))
            term_fraction(i) = fraction(term_fraction(i))
        end do product_loop

        sum_fraction = 0.0_wp
        sum_exponent = 0
        sum_loop: do i = 1, ndim
            largest = 0
            largest_loop: do j = 1, ndim
                if (abs(term_fraction(j)) <= 0.0_wp) cycle largest_loop
                if (largest == 0) then
                    largest = j
                    cycle largest_loop
                end if
                if (term_exponent(j) > term_exponent(largest)) then
                    largest = j
                else if (term_exponent(j) == term_exponent(largest)) then
                    if (abs(term_fraction(j)) > abs(term_fraction(largest))) largest = j
                end if
            end do largest_loop
            if (largest == 0) exit sum_loop

            if (abs(sum_fraction) <= 0.0_wp) then
                sum_fraction = term_fraction(largest)
                sum_exponent = term_exponent(largest)
            else
                common_exponent = max(sum_exponent,term_exponent(largest))
                sum_fraction = scale(sum_fraction,sum_exponent-common_exponent)+scale(term_fraction(largest),term_exponent(largest)-common_exponent)
                sum_exponent = common_exponent+exponent(sum_fraction)
                sum_fraction = fraction(sum_fraction)
            end if
            term_fraction(largest) = 0.0_wp
        end do sum_loop
        flux = -scale(sum_fraction,sum_exponent)
    end function scaled_flux_component

    subroutine build_tensor(self, diffusivity, status)
        !! Form a symmetric tensor and validate scaled Cholesky pivots.
        !!
        !! Scaling prevents large finite coefficients from overflowing
        !! determinant products. Strictly positive computed pivots are
        !! required. No diagonal stabilization or clipping is performed.
        class(fick_anisotropic_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        real(wp), intent(out) :: diffusivity(ndim,ndim)
            !! Symmetric spatial diffusivity tensor assembled from the model coefficients.
        integer, intent(out) :: status
            !! Operation status. `cf_ok` indicates a valid model or successful initialization.
        real(wp) :: scaled(ndim,ndim)
        real(wp) :: scale
        real(wp) :: l11
        real(wp) :: l21
        real(wp) :: l31
        real(wp) :: l22
        real(wp) :: l32
        real(wp) :: pivot

        status = cf_invalid_argument
        diffusivity(1,:) = [self%diffusivity_11, self%diffusivity_12, self%diffusivity_13]
        diffusivity(2,:) = [self%diffusivity_12, self%diffusivity_22, self%diffusivity_23]
        diffusivity(3,:) = [self%diffusivity_13, self%diffusivity_23, self%diffusivity_33]
        if (.not. all(ieee_is_finite(diffusivity))) return
        scale = maxval(abs(diffusivity))
        if (scale <= 0.0_wp) return
        scaled = diffusivity/scale

        pivot = scaled(1,1)
        if (pivot <= 0.0_wp) return
        l11 = sqrt(pivot)
        l21 = scaled(2,1)/l11
        l31 = scaled(3,1)/l11
        pivot = scaled(2,2)-l21*l21
        if (.not. ieee_is_finite(pivot)) return
        if (pivot <= 0.0_wp) return
        l22 = sqrt(pivot)
        l32 = (scaled(3,2)-l31*l21)/l22
        pivot = scaled(3,3)-l31*l31-l32*l32
        if (.not. ieee_is_finite(pivot)) return
        if (pivot <= 0.0_wp) return
        status = cf_ok
    end subroutine build_tensor

end module constituf_fick_anisotropic

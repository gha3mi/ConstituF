! SPDX-License-Identifier: Apache-2.0

module constituf_ohmic_anisotropic
    !! Homogeneous symmetric anisotropic Ohmic conduction in 3D.
    !!
    !! The prescribed electric field \(\boldsymbol{E}\) (electric potential per length) produces
    !! conduction current density \(\boldsymbol{J}\) (current per area) through
    !! \[\boldsymbol{J}=\boldsymbol{\sigma}\boldsymbol{E},\qquad
    !! \boldsymbol{\sigma}=\begin{bmatrix}
    !! \sigma_{11}&\sigma_{12}&\sigma_{13}\\[2pt]
    !! \sigma_{12}&\sigma_{22}&\sigma_{23}\\[2pt]
    !! \sigma_{13}&\sigma_{23}&\sigma_{33}
    !! \end{bmatrix}.\]
    !! All entries are in current density per electric field. The six independently selectable scalars
    !! are ordered by tensor indices \((11,22,33,12,13,23)\). This model *requires* finite, real,
    !! strictly symmetric positive-definite conductivity, so the Joule
    !! dissipation density is positive for all nonzero physical fields:
    !! \[\boldsymbol{E}^{\mathsf T}\boldsymbol{\sigma}
    !! \boldsymbol{E}>0,\qquad
    !! \boldsymbol{E}\ne\boldsymbol{0}.\]
    !! Symmetry is an explicit model restriction. Nonnegative Joule
    !! dissipation alone constrains only the symmetric part of a tensor.
    !! Exact local derivatives at fixed \(\boldsymbol{E}\) are
    !! \[\frac{\partial\boldsymbol{J}}{\partial\boldsymbol{E}}
    !! =\boldsymbol{\sigma},\quad
    !! \frac{\partial J_i}{\partial\sigma_{aa}}=\delta_{ia}E_a,\quad
    !! \frac{\partial J_i}{\partial\sigma_{ab}}=
    !! \delta_{ia}E_b+\delta_{ib}E_a\ (a\ne b).\]
    !! Both fields and tensor use matching spatial Cartesian axes \([1,2,3]\).
    !! \(\boldsymbol{E}\) is the physical electric field. \(\boldsymbol{J}\) is conventional conduction
    !! current density, with no potential-gradient sign change applied.
    !! The medium is stationary in the supplied spatial frame. Motional
    !! electromotive terms and convective charge transport are excluded.
    !! The caller must rotate field and tensor together when changing axes.
    !! No automatic unit/frame transforms or global solution derivatives.
    !! Scaled Cholesky checks strict SPD. Almost singular tensors can fail.
    !! No time-step dependence, internal state, magnetic Hall current,
    !! displacement current, thermoelectric effects, temperature coupling,
    !! nonlinear conduction or electric current conservation solve.
    !! The point loop is allocation-free and rejects nonfinite responses.
    !!
    !! Scientific basis: G. S. Ohm, *Die galvanische Kette, mathematisch
    !! bearbeitet*, T. H. Riemann, Berlin (1827). Digitized edition DOI
    !! 10.5479/sil.354716.39088005838644. The standard tensor Ohm law
    !! is used here with the chosen symmetric-SPD conductivity restriction.
    !! Direct tensor reference: COMSOL AC/DC Module User's Guide 6.4,
    !! "Current Conservation", subsection "Constitutive Relation Jc-E":
    !! https://doc.comsol.com/6.4/doc/com.comsol.help.acdc/acdc_ug_electric_fields.07.050.html .
    !! Code provenance: this tensor law
    !! and its analytical derivatives directly implement the stated
    !! equations. No third-party software code is reused in this revision.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: ndim = 3
    integer, parameter :: nparameters = 6

    public :: ohmic_anisotropic_t

    type, extends(material_model_t) :: ohmic_anisotropic_t
        !! Create with `ohmic_anisotropic_t(sigma11,sigma22,sigma33,sigma12,sigma13,sigma23)`.
        !! The coefficients are in current density per electric field. `material_t%init` clones the model.
        real(wp) :: conductivity_11 = 0.0_wp
            !! Spatial tensor entry \(\sigma_{11}\), current density per electric field.
        real(wp) :: conductivity_22 = 0.0_wp
            !! Spatial tensor entry \(\sigma_{22}\), current density per electric field.
        real(wp) :: conductivity_33 = 0.0_wp
            !! Spatial tensor entry \(\sigma_{33}\), current density per electric field.
        real(wp) :: conductivity_12 = 0.0_wp
            !! Symmetric entry \(\sigma_{12}=\sigma_{21}\), current density per electric field.
        real(wp) :: conductivity_13 = 0.0_wp
            !! Symmetric entry \(\sigma_{13}=\sigma_{31}\), current density per electric field.
        real(wp) :: conductivity_23 = 0.0_wp
            !! Symmetric entry \(\sigma_{23}=\sigma_{32}\), current density per electric field.
    contains
        procedure :: describe => describe_ohmic_anisotropic
            !! Declare the constitutive fields, parameter metadata and supported capabilities.
        procedure :: eval_batch => evaluate_ohmic_anisotropic
            !! Evaluate the batched constitutive response, optional tangent and trial history.
        procedure :: eval_sensitivity_batch => sensitivity_ohmic_anisotropic
            !! Evaluate exact parameter-path derivatives and corresponding trial history.
    end type ohmic_anisotropic_t

contains

    subroutine describe_ohmic_anisotropic(self, spec, status)
        !! Define tensor material with six selectable scalar sensitivities.
        class(ohmic_anisotropic_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema receiving field layouts, parameters and supported operations.
        integer, intent(out) :: status
            !! Operation status. `cf_ok` indicates a valid model or successful initialization.
        real(wp) :: conductivity(ndim,ndim)

        call build_tensor(self,conductivity,status)
        if (status /= cf_ok) return

        call spec%add_input("electric.field", ndim, "vector", "electric_field", "spatial", [ndim])
        call spec%add_output("electric.current_density", ndim, "vector", "electric_current_density", "spatial", [ndim])
        call spec%add_parameter("conductivity_11")
        call spec%add_parameter("conductivity_22")
        call spec%add_parameter("conductivity_33")
        call spec%add_parameter("conductivity_12")
        call spec%add_parameter("conductivity_13")
        call spec%add_parameter("conductivity_23")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_ohmic_anisotropic

    subroutine evaluate_ohmic_anisotropic(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Evaluate \(\boldsymbol{J}=\boldsymbol{\sigma}\boldsymbol{E}\) over a caller-owned batch, without state.
        class(ohmic_anisotropic_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed input `(ninput,npoints)` if required, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `electric.field` input `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued state `(nstate,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Trial real-valued state `(nstate,npoints)` returned by this evaluation.
        real(wp), intent(out) :: response(:,:)
            !! Packed `electric.current_density` response `(noutput,npoints)` for the evaluated material points.
        real(wp), intent(in) :: dt
            !! Prescribed time increment. It is positive when the material requires time integration.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact discrete response Jacobian `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(inout) :: point_status(:)
            !! Per-point status values; successful evaluations must set each entry to `cf_ok`.

        state_trial = state_old
        call evaluate_points(self,x_new,response,point_status,tangent)
    end subroutine evaluate_ohmic_anisotropic

    subroutine sensitivity_ohmic_anisotropic(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Exact local derivatives at fixed field in requested parameter order.
        !! `dsdp` has shape `(3,nselected,npoints)`. History buffers are empty.
        class(ohmic_anisotropic_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed input `(ninput,npoints)` if required, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `electric.field` input `(ninput,npoints)`.
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
            !! Packed `electric.current_density` response `(noutput,npoints)` for the evaluated material points.
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
        call evaluate_points(self,x_new,response,point_status,tangent, selected_parameters,dsdp)
    end subroutine sensitivity_ohmic_anisotropic

    subroutine evaluate_points(self, electric, response, point_status, tangent, selected_parameters, dsdp)
        !! Validate tensor once per batch. No point-local allocations or I/O.
        class(ohmic_anisotropic_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        real(wp), intent(in) :: electric(:,:)
            !! Prescribed electric field `electric.field` at each point.
        real(wp), intent(out) :: response(:,:)
            !! Packed `electric.current_density` response `(noutput,npoints)` for the evaluated material points.
        integer, intent(inout) :: point_status(:)
            !! Per-point status values; successful evaluations must set each entry to `cf_ok`.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact discrete response Jacobian `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(in), optional :: selected_parameters(:)
            !! One-based IDs of selected model parameters, in sensitivity-axis order.
        real(wp), intent(out), optional :: dsdp(:,:,:)
            !! Exact total response derivatives `(noutput,nselected,npoints)` at fixed prescribed input history.
        real(wp) :: conductivity(ndim,ndim)
        integer :: ipoint
        integer :: iparam
        integer :: id
        integer :: i
        integer :: status

        point_status = cf_invalid_argument
        call build_tensor(self,conductivity,status)
        if (status /= cf_ok) return
        if (present(selected_parameters)) then
            if (any(selected_parameters < 1) .or. any(selected_parameters > nparameters)) return
        end if

        point_loop: do ipoint = 1, size(electric,2)
            if (.not. all(ieee_is_finite(electric(:,ipoint)))) cycle point_loop
            do i = 1, ndim
                response(i,ipoint) = dot_product(conductivity(i,:),electric(:,ipoint))
                if (.not. ieee_is_finite(response(i,ipoint))) then
                    response(i,ipoint) = scaled_current_component(conductivity(i,:),electric(:,ipoint))
                end if
            end do
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop

            if (present(selected_parameters)) then
                do iparam = 1, size(selected_parameters)
                    id = selected_parameters(iparam)
                    dsdp(:,iparam,ipoint) = 0.0_wp
                    select case (id)
                    case (1:3)
                        dsdp(id,iparam,ipoint) = electric(id,ipoint)
                    case (4)
                        dsdp(1,iparam,ipoint) = electric(2,ipoint)
                        dsdp(2,iparam,ipoint) = electric(1,ipoint)
                    case (5)
                        dsdp(1,iparam,ipoint) = electric(3,ipoint)
                        dsdp(3,iparam,ipoint) = electric(1,ipoint)
                    case (6)
                        dsdp(2,iparam,ipoint) = electric(3,ipoint)
                        dsdp(3,iparam,ipoint) = electric(2,ipoint)
                    case default
                        cycle point_loop
                    end select
                end do
                if (.not. all(ieee_is_finite(dsdp(:,:,ipoint)))) cycle point_loop
            end if
            if (present(tangent)) tangent(:,:,ipoint) = conductivity
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_points

    pure function scaled_current_component(conductivity_row, electric) result(current_density)
        !! Retry an overflowing dot product using bounded product fractions.
        !! Sum larger terms first so cancellation preserves small finite terms.
        real(wp), intent(in) :: conductivity_row(:)
            !! Finite coefficients of one spatial electrical-conductivity row.
        real(wp), intent(in) :: electric(ndim)
            !! Finite prescribed spatial electric field.
        real(wp) :: current_density
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
            if (abs(conductivity_row(i)) <= 0.0_wp .or. abs(electric(i)) <= 0.0_wp) cycle product_loop
            term_fraction(i) = fraction(conductivity_row(i))*fraction(electric(i))
            term_exponent(i) = exponent(conductivity_row(i))+exponent(electric(i))+exponent(term_fraction(i))
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
        current_density = scale(sum_fraction,sum_exponent)
    end function scaled_current_component

    subroutine build_tensor(self, conductivity, status)
        !! Check SPD by Cholesky pivots after positive coefficient scaling.
        !! No numerical clipping or regularization is applied.
        class(ohmic_anisotropic_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        real(wp), intent(out) :: conductivity(ndim,ndim)
            !! Symmetric spatial electrical conductivity tensor assembled from the model coefficients.
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
        conductivity(1,:) = [self%conductivity_11,self%conductivity_12, self%conductivity_13]
        conductivity(2,:) = [self%conductivity_12,self%conductivity_22, self%conductivity_23]
        conductivity(3,:) = [self%conductivity_13,self%conductivity_23, self%conductivity_33]
        if (.not. all(ieee_is_finite(conductivity))) return
        scale = maxval(abs(conductivity))
        if (scale <= 0.0_wp) return
        scaled = conductivity/scale

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

end module constituf_ohmic_anisotropic

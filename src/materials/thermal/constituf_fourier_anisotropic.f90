! SPDX-License-Identifier: Apache-2.0

module constituf_fourier_anisotropic
    !! Homogeneous, symmetric anisotropic Fourier conduction in three dimensions.
    !!
    !! The prescribed spatial temperature gradient \(\boldsymbol{g}=\nabla T\)
    !! produces the spatial heat-flux vector \(\boldsymbol{q}\):
    !! \[\boldsymbol{q}=-\boldsymbol{K}\boldsymbol{g},\qquad
    !! \boldsymbol{K}=
    !! \begin{bmatrix}
    !! k_{11}&k_{12}&k_{13}\\[2pt]
    !! k_{12}&k_{22}&k_{23}\\[2pt]
    !! k_{13}&k_{23}&k_{33}
    !! \end{bmatrix}.\]
    !! All six independent coefficients have units of power per length per temperature interval and
    !! \(\boldsymbol{K}\) must be finite, symmetric and strictly positive
    !! definite. The exact local Jacobian at fixed conductivity is
    !! \[\frac{\partial\boldsymbol{q}}{\partial\boldsymbol{g}}
    !! =-\boldsymbol{K}.\]
    !! For \(p=k_{ab}\), the analytical local sensitivity at fixed
    !! prescribed gradient is \(\partial\boldsymbol{q}/\partial p
    !! =-(\partial\boldsymbol{K}/\partial p)\boldsymbol{g}\).
    !! An off-diagonal coefficient changes both symmetric matrix entries.
    !! The down-gradient property is
    !! \[-\boldsymbol{q}\cdot\nabla T
    !! =\boldsymbol{g}^{\mathsf T}\boldsymbol{K}\boldsymbol{g}\geq0.\]
    !! At absolute \(T>0\), conductive entropy production is this quantity
    !! divided by \(T^2\), not the quadratic form alone.
    !!
    !! Input `temperature.gradient` is a spatial vector in temperature interval per length.
    !! Output `heat.flux` is a spatial vector in power-per-area units. The six parameters
    !! refer to fixed axes aligned with the solver's spatial coordinates.
    !! Positive components denote heat transport toward positive coordinate
    !! directions. The outward flux through a surface with outward unit
    !! normal \(\boldsymbol{n}\) is \(\boldsymbol{q}\cdot\boldsymbol{n}\).
    !! The caller must rotate \(\boldsymbol{K}\) when changing the frame.
    !! The library does not perform material-axis or unit transformations.
    !! There is no history, time-step dependence, heat capacity, temperature
    !! dependence or heat-equation solve. Near-singular tensors can be
    !! rejected when the Cholesky pivots lose numerical significance.
    !!
    !! Original law: J. Fourier, *Théorie analytique de la chaleur*,
    !! Paris, 1822. Anisotropic tensor treatment: H. S. Carslaw and
    !! J. C. Jaeger, *Conduction of Heat in Solids*, 2nd ed.,
    !! Clarendon Press, Oxford, 1959, ISBN 978-0-19-853368-9.
    !! Direct tensor law and Onsager symmetry: G. P. Beretta, MIT 2.43
    !! Advanced Thermodynamics (2024), Lecture 21, slides 21.17--21.18:
    !! https://ocw.mit.edu/courses/2-43-advanced-thermodynamics-spring-2024/mit2_43_s24_lec21.pdf .
    !! Code provenance: the tensor flux and its analytical derivatives
    !! directly implement the stated law. No third-party software code is
    !! reused in this revision.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: ndim = 3
    integer, parameter :: nparameters = 6

    public :: fourier_anisotropic_t

    type, extends(material_model_t) :: fourier_anisotropic_t
        !! Spatial conductivity tensor with six independent coefficients in consistent units.
        !!
        !! Construct as `fourier_anisotropic_t(k11,k22,k33,k12,k13,k23)`.
        !! `material_t%init` clones this definition. Invalid or non-SPD
        !! tensors are rejected at initialization, not silently regularized.
        real(wp) :: conductivity_11 = 0.0_wp
            !! Diagonal coefficient \(k_{11}\), power per length per temperature interval.
        real(wp) :: conductivity_22 = 0.0_wp
            !! Diagonal coefficient \(k_{22}\), power per length per temperature interval.
        real(wp) :: conductivity_33 = 0.0_wp
            !! Diagonal coefficient \(k_{33}\), power per length per temperature interval.
        real(wp) :: conductivity_12 = 0.0_wp
            !! Symmetric coefficient \(k_{12}=k_{21}\), power per length per temperature interval.
        real(wp) :: conductivity_13 = 0.0_wp
            !! Symmetric coefficient \(k_{13}=k_{31}\), power per length per temperature interval.
        real(wp) :: conductivity_23 = 0.0_wp
            !! Symmetric coefficient \(k_{23}=k_{32}\), power per length per temperature interval.
    contains
        procedure :: describe => describe_fourier
            !! Declare material fields, parameters and supported operations.
        procedure :: eval_batch => evaluate_fourier
            !! Evaluate the batched material response, optional tangent and trial state.
        procedure :: eval_sensitivity_batch => sensitivity_fourier
            !! Evaluate analytical parameter derivatives through the discrete material update.
    end type fourier_anisotropic_t

contains

    subroutine describe_fourier(self, spec, status)
        !! Declare the spatial gradient, flux and six scalar design parameters.
        !!
        !! An admissible tensor has a strictly positive Cholesky factorization.
        !! Supports exact \(3\times3\) tangent and optional local parameter derivatives.
        class(fourier_anisotropic_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema populated with fields, physical metadata and model capabilities.
        integer, intent(out) :: status
            !! Operation status code. `cf_ok` indicates success.
        real(wp) :: conductivity(ndim,ndim)

        call build_tensor(self, conductivity, status)
        if (status /= cf_ok) return

        call spec%add_input("temperature.gradient", ndim, "vector", "temperature_gradient", "spatial", [ndim])
        call spec%add_output("heat.flux", ndim, "vector", "heat_flux", "spatial", [ndim])
        call spec%add_parameter("conductivity_11")
        call spec%add_parameter("conductivity_22")
        call spec%add_parameter("conductivity_33")
        call spec%add_parameter("conductivity_12")
        call spec%add_parameter("conductivity_13")
        call spec%add_parameter("conductivity_23")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
    end subroutine describe_fourier

    subroutine evaluate_fourier(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Evaluate the flux and optional input Jacobian for a point batch.
        !!
        !! `x_new` and response have shape `(3,npoints)`. Optional tangent has
        !! shape `(3,3,npoints)`. Empty state remains empty. Failed trials
        !! cannot be committed and output contents on failure are unspecified.
        class(fourier_anisotropic_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previous packed inputs `(ninput,npoints)` when required by the model, otherwise a zero-sized array.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `temperature.gradient` inputs `(ninput,npoints)` at the current time step.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued history `(nstate,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Updated trial real-valued history `(nstate,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `heat.flux` outputs `(noutput,npoints)`.
        real(wp), intent(in) :: dt
            !! Prescribed time increment, strictly positive for time-dependent materials.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact discrete input Jacobian `(noutput,ninput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Per-point status codes. successful points must be assigned `cf_ok`.

        state_trial = state_old
        call evaluate_points(self, x_new, response, point_status, tangent)
    end subroutine evaluate_fourier

    subroutine sensitivity_fourier(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Evaluate selected local analytical flux sensitivities.
        !!
        !! `dsdp(3,nselected,npoints)` follows the exact setup selection order.
        !! The gradient is held fixed. No parameter-dependent prescribed
        !! histories or global solution derivatives are included. Neither
        !! ordinary nor sensitivity evaluation retains history or allocates.
        class(fourier_anisotropic_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previous packed inputs `(ninput,npoints)` when required by the model, otherwise a zero-sized array.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `temperature.gradient` inputs `(ninput,npoints)` at the current time step.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued history `(nstate,npoints)`.
        integer, intent(in) :: selected_parameters(:)
            !! One-based IDs of scalar parameters in the requested derivative-axis order.
        real(wp), intent(in) :: sensitivity_old(:,:,:)
            !! Committed history derivatives `(nstate,nselected,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Updated trial real-valued history `(nstate,npoints)`.
        real(wp), intent(out) :: sensitivity_trial(:,:,:)
            !! Updated trial history derivatives `(nstate,nselected,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `heat.flux` outputs `(noutput,npoints)`.
        real(wp), intent(out) :: dsdp(:,:,:)
            !! Exact response derivatives `(noutput,nselected,npoints)` at fixed prescribed input histories.
        real(wp), intent(in) :: dt
            !! Prescribed time increment, strictly positive for time-dependent materials.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact discrete input Jacobian `(noutput,ninput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Per-point status codes. successful points must be assigned `cf_ok`.

        state_trial = state_old
        sensitivity_trial = sensitivity_old
        call evaluate_points(self, x_new, response, point_status, tangent, selected_parameters, dsdp)
    end subroutine sensitivity_fourier

    subroutine evaluate_points(self, gradient, response, point_status, tangent, selected_parameters, dsdp)
        !! Build conductivity once per batch, with no per-point allocations.
        class(fourier_anisotropic_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        real(wp), intent(in) :: gradient(:,:)
            !! Spatial `temperature.gradient` values at each material point.
        real(wp), intent(out) :: response(:,:)
            !! Packed `heat.flux` outputs `(noutput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Per-point status codes. successful points must be assigned `cf_ok`.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact discrete input Jacobian `(noutput,ninput,npoints)`.
        integer, intent(in), optional :: selected_parameters(:)
            !! One-based IDs of scalar parameters in the requested derivative-axis order.
        real(wp), intent(out), optional :: dsdp(:,:,:)
            !! Exact response derivatives `(noutput,nselected,npoints)` at fixed prescribed input histories.
        real(wp) :: conductivity(ndim,ndim)
        integer :: icomponent
        integer :: ipoint
        integer :: iparam
        integer :: id
        integer :: status

        point_status = cf_invalid_argument
        call build_tensor(self, conductivity, status)
        if (status /= cf_ok) return
        if (present(selected_parameters)) then
            if (any(selected_parameters < 1) .or. any(selected_parameters > nparameters)) return
        end if

        point_loop: do ipoint = 1, size(gradient,2)
            if (.not. all(ieee_is_finite(gradient(:,ipoint)))) cycle point_loop
            response(:,ipoint) = -matmul(conductivity, gradient(:,ipoint))
            if (.not. all(ieee_is_finite(response(:,ipoint)))) then
                do icomponent = 1, ndim
                    if (.not. ieee_is_finite(response(icomponent,ipoint))) then
                        response(icomponent,ipoint) = scaled_flux_component(conductivity(icomponent,:), gradient(:,ipoint))
                    end if
                end do
            end if
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
            if (present(tangent)) tangent(:,:,ipoint) = -conductivity
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_points

    pure function scaled_flux_component(conductivity_row, gradient) result(flux)
        !! Retry an overflowing dot product using bounded product fractions.
        !! Sum larger terms first so cancellation preserves small finite terms.
        real(wp), intent(in) :: conductivity_row(:)
            !! Finite coefficients of one spatial conductivity row.
        real(wp), intent(in) :: gradient(ndim)
            !! Finite prescribed spatial temperature gradient.
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
            if (abs(conductivity_row(i)) <= 0.0_wp .or. abs(gradient(i)) <= 0.0_wp) cycle product_loop
            term_fraction(i) = fraction(conductivity_row(i))*fraction(gradient(i))
            term_exponent(i) = exponent(conductivity_row(i))+exponent(gradient(i))+exponent(term_fraction(i))
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

    subroutine build_tensor(self, conductivity, status)
        !! Form a symmetric tensor and validate scaled Cholesky pivots.
        !!
        !! Scaling prevents large finite coefficients from overflowing
        !! determinant products. Strictly positive computed pivots are
        !! required. No diagonal stabilization or clipping is performed.
        class(fourier_anisotropic_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        real(wp), intent(out) :: conductivity(ndim,ndim)
            !! Symmetric three-by-three conductivity tensor for Fourier heat conduction.
        integer, intent(out) :: status
            !! Operation status code. `cf_ok` indicates success.
        real(wp) :: scaled(ndim,ndim)
        real(wp) :: scale
        real(wp) :: l11
        real(wp) :: l21
        real(wp) :: l31
        real(wp) :: l22
        real(wp) :: l32
        real(wp) :: pivot

        status = cf_invalid_argument
        conductivity(1,:) = [self%conductivity_11, self%conductivity_12, self%conductivity_13]
        conductivity(2,:) = [self%conductivity_12, self%conductivity_22, self%conductivity_23]
        conductivity(3,:) = [self%conductivity_13, self%conductivity_23, self%conductivity_33]
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

end module constituf_fourier_anisotropic

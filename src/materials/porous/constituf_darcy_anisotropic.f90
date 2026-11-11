! SPDX-License-Identifier: Apache-2.0

module constituf_darcy_anisotropic
    !! Three-dimensional anisotropic saturated single-phase Darcy flow.
    !!
    !! A fixed symmetric intrinsic permeability tensor \(\boldsymbol{K}\)
    !! (length squared) and scalar dynamic viscosity \(\mu>0\) (stress times time) define
    !! \[\boldsymbol{v}=-\mu^{-1}\boldsymbol{K}\boldsymbol{g},\quad
    !! \boldsymbol{g}=\nabla p.\]
    !! The prescribed spatial pressure gradient \(\boldsymbol{g}\) has
    !! units of pressure per length and the superficial (Darcy) velocity \(\boldsymbol{v}\)
    !! has units of length per time. This is specific discharge relative to a stationary
    !! porous skeleton, not interstitial pore velocity. The independent
    !! permeability coefficients are
    !! \(K_{11},K_{22},K_{33},K_{12},K_{13},K_{23}\), all in squared-length units. Off-diagonal
    !! entries are symmetric, \(K_{ij}=K_{ji}\), with no Mandel factors.
    !! The tensor must be finite,
    !! strictly positive definite in the prescribed spatial axes. The
    !! dynamic viscosity must be finite and positive. The mobility tensor
    !! \(\boldsymbol{K}/\mu\) must have finite entries and at least one nonzero entry
    !! in the working precision. All parameters are constant in space,
    !! pressure and time for this material.
    !!
    !! The exact input Jacobian is
    !! \[\frac{\partial\boldsymbol{v}}{\partial\boldsymbol{g}}
    !! =-\frac{\boldsymbol{K}}{\mu}.\]
    !! For a selected scalar parameter \(K_{ab}\), the local derivative is
    !! \[\frac{\partial v_i}{\partial K_{aa}}=
    !! -\frac{\delta_{ia}g_a}{\mu},\qquad
    !! \frac{\partial v_i}{\partial K_{ab}}=
    !! -\frac{\delta_{ia}g_b+\delta_{ib}g_a}{\mu}\quad(a\ne b).\]
    !! The viscosity sensitivity at fixed prescribed gradient is
    !! \[\frac{\partial\boldsymbol{v}}{\partial\mu}
    !! =\frac{\boldsymbol{K}\boldsymbol{g}}{\mu^2}
    !! =-\frac{\boldsymbol{v}}{\mu}.\]
    !! The pressure-work dissipation rate per unit bulk volume (power per volume)
    !! is nonnegative by tensor positive definiteness:
    !! \[-\boldsymbol{g}\cdot\boldsymbol{v}
    !! =\mu^{-1}\boldsymbol{g}^{T}\boldsymbol{K}\boldsymbol{g}\geq0.\]
    !! This derived rate is not a returned response field.
    !!
    !! Input `pressure.gradient` is a spatial 3-vector in pressure per length. Output
    !! `fluid.darcy_velocity` is a spatial 3-vector in length per time. The field
    !! ordering is Cartesian \([1,2,3]\) and no axis rotation or unit conversion
    !! is performed. Coefficients are defined in the caller's spatial
    !! frame. This is stateless, single-phase and isothermal with no `dt`
    !! dependence. Gravity is omitted for the stated pressure-gradient
    !! input. A caller including a body acceleration \(\boldsymbol b\)
    !! may supply \(\boldsymbol g_{\mathrm{eff}}=\nabla p-\rho\boldsymbol b\)
    !! externally in the same field and units. Both the flux law and
    !! dissipation then use this effective driving gradient. No porosity,
    !! saturation changes, multiphase effects, deformation, temperature
    !! dependence, inertial resistance or global flow-equation solve is
    !! included. Derivatives are local with fixed prescribed fields,
    !! not global solver sensitivities. Cholesky
    !! validation may reject near-singular tensors owing to roundoff.
    !! Failures leave committed history untouched. The point loop
    !! performs no heap allocation, file I/O or string matching.
    !!
    !! Scientific basis (not code provenance): H. Darcy, *Les fontaines
    !! publiques de la ville de Dijon*, Dalmont, Paris (1856), Appendix D.
    !! J. Bear, *Dynamics of Fluids in Porous Media*, American Elsevier,
    !! New York (1972).
    !! The standard tensor law is restricted here to symmetric positive
    !! definite permeability. The modern pressure form and specific-discharge
    !! convention appear in the COMSOL Porous Media Flow Module User's
    !! Guide 6.1 (2022), "Darcy's Law--Equation Formulation", Eqs. (4-1)
    !! and (4-3), pp. 114-116. Symmetric tensor permeability is described
    !! in "Porous Matrix--Darcian Flow", p. 126:
    !! https://doc.comsol.com/6.1/doc/com.comsol.help.porous/PorousMediaFlowModuleUsersGuide.pdf .
    !! No third-party software code is reused in this revision.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: ndim = 3
    integer, parameter :: nparameters = 7

    public :: darcy_anisotropic_t

    type, extends(material_model_t) :: darcy_anisotropic_t
        !! Construct as `darcy_anisotropic_t(k11,k22,k33,k12,k13,k23,mu)`.
        !! Scalar permeability components have units of length squared, viscosity stress times time.
        real(wp) :: permeability_11 = 0.0_wp
            !! Diagonal \(K_{11}\), in squared-length units.
        real(wp) :: permeability_22 = 0.0_wp
            !! Diagonal \(K_{22}\), in squared-length units.
        real(wp) :: permeability_33 = 0.0_wp
            !! Diagonal \(K_{33}\), in squared-length units.
        real(wp) :: permeability_12 = 0.0_wp
            !! Symmetric off-diagonal \(K_{12}=K_{21}\), in squared-length units.
        real(wp) :: permeability_13 = 0.0_wp
            !! Symmetric off-diagonal \(K_{13}=K_{31}\), in squared-length units.
        real(wp) :: permeability_23 = 0.0_wp
            !! Symmetric off-diagonal \(K_{23}=K_{32}\), in squared-length units.
        real(wp) :: dynamic_viscosity = 0.0_wp
            !! Positive dynamic viscosity \(\mu\), stress times time.
    contains
        procedure :: describe => describe_darcy_anisotropic
            !! Describe the material input, response, state and scalar parameter schema.
        procedure :: eval_batch => evaluate_darcy_anisotropic
            !! Evaluate the batched material response and trial history.
        procedure :: eval_sensitivity_batch => sensitivity_darcy_anisotropic
            !! Evaluate exact selected-parameter sensitivities through the discrete material update.
    end type darcy_anisotropic_t

contains

    subroutine describe_darcy_anisotropic(self, spec, status)
        !! Register a \(3\times3\) exact tangent and seven optional scalar derivatives.
        class(darcy_anisotropic_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema receiving the packed field descriptions and supported capabilities.
        integer, intent(out) :: status
            !! Status code; `cf_ok` indicates success.
        real(wp) :: permeability(ndim,ndim)
        real(wp) :: mobility(ndim,ndim)

        status = cf_invalid_argument
        if (.not. ieee_is_finite(self%dynamic_viscosity)) return
        if (self%dynamic_viscosity <= 0.0_wp) return
        call build_tensor(self, permeability, status)
        if (status /= cf_ok) return
        status = cf_invalid_argument
        mobility = permeability/self%dynamic_viscosity
        if (.not. all(ieee_is_finite(mobility))) return
        if (maxval(abs(mobility)) <= 0.0_wp) return

        call spec%add_input("pressure.gradient", ndim, "vector", "pressure_gradient", "spatial", [ndim])
        call spec%add_output("fluid.darcy_velocity", ndim, "vector", "specific_discharge", "spatial", [ndim])
        call spec%add_parameter("permeability_11")
        call spec%add_parameter("permeability_22")
        call spec%add_parameter("permeability_33")
        call spec%add_parameter("permeability_12")
        call spec%add_parameter("permeability_13")
        call spec%add_parameter("permeability_23")
        call spec%add_parameter("dynamic_viscosity")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_darcy_anisotropic

    subroutine evaluate_darcy_anisotropic(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Stateless spatial velocity evaluation over preallocated point arrays.
        class(darcy_anisotropic_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed inputs `(ninput,npoints)` when required, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `pressure.gradient` inputs `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued history `(nstate,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Computed uncommitted real-valued history `(nstate,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `fluid.darcy_velocity` responses `(noutput,npoints)`.
        real(wp), intent(in) :: dt
            !! Prescribed time increment; a positive value is required for time-dependent laws.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic tangent `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(inout) :: point_status(:)
            !! Per-point outcome; write `cf_ok` only for a successful constitutive update.

        state_trial = state_old
        call evaluate_points(self, x_new, response, point_status, tangent)
    end subroutine evaluate_darcy_anisotropic

    subroutine sensitivity_darcy_anisotropic(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Exact local derivatives for selected parameter names and order.
        !! `dsdp` has shape `(3,nselected,npoints)`, with zero state components.
        class(darcy_anisotropic_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed inputs `(ninput,npoints)` when required, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `pressure.gradient` inputs `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued history `(nstate,npoints)`.
        integer, intent(in) :: selected_parameters(:)
            !! One-based parameter IDs in the caller's selected derivative-axis order.
        real(wp), intent(in) :: sensitivity_old(:,:,:)
            !! Committed real-history parameter derivatives `(nstate,nselected,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Computed uncommitted real-valued history `(nstate,npoints)`.
        real(wp), intent(out) :: sensitivity_trial(:,:,:)
            !! Updated real-history parameter derivatives `(nstate,nselected,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `fluid.darcy_velocity` responses `(noutput,npoints)`.
        real(wp), intent(out) :: dsdp(:,:,:)
            !! Exact total response derivatives `(noutput,nselected,npoints)` at fixed prescribed driving histories.
        real(wp), intent(in) :: dt
            !! Prescribed time increment; a positive value is required for time-dependent laws.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic tangent `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(inout) :: point_status(:)
            !! Per-point outcome; write `cf_ok` only for a successful constitutive update.

        state_trial = state_old
        sensitivity_trial = sensitivity_old
        call evaluate_points(self, x_new, response, point_status, tangent, selected_parameters, dsdp)
    end subroutine sensitivity_darcy_anisotropic

    subroutine evaluate_points(self, gradient, response, point_status, tangent, selected_parameters, dsdp)
        !! Validate and form mobility once per batch. No allocation per point.
        class(darcy_anisotropic_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        real(wp), intent(in) :: gradient(:,:)
            !! Spatial `pressure.gradient` field at the material points.
        real(wp), intent(out) :: response(:,:)
            !! Packed `fluid.darcy_velocity` responses `(noutput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Per-point outcome; write `cf_ok` only for a successful constitutive update.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic tangent `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(in), optional :: selected_parameters(:)
            !! One-based parameter IDs in the caller's selected derivative-axis order.
        real(wp), intent(out), optional :: dsdp(:,:,:)
            !! Exact total response derivatives `(noutput,nselected,npoints)` at fixed prescribed driving histories.
        real(wp) :: permeability(ndim,ndim)
        real(wp) :: mobility(ndim,ndim)
        logical :: rescale_row(ndim)
        integer :: status
        integer :: ipoint
        integer :: iparam
        integer :: id
        integer :: i

        point_status = cf_invalid_argument
        call build_tensor(self,permeability,status)
        if (status /= cf_ok) return
        mobility = permeability/self%dynamic_viscosity
        if (.not. all(ieee_is_finite(mobility))) return
        do i = 1, ndim
            rescale_row(i) = any(abs(permeability(i,:)) > 0.0_wp .and. abs(mobility(i,:)) < tiny(1.0_wp))
        end do
        if (present(selected_parameters)) then
            if (any(selected_parameters < 1) .or. any(selected_parameters > nparameters)) return
        end if

        point_loop: do ipoint = 1, size(gradient,2)
            if (.not. all(ieee_is_finite(gradient(:,ipoint)))) cycle point_loop
            do i = 1, ndim
                response(i,ipoint) = -dot_product(mobility(i,:),gradient(:,ipoint))
                if (rescale_row(i) .or. .not. ieee_is_finite(response(i,ipoint)) .or. abs(response(i,ipoint)) < tiny(1.0_wp)) then
                    response(i,ipoint) = -scaled_darcy_product(permeability(i,:),gradient(:,ipoint),self%dynamic_viscosity,1)
                end if
            end do
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop

            if (present(tangent)) tangent(:,:,ipoint) = -mobility

            if (present(selected_parameters)) then
                do iparam = 1, size(selected_parameters)
                    id = selected_parameters(iparam)
                    dsdp(:,iparam,ipoint) = 0.0_wp
                    select case (id)
                    case (1:3)
                        dsdp(id,iparam,ipoint) = -gradient(id,ipoint)/self%dynamic_viscosity
                    case (4)
                        dsdp(1,iparam,ipoint) = -gradient(2,ipoint)/self%dynamic_viscosity
                        dsdp(2,iparam,ipoint) = -gradient(1,ipoint)/self%dynamic_viscosity
                    case (5)
                        dsdp(1,iparam,ipoint) = -gradient(3,ipoint)/self%dynamic_viscosity
                        dsdp(3,iparam,ipoint) = -gradient(1,ipoint)/self%dynamic_viscosity
                    case (6)
                        dsdp(2,iparam,ipoint) = -gradient(3,ipoint)/self%dynamic_viscosity
                        dsdp(3,iparam,ipoint) = -gradient(2,ipoint)/self%dynamic_viscosity
                    case (7)
                        do i = 1, ndim
                            dsdp(i,iparam,ipoint) = scaled_darcy_product(permeability(i,:),gradient(:,ipoint),self%dynamic_viscosity,2)
                        end do
                    case default
                        cycle point_loop
                    end select
                end do
                if (.not. all(ieee_is_finite(dsdp(:,:,ipoint)))) cycle point_loop
            end if
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_points

    pure function scaled_darcy_product(permeability_row, gradient, dynamic_viscosity, viscosity_power) result(value)
        !! Evaluate \((\boldsymbol{K}\boldsymbol{g})_i/\mu^p\), with exponent \(p\) given by `viscosity_power`, from original factors, preserving finite cancellation.
        real(wp), intent(in) :: permeability_row(:)
            !! Three coefficients in one spatial permeability-tensor row.
        real(wp), intent(in) :: gradient(ndim)
            !! Prescribed spatial pressure-gradient components.
        real(wp), intent(in) :: dynamic_viscosity
            !! Positive dynamic viscosity in stress times time.
        integer, intent(in) :: viscosity_power
            !! One for the response product or two for its viscosity derivative.
        real(wp) :: value
        real(wp) :: term_fraction(ndim)
        real(wp) :: sum_fraction
        integer :: term_exponent(ndim)
        integer :: sum_exponent
        integer :: common_exponent
        integer :: largest
        integer :: i
        integer :: j

        term_fraction = (fraction(permeability_row)*fraction(gradient))/fraction(dynamic_viscosity)**viscosity_power
        term_exponent = exponent(permeability_row)+exponent(gradient)-viscosity_power*exponent(dynamic_viscosity)+exponent(term_fraction)
        term_fraction = fraction(term_fraction)

        ! Cancel larger terms before adding smaller ones and restore the
        ! exponent only after the complete row sum has been formed.
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
        value = scale(sum_fraction,sum_exponent)
    end function scaled_darcy_product

    subroutine build_tensor(self, permeability, status)
        !! Scaled Cholesky validation rejects non-SPD tensors without clipping.
        !! Positive scaling avoids overflow in determinant products.
        class(darcy_anisotropic_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        real(wp), intent(out) :: permeability(ndim,ndim)
            !! Symmetric spatial permeability tensor assembled from model coefficients.
        integer, intent(out) :: status
            !! Status code; `cf_ok` indicates success.
        real(wp) :: scaled(ndim,ndim)
        real(wp) :: scale
        real(wp) :: l11
        real(wp) :: l21
        real(wp) :: l31
        real(wp) :: l22
        real(wp) :: l32
        real(wp) :: pivot

        status = cf_invalid_argument
        permeability(1,:) = [self%permeability_11,self%permeability_12, self%permeability_13]
        permeability(2,:) = [self%permeability_12,self%permeability_22, self%permeability_23]
        permeability(3,:) = [self%permeability_13,self%permeability_23, self%permeability_33]
        if (.not. all(ieee_is_finite(permeability))) return
        scale = maxval(abs(permeability))
        if (scale <= 0.0_wp) return
        scaled = permeability/scale

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

end module constituf_darcy_anisotropic

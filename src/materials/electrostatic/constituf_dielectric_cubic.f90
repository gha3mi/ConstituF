! SPDX-License-Identifier: Apache-2.0

module constituf_dielectric_cubic
    !! Isotropic cubic-field dielectric response in three spatial dimensions.
    !!
    !! This is a phenomenological, reversible *electrostatic* model with
    !! electric coenergy density in energy-per-volume units
    !! \[\mathcal{W}^{*}(\boldsymbol{E})=
    !! \frac{\varepsilon}{2}\boldsymbol{E}\cdot\boldsymbol{E}
    !! +\frac{\alpha}{4}(\boldsymbol{E}\cdot\boldsymbol{E})^2,\]
    !! where the prescribed electric field \(\boldsymbol{E}\) has dimensions
    !! of electric potential per length. Absolute linear permittivity
    !! \(\varepsilon>0\) has units of electric displacement per electric
    !! field. The cubic coefficient \(\alpha\geq0\) has units of electric
    !! displacement per electric field cubed. Electric displacement
    !! \(\boldsymbol{D}\), in charge per area, is
    !! \[\boldsymbol{D}=\partial_{\boldsymbol{E}}\mathcal{W}^{*}
    !! =(\varepsilon+\alpha\|\boldsymbol{E}\|^2)\boldsymbol{E}.\]
    !! The exact, symmetric constitutive Jacobian is
    !! \[\frac{\partial D_i}{\partial E_j}
    !! =(\varepsilon+\alpha\|\boldsymbol{E}\|^2)\delta_{ij}
    !! +2\alpha E_iE_j.\]
    !! The positive-coefficient domain makes the coenergy strictly
    !! convex and the true reversible stored energy (Legendre transform)
    !! is \(u=\boldsymbol{E}\cdot\boldsymbol{D}-\mathcal W^{*}
    !! =\varepsilon\|\boldsymbol{E}\|^2/2
    !! +3\alpha\|\boldsymbol{E}\|^4/4\geq0\).
    !! These energies define the constitutive law. The output field is
    !! electric displacement, without an additional energy output.
    !!
    !! Optional analytical *local* parameter sensitivities at fixed \(\boldsymbol{E}\) are
    !! \[\frac{\partial\boldsymbol{D}}{\partial\varepsilon}
    !! =\boldsymbol{E},\qquad
    !! \frac{\partial\boldsymbol{D}}{\partial\alpha}
    !! =\|\boldsymbol{E}\|^2\boldsymbol{E}.\]
    !! The input `electric.field` and output `electric.displacement` are
    !! spatial Cartesian 3-vectors, with no implicit unit conversions.
    !! The linear limit \(\alpha=0\) matches `dielectric_isotropic_t`.
    !!
    !! The quartic field coenergy is an explicit ConstituF *model choice*,
    !! not attributed to Jackson as an original formula and not a claim
    !! of an optical Kerr model. This model has no saturation, hysteresis,
    !! dynamics, polarization memory, conductivity, Maxwell stress, or
    !! electrostatic field solve. It assumes isothermal prescribed fields
    !! independent of selected model parameters. Nonfinite input, response
    !! or requested derivatives reject the full trial without clipping.
    !! Range checks apply to returned quantities, including only the
    !! derivatives selected by the caller. It is stateless and requires no `dt`.
    !!
    !! General electrostatic background (not copied code): J. D. Jackson,
    !! *Classical Electrodynamics*, 3rd ed., Wiley (1998),
    !! ISBN 978-0-471-30932-1, Chapter 4.
    !! Direct energy/coenergy convention: MIT, *Electromagnetic Fields
    !! and Energy*, section 11.4:
    !! https://web.mit.edu/6.013_book/www/chapter11/11.4.html .
    !! Code provenance: the polynomial response and analytical derivatives
    !! are independently derived from the stated coenergy. No third-party
    !! software code is reused in this revision.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: ndim = 3

    public :: dielectric_cubic_t

    type, extends(material_model_t) :: dielectric_cubic_t
        !! Create as `dielectric_cubic_t(permittivity,cubic_coefficient)`.
        real(wp) :: permittivity = 0.0_wp
            !! Absolute small-field permittivity \(\varepsilon>0\), electric displacement per electric field.
        real(wp) :: cubic_coefficient = 0.0_wp
            !! Nonnegative nonlinear coefficient \(\alpha\), electric displacement per electric field cubed.
    contains
        procedure :: describe => describe_dielectric_cubic
            !! Declare input and output fields, parameter names and constitutive capabilities.
        procedure :: eval_batch => evaluate_dielectric_cubic
            !! Evaluate the batched constitutive response and optional consistent tangent.
        procedure :: eval_sensitivity_batch => sensitivity_dielectric_cubic
            !! Evaluate exact scalar-parameter derivatives of the constitutive response.
    end type dielectric_cubic_t

contains

    subroutine describe_dielectric_cubic(self, spec, status)
        !! Register two named parameters and local exact derivatives.
        class(dielectric_cubic_t), intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        type(material_spec_t), intent(out) :: spec
            !! Packed field schema, material parameter metadata and declared capabilities.
        integer, intent(out) :: status
            !! Operation status; `cf_ok` indicates success.

        status = cf_invalid_argument
        if (.not. ieee_is_finite(self%permittivity)) return
        if (.not. ieee_is_finite(self%cubic_coefficient)) return
        if (self%permittivity <= 0.0_wp) return
        if (self%cubic_coefficient < 0.0_wp) return

        call spec%add_input("electric.field", ndim, "vector", "electric_field", "spatial", [ndim])
        call spec%add_output("electric.displacement", ndim, "vector", "electric_displacement", "spatial", [ndim])
        call spec%add_parameter("permittivity")
        call spec%add_parameter("cubic_coefficient")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_dielectric_cubic

    subroutine evaluate_dielectric_cubic(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Stateless point-batch evaluation, with optional \(3\times3\) tangent.
        class(dielectric_cubic_t), intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed input `(ninput,npoints)` when needed, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `electric.field` input `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued material history `(nstate,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Trial real-valued material history `(nstate,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `electric.displacement` output `(noutput,npoints)`.
        real(wp), intent(in) :: dt
            !! Prescribed time increment; required to be positive only for time-dependent models.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic response Jacobian `(noutput,ninput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Status of each point; a successful update must set its code to `cf_ok`.

        state_trial = state_old
        call evaluate_points(self,x_new,response,point_status,tangent)
    end subroutine evaluate_dielectric_cubic

    subroutine sensitivity_dielectric_cubic(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Exact selected local derivatives, with no stored derivative history.
        !! Parameter IDs are `permittivity=1` and `cubic_coefficient=2`.
        class(dielectric_cubic_t), intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed input `(ninput,npoints)` when needed, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `electric.field` input `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued material history `(nstate,npoints)`.
        integer, intent(in) :: selected_parameters(:)
            !! One-based IDs of scalar parameters in the requested sensitivity-axis order.
        real(wp), intent(in) :: sensitivity_old(:,:,:)
            !! Committed derivatives of real material history `(nstate,nselected,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Trial real-valued material history `(nstate,npoints)`.
        real(wp), intent(out) :: sensitivity_trial(:,:,:)
            !! Trial derivatives of real material history `(nstate,nselected,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `electric.displacement` output `(noutput,npoints)`.
        real(wp), intent(out) :: dsdp(:,:,:)
            !! Exact discrete response derivatives `(noutput,nselected,npoints)` at fixed input histories.
        real(wp), intent(in) :: dt
            !! Prescribed time increment; required to be positive only for time-dependent models.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic response Jacobian `(noutput,ninput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Status of each point; a successful update must set its code to `cf_ok`.

        state_trial = state_old
        sensitivity_trial = sensitivity_old
        call evaluate_points(self,x_new,response,point_status,tangent, selected_parameters,dsdp)
    end subroutine sensitivity_dielectric_cubic

    subroutine evaluate_points(self, electric, response, point_status, tangent, selected_parameters, dsdp)
        !! Fused allocation-free kernel. Derivatives computed only when asked.
        class(dielectric_cubic_t), intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        real(wp), intent(in) :: electric(:,:)
            !! Prescribed electric field `electric.field` at each point.
        real(wp), intent(out) :: response(:,:)
            !! Packed `electric.displacement` output `(noutput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Status of each point; a successful update must set its code to `cf_ok`.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic response Jacobian `(noutput,ninput,npoints)`.
        integer, intent(in), optional :: selected_parameters(:)
            !! One-based IDs of scalar parameters in the requested sensitivity-axis order.
        real(wp), intent(out), optional :: dsdp(:,:,:)
            !! Exact discrete response derivatives `(noutput,nselected,npoints)` at fixed input histories.
        real(wp) :: norm_squared
        real(wp) :: field_scale
        real(wp) :: scaled_norm_squared
        real(wp) :: diagonal_norm_squared
        real(wp) :: effective_permittivity
        real(wp) :: nonlinear_factor
        logical :: scaled_evaluation
        integer :: ipoint
        integer :: i
        integer :: j
        integer :: iparam

        point_status = cf_invalid_argument
        if (present(selected_parameters)) then
            if (any(selected_parameters < 1) .or. any(selected_parameters > 2)) return
        end if

        point_loop: do ipoint = 1, size(electric,2)
            if (.not. all(ieee_is_finite(electric(:,ipoint)))) cycle point_loop
            norm_squared = dot_product(electric(:,ipoint),electric(:,ipoint))
            scaled_evaluation = .not. ieee_is_finite(norm_squared)
            if (.not. scaled_evaluation) then
                if (norm_squared < tiny(1.0_wp)) scaled_evaluation = maxval(abs(electric(:,ipoint))) > 0.0_wp
            end if
            if (.not. scaled_evaluation) then
                nonlinear_factor = self%cubic_coefficient*norm_squared
                effective_permittivity = self%permittivity+nonlinear_factor
                scaled_evaluation = .not. ieee_is_finite(effective_permittivity)
                if (effective_permittivity < tiny(1.0_wp)) scaled_evaluation = .true.
            end if
            if (scaled_evaluation) then
                field_scale = maxval(abs(electric(:,ipoint)))
                scaled_norm_squared = 0.0_wp
                if (field_scale > 0.0_wp) scaled_norm_squared = sum((electric(:,ipoint)/field_scale)**2)
                response(:,ipoint) = scaled_norm_product(self%cubic_coefficient,field_scale,scaled_norm_squared,electric(:,ipoint),self%permittivity)
            else
                response(:,ipoint) = effective_permittivity*electric(:,ipoint)
            end if
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop

            if (present(tangent)) then
                do j = 1, ndim
                    do i = 1, ndim
                        tangent(i,j,ipoint) = scaled_outer_product(self%cubic_coefficient,electric(i,ipoint),electric(j,ipoint))
                        if (i == j) then
                            if (scaled_evaluation) then
                                diagonal_norm_squared = scaled_norm_squared
                                if (field_scale > 0.0_wp) diagonal_norm_squared = diagonal_norm_squared+2.0_wp*(electric(i,ipoint)/field_scale)**2
                                tangent(i,j,ipoint) = scaled_norm_product(self%cubic_coefficient,field_scale,diagonal_norm_squared,1.0_wp,self%permittivity)
                            else
                                tangent(i,j,ipoint) = tangent(i,j,ipoint)+effective_permittivity
                            end if
                        end if
                    end do
                end do
                if (.not. all(ieee_is_finite(tangent(:,:,ipoint)))) cycle point_loop
            end if

            if (present(selected_parameters)) then
                do iparam = 1, size(selected_parameters)
                    select case (selected_parameters(iparam))
                    case (1)
                        dsdp(:,iparam,ipoint) = electric(:,ipoint)
                    case (2)
                        if (scaled_evaluation) then
                            dsdp(:,iparam,ipoint) = scaled_norm_product(1.0_wp,field_scale,scaled_norm_squared,electric(:,ipoint))
                        else
                            dsdp(:,iparam,ipoint) = norm_squared*electric(:,ipoint)
                        end if
                    case default
                        cycle point_loop
                    end select
                end do
                if (.not. all(ieee_is_finite(dsdp(:,:,ipoint)))) cycle point_loop
            end if
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_points

    pure elemental function scaled_norm_product(coefficient, field_scale, scaled_norm_squared, component, linear_permittivity) result(product)
        !! Form \(c\|\boldsymbol{E}\|^2u\), with `coefficient` \(c\) and `component` \(u\), using one final exponent restoration.
        !! Retaining the original component preserves small nonzero entries.
        !! Combine an optional linear permittivity before subnormal rounding.
        real(wp), intent(in) :: coefficient
            !! Cubic coefficient, or unity for its parameter derivative.
        real(wp), intent(in) :: field_scale
            !! Largest absolute prescribed electric-field component.
        real(wp), intent(in) :: scaled_norm_squared
            !! Squared normalized field norm, including a diagonal tangent weight when needed.
        real(wp), intent(in) :: component
            !! Original field component, or unity for the quadratic factor.
        real(wp), intent(in), optional :: linear_permittivity
            !! Positive linear coefficient added to the quadratic factor.
        real(wp) :: product
        real(wp) :: quadratic_fraction
        real(wp) :: coefficient_fraction
        integer :: quadratic_exponent
        integer :: coefficient_exponent

        quadratic_fraction = fraction(coefficient)*fraction(field_scale)**2
        quadratic_exponent = exponent(coefficient)+2*exponent(field_scale)
        if (present(linear_permittivity)) then
            if (coefficient <= 0.0_wp .or. field_scale <= 0.0_wp .or. scaled_norm_squared <= 0.0_wp) then
                product = linear_permittivity*component
                return
            end if
            coefficient_exponent = max(exponent(linear_permittivity),quadratic_exponent)
            coefficient_fraction = scale(fraction(linear_permittivity),exponent(linear_permittivity)-coefficient_exponent)+scale(quadratic_fraction*scaled_norm_squared,quadratic_exponent-coefficient_exponent)
            product = scale(coefficient_fraction*fraction(component),coefficient_exponent+exponent(component))
        else
            product = scale((quadratic_fraction*fraction(component))*scaled_norm_squared,quadratic_exponent+exponent(component))
        end if
    end function scaled_norm_product

    pure elemental function scaled_outer_product(coefficient, field_i, field_j) result(product)
        !! Form the exact tangent contribution \(2\alpha E_iE_j\) without
        !! requiring either \(2\alpha\) or a field product to be representable.
        real(wp), intent(in) :: coefficient
            !! Nonnegative cubic coefficient.
        real(wp), intent(in) :: field_i
            !! Original electric-field component \(E_i\).
        real(wp), intent(in) :: field_j
            !! Original electric-field component \(E_j\).
        real(wp) :: product

        product = scale((2.0_wp*fraction(coefficient)*fraction(field_i))*fraction(field_j),exponent(coefficient)+exponent(field_i)+exponent(field_j))
    end function scaled_outer_product

end module constituf_dielectric_cubic

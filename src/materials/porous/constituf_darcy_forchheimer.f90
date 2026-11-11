! SPDX-License-Identifier: Apache-2.0

module constituf_darcy_forchheimer
    !! Isotropic Darcy--Forchheimer flow through a rigid saturated medium.
    !!
    !! A prescribed spatial pressure gradient \(\boldsymbol{g}=\nabla p\)
    !! in pressure per length drives the superficial velocity \(\boldsymbol{v}\) in length per time:
    !! \[\boldsymbol{g}=-(a+b\|\boldsymbol{v}\|)\boldsymbol{v}.\]
    !! The constant linear resistance \(a>0\) has dimensions of pressure
    !! times time per length squared. The nonnegative quadratic resistance
    !! \(b\geq0\) has dimensions of pressure times time squared per length
    !! cubed. Both parameters must be finite. Initialization also
    !! requires \(1/a\) to be finite and strictly positive in the working
    !! precision. In the usual fluid/porous parametrization,
    !! \(a=\mu/k\) and \(b=\rho\beta\), with dynamic viscosity
    !! \(\mu\), intrinsic permeability \(k\), density \(\rho\)
    !! and inertial coefficient \(\beta\) in inverse length. This implementation accepts
    !! the two *lumped resistances* directly and does not infer them.
    !!
    !! For \(r=\|\boldsymbol{g}\|\), solve the scalar quadratic exactly,
    !! \[b s^2+a s=r,\quad s=\|\boldsymbol{v}\|,\qquad
    !! \boldsymbol{v}=-m(r)\boldsymbol{g},\quad
    !! m(r)=\frac{2}{a+\sqrt{a^2+4br}}.\]
    !! The implementation carries the gradient norm and
    !! \(d=\sqrt{a^2+4br}\) as scaled factors and evaluates their inverse
    !! products without requiring either magnitude to be representable.
    !! The positive root avoids subtractive quadratic cancellation.
    !! At zero gradient, the velocity is zero and the tangent is
    !! \(-\boldsymbol{I}/a\), the continuous zero-field limit.
    !!
    !! At nonzero \(r\) let \(\boldsymbol{n}=\boldsymbol{g}/r\).
    !! The exact local algorithmic tangent at fixed parameters is
    !! \[\frac{\partial\boldsymbol{v}}{\partial\boldsymbol{g}}
    !! =-m\boldsymbol{I}+(m-d^{-1})
    !! \boldsymbol{n}\otimes\boldsymbol{n}.\]
    !! The coefficient \(m-d^{-1}=brm^2/d\) is evaluated without
    !! subtracting nearly equal mobilities. Original field components
    !! are retained in the tangent and sensitivity products.
    !! Optional exact local parameter sensitivities, holding the supplied
    !! pressure gradient fixed, are
    !! \[\frac{\partial\boldsymbol{v}}{\partial a}
    !! =-\frac{\boldsymbol{v}}{d},\qquad
    !! \frac{\partial\boldsymbol{v}}{\partial b}
    !! =-\frac{s}{d}\boldsymbol{v}.\]
    !! At \(b=0\) the quadratic-resistance sensitivity is the right derivative.
    !! The pressure-work dissipation rate per unit bulk volume (power per volume) obeys
    !! \[-\boldsymbol{g}\cdot\boldsymbol{v}
    !! =a s^2+b s^3\geq0.\]
    !! This derived rate is not a returned response field.
    !! The \(b=0\) limit is exactly the Darcy law
    !! \(\boldsymbol{v}=-(1/a)\boldsymbol{g}\).
    !!
    !! Fields `pressure.gradient` (pressure per length) and `fluid.darcy_velocity`
    !! (length per time) are spatial 3-vectors in the same fixed frame. The velocity is
    !! specific discharge relative to the stationary porous skeleton.
    !! No pore velocity, porosity or local acceleration is computed.
    !! Gravity is omitted for the stated pressure-gradient input. A caller
    !! including a body acceleration \(\boldsymbol b_{\mathrm{acc}}\) may
    !! supply \(\boldsymbol g_{\mathrm{eff}}=\nabla p-\rho\boldsymbol b_{\mathrm{acc}}\)
    !! externally in the same field and units. Both the resistance law and
    !! dissipation then use this effective driving gradient. Hydrostatic
    !! head, relative permeability, saturation changes, temperature,
    !! transients, mechanical deformation and global mass conservation
    !! are not implemented.
    !! This is a stateless algebraic closure, not a flow solver.
    !!
    !! Parameter-independent prescribed fields are assumed for local
    !! sensitivities. Nonfinite inputs or nonrepresentable returned
    !! velocities, tangents or requested derivatives fail the trial
    !! instead of silently clipping. All work is allocation-free in
    !! the batch kernel. No nonlinear iteration is required.
    !!
    !! Scientific sources (not code): P. Forchheimer, "Wasserbewegung
    !! durch Boden", Zeitschrift des Vereines deutscher Ingenieure,
    !! vol. 45, no. 50 (1901), pp. 1781--1788. H. Darcy,
    !! *Les fontaines publiques de la ville de Dijon* (1856).
    !! The modern isotropic vector law and the mapping
    !! \(b=\rho c_F/\sqrt{k}\) appear in the COMSOL Porous Media Flow Module
    !! User's Guide 6.1 (2022), "Permeability Models--Forchheimer",
    !! Eq. (2-6), p. 32:
    !! https://doc.comsol.com/6.1/doc/com.comsol.help.porous/PorousMediaFlowModuleUsersGuide.pdf .
    !! The algebraic inversion and derivatives are independently derived
    !! from this standard lumped-resistance law. No third-party software
    !! code is reused in this revision.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: ndim = 3

    public :: darcy_forchheimer_t

    type, extends(material_model_t) :: darcy_forchheimer_t
        !! Construct as `darcy_forchheimer_t(linear_resistance,quadratic_resistance)`.
        !! Positive linear resistance is essential for finite zero-field tangent.
        real(wp) :: linear_resistance = 0.0_wp
            !! Positive resistance \(a\), pressure times time per length squared.
        real(wp) :: quadratic_resistance = 0.0_wp
            !! Nonnegative resistance \(b\), pressure times time squared per length cubed.
    contains
        procedure :: describe => describe_forchheimer
            !! Declare material fields, parameters and supported operations.
        procedure :: eval_batch => evaluate_forchheimer
            !! Evaluate the batched material response, optional tangent and trial state.
        procedure :: eval_sensitivity_batch => sensitivity_forchheimer
            !! Evaluate analytical parameter derivatives through the discrete material update.
    end type darcy_forchheimer_t

contains

    subroutine describe_forchheimer(self, spec, status)
        !! Define pressure gradient, specific discharge and two named parameters.
        class(darcy_forchheimer_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema populated with fields, physical metadata and model capabilities.
        integer, intent(out) :: status
            !! Operation status code. `cf_ok` indicates success.
        real(wp) :: zero_field_mobility

        status = cf_invalid_argument
        if (.not. all(ieee_is_finite([self%linear_resistance, self%quadratic_resistance]))) return
        if (self%linear_resistance <= 0.0_wp) return
        if (self%quadratic_resistance < 0.0_wp) return
        zero_field_mobility = 1.0_wp/self%linear_resistance
        if (.not. ieee_is_finite(zero_field_mobility)) return
        if (zero_field_mobility <= 0.0_wp) return

        call spec%add_input("pressure.gradient", ndim, "vector", "pressure_gradient", "spatial", [ndim])
        call spec%add_output("fluid.darcy_velocity", ndim, "vector", "specific_discharge", "spatial", [ndim])
        call spec%add_parameter("linear_resistance")
        call spec%add_parameter("quadratic_resistance")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_forchheimer

    subroutine evaluate_forchheimer(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Stateless closed-form Forchheimer inversion for a point batch.
        class(darcy_forchheimer_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previous packed inputs `(ninput,npoints)` when required by the model, otherwise a zero-sized array.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `pressure.gradient` inputs `(ninput,npoints)` at the current time step.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued history `(nstate,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Updated trial real-valued history `(nstate,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `fluid.darcy_velocity` outputs `(noutput,npoints)`.
        real(wp), intent(in) :: dt
            !! Prescribed time increment, strictly positive for time-dependent materials.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact discrete input Jacobian `(noutput,ninput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Per-point status codes. successful points must be assigned `cf_ok`.

        state_trial = state_old
        call evaluate_points(self, x_new, response, point_status, tangent)
    end subroutine evaluate_forchheimer

    subroutine sensitivity_forchheimer(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Exact local derivatives in arbitrary selection order, no state.
        !! Parameter IDs: `linear_resistance=1` and `quadratic_resistance=2`.
        class(darcy_forchheimer_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previous packed inputs `(ninput,npoints)` when required by the model, otherwise a zero-sized array.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `pressure.gradient` inputs `(ninput,npoints)` at the current time step.
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
            !! Packed `fluid.darcy_velocity` outputs `(noutput,npoints)`.
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
    end subroutine sensitivity_forchheimer

    subroutine evaluate_points(self, gradient, response, point_status, tangent, selected_parameters, dsdp)
        !! Closed form without quadratic cancellation or iteration.
        class(darcy_forchheimer_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        real(wp), intent(in) :: gradient(:,:)
            !! Spatial `pressure.gradient` values at each material point.
        real(wp), intent(out) :: response(:,:)
            !! Packed `fluid.darcy_velocity` outputs `(noutput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Per-point status codes. successful points must be assigned `cf_ok`.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact discrete input Jacobian `(noutput,ninput,npoints)`.
        integer, intent(in), optional :: selected_parameters(:)
            !! One-based IDs of scalar parameters in the requested derivative-axis order.
        real(wp), intent(out), optional :: dsdp(:,:,:)
            !! Exact response derivatives `(noutput,nselected,npoints)` at fixed prescribed input histories.
        real(wp) :: a
        real(wp) :: b
        real(wp) :: gradient_scale
        real(wp) :: scaled_radius
        real(wp) :: scaled_a
        real(wp) :: scaled_quadratic_root
        real(wp) :: root_fraction
        real(wp) :: scaled_root
        real(wp) :: mobility
        real(wp) :: radial_mobility
        real(wp) :: value
        integer :: product_exponent
        integer :: odd_exponent
        integer :: root_exponent
        integer :: common_exponent
        integer :: ipoint
        integer :: i
        integer :: j
        integer :: iparam

        point_status = cf_invalid_argument
        if (present(selected_parameters)) then
            if (any(selected_parameters < 1) .or. any(selected_parameters > 2)) return
        end if

        a = self%linear_resistance
        b = self%quadratic_resistance
        point_loop: do ipoint = 1, size(gradient,2)
            if (.not. all(ieee_is_finite(gradient(:,ipoint)))) cycle point_loop
            gradient_scale = maxval(abs(gradient(:,ipoint)))
            if (gradient_scale <= 0.0_wp) then
                response(:,ipoint) = 0.0_wp
                if (present(tangent)) then
                    tangent(:,:,ipoint) = 0.0_wp
                    do i = 1, ndim
                        tangent(i,i,ipoint) = -1.0_wp/a
                    end do
                end if
                if (present(selected_parameters)) dsdp(:,:,ipoint) = 0.0_wp
                point_status(ipoint) = cf_ok
                cycle point_loop
            end if
            scaled_radius = norm2(gradient(:,ipoint)/gradient_scale)

            if (b > 0.0_wp) then
                ! Represent 2*sqrt(b*r) before aligning it with a. Neither
                ! the gradient norm nor the physical discriminant is formed.
                product_exponent = exponent(b)+exponent(gradient_scale)
                odd_exponent = modulo(product_exponent,2)
                root_exponent = (product_exponent-odd_exponent)/2
                root_fraction = 2.0_wp*sqrt(scale(fraction(b)*fraction(gradient_scale)*scaled_radius,odd_exponent))
                common_exponent = max(exponent(a),root_exponent)
                scaled_a = scale(fraction(a),exponent(a)-common_exponent)
                scaled_quadratic_root = scale(root_fraction,root_exponent-common_exponent)
                scaled_root = sqrt(scaled_a**2+scaled_quadratic_root**2)
                radial_mobility = scale(1.0_wp/scaled_root,-common_exponent)
                mobility = scale((2.0_wp/(1.0_wp+scaled_a/scaled_root))/scaled_root,-common_exponent)
            else
                mobility = 1.0_wp/a
                radial_mobility = mobility
            end if
            if (.not. ieee_is_finite(mobility)) cycle point_loop
            if (mobility <= 0.0_wp) cycle point_loop
            response(:,ipoint) = -mobility*gradient(:,ipoint)
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop

            if (present(tangent)) then
                tangent(:,:,ipoint) = 0.0_wp
                do i = 1, ndim
                    tangent(i,i,ipoint) = -mobility
                end do
                if (b > 0.0_wp) then
                    do j = 1, ndim
                        do i = 1, ndim
                            value = scaled_direction_product(b,mobility,radial_mobility,gradient(i,ipoint),gradient(j,ipoint),gradient_scale,scaled_radius)
                            tangent(i,j,ipoint) = tangent(i,j,ipoint)+value
                        end do
                    end do
                end if
                if (.not. all(ieee_is_finite(tangent(:,:,ipoint)))) cycle point_loop
            end if

            if (present(selected_parameters)) then
                do iparam = 1, size(selected_parameters)
                    select case (selected_parameters(iparam))
                    case (1)
                        dsdp(:,iparam,ipoint) = scaled_linear_derivative(mobility,radial_mobility,gradient(:,ipoint))
                    case (2)
                        dsdp(:,iparam,ipoint) = scaled_quadratic_product(gradient_scale,scaled_radius,mobility,radial_mobility,gradient(:,ipoint))
                    case default
                        cycle point_loop
                    end select
                end do
                if (.not. all(ieee_is_finite(dsdp(:,:,ipoint)))) cycle point_loop
            end if
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_points

    pure elemental function scaled_linear_derivative(mobility, radial_mobility, gradient) result(value)
        !! Evaluate \(m g_i/d\) without an underflowed or overflowing intermediate product.
        real(wp), intent(in) :: mobility
            !! Transverse mobility \(m\).
        real(wp), intent(in) :: radial_mobility
            !! Radial mobility \(1/d\).
        real(wp), intent(in) :: gradient
            !! One original pressure-gradient component.
        real(wp) :: value

        value = scale((fraction(mobility)*fraction(radial_mobility))*fraction(gradient), exponent(mobility)+exponent(radial_mobility)+exponent(gradient))
    end function scaled_linear_derivative

    pure elemental function scaled_quadratic_product(gradient_scale, scaled_radius, mobility, radial_mobility, factor) result(value)
        !! Evaluate \(r m^2 f/d\), where \(f\) is the supplied `factor`, from the original scaled radius and factor.
        real(wp), intent(in) :: gradient_scale
            !! Maximum absolute pressure-gradient component.
        real(wp), intent(in) :: scaled_radius
            !! Norm of the gradient divided by `gradient_scale`, between \(1\) and \(\sqrt{3}\).
        real(wp), intent(in) :: mobility
            !! Transverse mobility \(m\).
        real(wp), intent(in) :: radial_mobility
            !! Radial mobility \(1/d\).
        real(wp), intent(in) :: factor
            !! Original gradient component for the quadratic-resistance derivative.
        real(wp) :: value

        value = scale(((fraction(gradient_scale)*scaled_radius)*fraction(mobility)**2)*fraction(radial_mobility)*fraction(factor), exponent(gradient_scale)+2*exponent(mobility)+exponent(radial_mobility)+exponent(factor))
    end function scaled_quadratic_product

    pure elemental function scaled_direction_product(quadratic_resistance, mobility, radial_mobility, first, second, gradient_scale, scaled_radius) result(value)
        !! Evaluate \(b m^2 g_i g_j/(d r)\) before rounding the tangent weight.
        real(wp), intent(in) :: quadratic_resistance
            !! Nonnegative quadratic resistance \(b\).
        real(wp), intent(in) :: mobility
            !! Transverse mobility \(m\).
        real(wp), intent(in) :: radial_mobility
            !! Radial mobility \(1/d\).
        real(wp), intent(in) :: first
            !! Original pressure-gradient component \(g_i\).
        real(wp), intent(in) :: second
            !! Original pressure-gradient component \(g_j\).
        real(wp), intent(in) :: gradient_scale
            !! Maximum absolute pressure-gradient component.
        real(wp), intent(in) :: scaled_radius
            !! Norm of the gradient divided by `gradient_scale`.
        real(wp) :: value

        value = scale(((fraction(quadratic_resistance)*fraction(mobility)**2)*fraction(radial_mobility))*(fraction(first)*fraction(second))/(fraction(gradient_scale)*scaled_radius), exponent(quadratic_resistance)+2*exponent(mobility)+exponent(radial_mobility)+exponent(first)+exponent(second)-exponent(gradient_scale))
    end function scaled_direction_product

end module constituf_darcy_forchheimer

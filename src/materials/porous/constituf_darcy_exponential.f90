! SPDX-License-Identifier: Apache-2.0

module constituf_darcy_exponential
    !! Isotropic Darcy flow with phenomenological pressure-dependent permeability.
    !!
    !! At prescribed absolute pressure \(p\) (stress units) and spatial gradient
    !! \(\boldsymbol{g}=\nabla p\) (pressure per length), the superficial velocity is
    !! \[\boldsymbol{v}=-m(p)\boldsymbol{g},\qquad
    !! m(p)=\frac{k_0}{\mu}\exp[\beta(p-p_0)].\]
    !! The strictly positive reference intrinsic permeability \(k_0\)
    !! has units of length squared. Dynamic viscosity \(\mu>0\) has units of stress times time.
    !! \(p_0\geq0\) is the reference absolute pressure in the chosen stress unit and
    !! \(\beta\) is a finite pressure coefficient in inverse pressure (either sign).
    !! The reference mobility \(k_0/\mu\) must have a finite, positive
    !! normally representable value. The exponential dependence on
    !! pressure is a *ConstituF phenomenological choice*, not claimed
    !! to be part of Darcy's original constant-permeability experiments.
    !! It sets \(k(p_0)=k_0\) and \(\partial\ln k/\partial p=\beta\).
    !! Pressure and reference pressure use the same absolute datum.
    !!
    !! The packed independent input
    !! \([\partial_1p,\partial_2p,\partial_3p,p]\) comprises two fields,
    !! `pressure.gradient` (spatial 3-vector in pressure per length) and
    !! `pressure.absolute` (material scalar in the chosen stress unit), in that order.
    !! The output `fluid.darcy_velocity` is a spatial 3-vector in length per time,
    !! representing specific discharge relative to a stationary porous
    !! skeleton, not interstitial pore velocity.
    !! For the \(3\times4\) exact local input Jacobian at fixed parameters:
    !! \[\frac{\partial v_i}{\partial g_j}=-m(p)\delta_{ij},\quad
    !! \frac{\partial\boldsymbol{v}}{\partial p}=\beta\boldsymbol{v}.\]
    !! Optional exact local parameter sensitivities at fixed prescribed
    !! pressure and pressure gradient are
    !! \[\frac{\partial\boldsymbol{v}}{\partial k_0}
    !! =\frac{\boldsymbol{v}}{k_0},\quad
    !! \frac{\partial\boldsymbol{v}}{\partial\mu}
    !! =-\frac{\boldsymbol{v}}{\mu},\quad
    !! \frac{\partial\boldsymbol{v}}{\partial\beta}
    !! =(p-p_0)\boldsymbol{v},\quad
    !! \frac{\partial\boldsymbol{v}}{\partial p_0}
    !! =-\beta\boldsymbol{v}.\]
    !! Positive mobility ensures a nonnegative pressure-work dissipation
    !! rate per unit bulk volume (power per volume):
    !! \[-\boldsymbol{g}\cdot\boldsymbol{v}
    !! =m(p)\|\boldsymbol{g}\|^2\geq0.\]
    !! This derived rate is not a returned response field.
    !! At \(\beta=0\), this is the constant-permeability isotropic
    !! Darcy model. There is no history or time-increment dependence.
    !!
    !! Absolute pressure must be nonnegative and finite. Exponentials
    !! outside the normally representable mobility range are rejected,
    !! never clipped to a physical or numerical cap. Nonfinite flux,
    !! tangent or selected parameter derivatives invalidate the trial.
    !! Gravity is omitted for the stated pressure-gradient input. A caller
    !! including a body acceleration \(\boldsymbol b\) may supply
    !! \(\boldsymbol g_{\mathrm{eff}}=\nabla p-\rho\boldsymbol b\)
    !! externally in the gradient field and the same units. The flux law
    !! and dissipation then use this effective driving gradient, while
    !! the separate absolute pressure controls mobility. The material
    !! models neither saturation nor porosity evolution, fluid
    !! compressibility, hysteresis, multiphase effects,
    !! mechanical coupling or global mass conservation. No global solver
    !! or implicit field-parameter sensitivity is provided. The
    !! batch integration loop is allocation-free, I/O-free and stateless.
    !!
    !! Original Darcy flow: H. Darcy, *Les fontaines publiques de la
    !! ville de Dijon*, Dalmont, Paris (1856), Appendix D.
    !! The modern pressure-gradient form and specific-discharge convention
    !! appear in the COMSOL Porous Media Flow Module User's Guide 6.1
    !! (2022), "Darcy's Law--Equation Formulation", Eqs. (4-1) and (4-3),
    !! pp. 114-116:
    !! https://doc.comsol.com/6.1/doc/com.comsol.help.porous/PorousMediaFlowModuleUsersGuide.pdf .
    !! The exponential pressure parametrization is the stated independent
    !! adaptation. No third-party software code is reused in this revision.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: ndim = 3
    integer, parameter :: ninput = ndim + 1
    integer, parameter :: pressure_component = ninput

    public :: darcy_exponential_t

    type, extends(material_model_t) :: darcy_exponential_t
        !! Construct as `darcy_exponential_t(k0,mu,beta,p0)`, using a coherent unit system.
        real(wp) :: reference_permeability = 0.0_wp
            !! Intrinsic permeability \(k_0>0\), in squared-length units.
        real(wp) :: dynamic_viscosity = 0.0_wp
            !! Positive dynamic viscosity \(\mu\), stress times time.
        real(wp) :: pressure_coefficient = 0.0_wp
            !! Exponential pressure coefficient \(\beta\), inverse pressure.
        real(wp) :: reference_pressure = 0.0_wp
            !! Reference absolute pressure \(p_0\geq0\), in the chosen stress unit.
    contains
        procedure :: describe => describe_darcy_exponential
            !! Describe the material input, response, state and scalar parameter schema.
        procedure :: eval_batch => evaluate_darcy_exponential
            !! Evaluate the batched material response and trial history.
        procedure :: eval_sensitivity_batch => sensitivity_darcy_exponential
            !! Evaluate exact selected-parameter sensitivities through the discrete material update.
    end type darcy_exponential_t

contains

    subroutine describe_darcy_exponential(self, spec, status)
        !! Declare the exact \(3\times4\) Jacobian and four named scalar parameters.
        class(darcy_exponential_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema receiving the packed field descriptions and supported capabilities.
        integer, intent(out) :: status
            !! Status code; `cf_ok` indicates success.
        real(wp) :: reference_mobility

        status = cf_invalid_argument
        if (.not. all(ieee_is_finite([self%reference_permeability, self%dynamic_viscosity,self%pressure_coefficient, self%reference_pressure]))) return
        if (self%reference_permeability <= 0.0_wp) return
        if (self%dynamic_viscosity <= 0.0_wp) return
        if (self%reference_pressure < 0.0_wp) return
        reference_mobility = self%reference_permeability/self%dynamic_viscosity
        if (.not. ieee_is_finite(reference_mobility)) return
        if (reference_mobility < tiny(1.0_wp)) return

        call spec%add_input("pressure.gradient", ndim, "vector", "pressure_gradient", "spatial", [ndim])
        call spec%add_input("pressure.absolute", 1, "scalar", "absolute_pressure", "material")
        call spec%add_output("fluid.darcy_velocity", ndim, "vector", "specific_discharge", "spatial", [ndim])
        call spec%add_parameter("reference_permeability")
        call spec%add_parameter("dynamic_viscosity")
        call spec%add_parameter("pressure_coefficient")
        call spec%add_parameter("reference_pressure")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_darcy_exponential

    subroutine evaluate_darcy_exponential(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Evaluate pressure-dependent Darcy velocity on preallocated batches.
        class(darcy_exponential_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed inputs `(ninput,npoints)` when required, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `pressure.gradient` and `pressure.absolute` inputs `(ninput,npoints)`.
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
        call evaluate_points(self,x_new,response,point_status,tangent)
    end subroutine evaluate_darcy_exponential

    subroutine sensitivity_darcy_exponential(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Exact local derivatives for four independently selected coefficients.
        !! `dsdp(3,nselected,npoints)` has no retained history components.
        class(darcy_exponential_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed inputs `(ninput,npoints)` when required, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `pressure.gradient` and `pressure.absolute` inputs `(ninput,npoints)`.
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
        call evaluate_points(self,x_new,response,point_status,tangent, selected_parameters,dsdp)
    end subroutine sensitivity_darcy_exponential

    subroutine evaluate_points(self, x_new, response, point_status, tangent, selected_parameters, dsdp)
        !! Bounded exponential factors preserve the reference ratio and normal mobility domain.
        !! Differentiation skips all unrequested analytical parameter work.
        class(darcy_exponential_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `pressure.gradient` and `pressure.absolute` inputs `(ninput,npoints)`.
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
        real(wp) :: reference_mobility
        real(wp) :: beta
        real(wp) :: delta_pressure
        real(wp) :: pressure_exponent
        real(wp) :: exponent_limit
        real(wp) :: exponential_factor
        real(wp) :: mobility
        integer :: ipoint
        integer :: iparam
        integer :: i

        point_status = cf_invalid_argument
        if (present(selected_parameters)) then
            if (any(selected_parameters < 1) .or. any(selected_parameters > 4)) return
        end if

        reference_mobility = self%reference_permeability/self%dynamic_viscosity
        exponent_limit = log(huge(1.0_wp))-log(tiny(1.0_wp))
        beta = self%pressure_coefficient

        point_loop: do ipoint = 1, size(x_new,2)
            if (.not. all(ieee_is_finite(x_new(:,ipoint)))) cycle point_loop
            if (x_new(pressure_component,ipoint) < 0.0_wp) cycle point_loop
            delta_pressure = x_new(pressure_component,ipoint) - self%reference_pressure
            if (abs(beta) > 1.0_wp) then
                if (abs(delta_pressure) > huge(1.0_wp)/abs(beta)) cycle point_loop
            end if
            pressure_exponent = beta*delta_pressure
            if (.not. ieee_is_finite(pressure_exponent)) cycle point_loop
            if (abs(pressure_exponent) > exponent_limit) cycle point_loop
            exponential_factor = exp(pressure_exponent/4.0_wp)
            mobility = scale(fraction(reference_mobility)*fraction(exponential_factor)**4, exponent(reference_mobility)+4*exponent(exponential_factor))
            if (.not. ieee_is_finite(mobility)) cycle point_loop
            if (mobility < tiny(1.0_wp)) cycle point_loop
            response(:,ipoint) = -mobility*x_new(1:ndim,ipoint)
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop

            if (present(tangent)) then
                tangent(:,:,ipoint) = 0.0_wp
                do i = 1, ndim
                    tangent(i,i,ipoint) = -mobility
                end do
                tangent(:,pressure_component,ipoint) = scaled_derivative_product(-mobility,beta,x_new(1:ndim,ipoint))
                if (.not. all(ieee_is_finite(tangent(:,:,ipoint)))) cycle point_loop
            end if

            if (present(selected_parameters)) then
                do iparam = 1, size(selected_parameters)
                    select case (selected_parameters(iparam))
                    case (1)
                        dsdp(:,iparam,ipoint) = scaled_parameter_derivative(-mobility,self%reference_permeability,x_new(1:ndim,ipoint))
                    case (2)
                        dsdp(:,iparam,ipoint) = scaled_parameter_derivative(mobility,self%dynamic_viscosity,x_new(1:ndim,ipoint))
                    case (3)
                        dsdp(:,iparam,ipoint) = scaled_derivative_product(-mobility,delta_pressure,x_new(1:ndim,ipoint))
                    case (4)
                        dsdp(:,iparam,ipoint) = scaled_derivative_product(mobility,beta,x_new(1:ndim,ipoint))
                    case default
                        cycle point_loop
                    end select
                end do
                if (.not. all(ieee_is_finite(dsdp(:,:,ipoint)))) cycle point_loop
            end if
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_points

    pure elemental function scaled_derivative_product(mobility, factor, gradient) result(derivative)
        !! Preserve a finite three-factor derivative even when velocity underflows.
        real(wp), intent(in) :: mobility
            !! Signed evaluated mobility.
        real(wp), intent(in) :: factor
            !! Pressure coefficient or prescribed pressure difference.
        real(wp), intent(in) :: gradient
            !! One original pressure-gradient component.
        real(wp) :: derivative

        derivative = scale((fraction(mobility)*fraction(factor))*fraction(gradient), exponent(mobility)+exponent(factor)+exponent(gradient))
    end function scaled_derivative_product

    pure elemental function scaled_parameter_derivative(mobility, parameter, gradient) result(derivative)
        !! Form signed \(m g_i/\theta\), where \(\theta\) is the selected `parameter`, without first rounding the velocity or ratio.
        real(wp), intent(in) :: mobility
            !! Signed evaluated mobility.
        real(wp), intent(in) :: parameter
            !! Strictly positive reference permeability or dynamic viscosity.
        real(wp), intent(in) :: gradient
            !! One original pressure-gradient component.
        real(wp) :: derivative

        derivative = scale((fraction(mobility)*fraction(gradient))/fraction(parameter), exponent(mobility)+exponent(gradient)-exponent(parameter))
    end function scaled_parameter_derivative

end module constituf_darcy_exponential

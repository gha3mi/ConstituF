! SPDX-License-Identifier: Apache-2.0

module constituf_fourier_exponential
    !! Three-dimensional isotropic Fourier conduction with exponential \(k(T)\).
    !!
    !! The prescribed absolute temperature \(T\) and spatial gradient
    !! \(\boldsymbol{g}=\nabla T\), with dimensions of temperature interval
    !! per length, give the spatial heat-flux vector
    !! \[\boldsymbol{q}=-k(T)\boldsymbol{g},\qquad
    !! k(T)=k_0\exp[\beta(T-T_0)].\]
    !! Here \(k_0>0\) is reference conductivity, with dimensions of power
    !! per length per temperature interval. The coefficient \(\beta\) is
    !! finite and has inverse-temperature units. The reference temperature
    !! \(T_0\geq0\) uses the same absolute scale as \(T\). Both temperatures
    !! must be nonnegative, finite and measured from absolute zero.
    !! This exponential dependence is a
    !! *phenomenological ConstituF parametrization*, not a claim that
    !! Fourier originally proposed it or that it fits every material.
    !! It assumes constant logarithmic slope \(\mathrm{d}\ln k/\mathrm{d}T=\beta\) and is
    !! not an Arrhenius dependence on reciprocal absolute temperature.
    !!
    !! The \(3\times4\) exact local Jacobian for packed input
    !! \([\partial_1T,\partial_2T,\partial_3T,T]\) is
    !! \[\frac{\partial q_i}{\partial g_j}=-k\delta_{ij},\qquad
    !! \frac{\partial q_i}{\partial T}=-k\beta g_i.\]
    !! With fixed prescribed temperature/gradient histories, optional
    !! parameter sensitivities are
    !! \[\frac{\partial\boldsymbol{q}}{\partial k_0}
    !! =-\frac{k}{k_0}\boldsymbol{g},\quad
    !! \frac{\partial\boldsymbol{q}}{\partial\beta}
    !! =-k(T-T_0)\boldsymbol{g},\quad
    !! \frac{\partial\boldsymbol{q}}{\partial T_0}
    !! =k\beta\boldsymbol{g}.\]
    !! The down-gradient condition
    !! \[-\boldsymbol{q}\cdot\boldsymbol{g}
    !! =k\|\boldsymbol{g}\|^2\geq0\]
    !! holds. At \(T>0\), dividing this quantity by \(T^2\) gives the
    !! conductive entropy-production rate per unit volume.
    !!
    !! The two inputs are separate schema fields `temperature.gradient`
    !! (three spatial components, temperature interval per length) and
    !! `temperature.absolute` (one material scalar), packed in that order.
    !! Output `heat.flux` has three spatial components and power-per-area units.
    !! Positive components denote transport toward positive coordinate
    !! directions. Outward flux through a surface with outward unit normal \(\boldsymbol{n}\) is \(\boldsymbol{q}\cdot\boldsymbol{n}\).
    !! This is a local constitutive closure, not a heat-equation solver.
    !! It is stateless and has no transient capacity, phase transition,
    !! cross effects, material-axis rotation, or global solution
    !! sensitivities. A nonpositive or nonfinite evaluated conductivity,
    !! nonfinite flux or nonfinite requested derivative rejects the
    !! complete trial. Positive subnormal conductivities are supported.
    !!
    !! Classical conduction origin: J. Fourier, *Théorie analytique de
    !! la chaleur*, Paris, 1822. Standard conduction treatment:
    !! H. S. Carslaw and J. C. Jaeger, *Conduction of Heat in Solids*,
    !! 2nd ed., Clarendon Press, 1959, ISBN 978-0-19-853368-9.
    !! Modern Fourier and entropy-production statements: G. P. Beretta,
    !! MIT 2.43 Advanced Thermodynamics (2024), Lecture 23, slides 23.14
    !! and 23.07:
    !! https://ocw.mit.edu/courses/2-43-advanced-thermodynamics-spring-2024/mit2_43_s24_lec23.pdf .
    !! Code provenance: the stated exponential closure and its analytical
    !! derivatives are independently derived here. No third-party software
    !! code is reused in this revision.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: ndim = 3
    integer, parameter :: ninput = ndim + 1
    integer, parameter :: temperature_component = ninput

    public :: fourier_exponential_t

    type, extends(material_model_t) :: fourier_exponential_t
        !! Stateless isotropic temperature-dependent heat conductor.
        real(wp) :: reference_conductivity = 0.0_wp
            !! Reference conductivity \(k_0>0\), power per length per temperature interval.
        real(wp) :: temperature_coefficient = 0.0_wp
            !! Exponential slope \(\beta\), inverse temperature interval. May be negative.
        real(wp) :: reference_temperature = 0.0_wp
            !! Absolute temperature \(T_0\geq0\), in the chosen absolute-temperature unit.
    contains
        procedure :: describe => describe_fourier_exponential
            !! Declare material fields, parameters and supported operations.
        procedure :: eval_batch => evaluate_fourier_exponential
            !! Evaluate the batched material response, optional tangent and trial state.
        procedure :: eval_sensitivity_batch => sensitivity_fourier_exponential
            !! Evaluate analytical parameter derivatives through the discrete material update.
    end type fourier_exponential_t

contains

    subroutine describe_fourier_exponential(self, spec, status)
        !! Declare the two independent inputs, flux and named scalar parameters.
        class(fourier_exponential_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema populated with fields, physical metadata and model capabilities.
        integer, intent(out) :: status
            !! Operation status code. `cf_ok` indicates success.

        status = cf_invalid_argument
        if (.not. all(ieee_is_finite([self%reference_conductivity, self%temperature_coefficient,self%reference_temperature]))) return
        if (self%reference_conductivity <= 0.0_wp) return
        if (self%reference_temperature < 0.0_wp) return

        call spec%add_input("temperature.gradient", ndim, "vector", "temperature_gradient", "spatial", [ndim])
        call spec%add_input("temperature.absolute", 1, "scalar", "absolute_temperature", "material")
        call spec%add_output("heat.flux", ndim, "vector", "heat_flux", "spatial", [ndim])
        call spec%add_parameter("reference_conductivity")
        call spec%add_parameter("temperature_coefficient")
        call spec%add_parameter("reference_temperature")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_fourier_exponential

    subroutine evaluate_fourier_exponential(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Evaluate \(\boldsymbol{q}(T,\nabla T)\) for every point without state or allocation.
        class(fourier_exponential_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previous packed inputs `(ninput,npoints)` when required by the model, otherwise a zero-sized array.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `temperature.gradient` and `temperature.absolute` inputs `(ninput,npoints)` at the current time step.
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
        call evaluate_points(self,x_new,response,point_status,tangent)
    end subroutine evaluate_fourier_exponential

    subroutine sensitivity_fourier_exponential(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Local derivatives for the selected named material parameters.
        !!
        !! Parameter indices: `reference_conductivity=1`,
        !! `temperature_coefficient=2` and `reference_temperature=3`.
        !! `dsdp(3,nselected,npoints)` follows the requested selection order.
        class(fourier_exponential_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previous packed inputs `(ninput,npoints)` when required by the model, otherwise a zero-sized array.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `temperature.gradient` and `temperature.absolute` inputs `(ninput,npoints)` at the current time step.
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
        call evaluate_points(self,x_new,response,point_status,tangent, selected_parameters,dsdp)
    end subroutine sensitivity_fourier_exponential

    subroutine evaluate_points(self, x_new, response, point_status, tangent, selected_parameters, dsdp)
        !! Evaluate \(k_0\exp(\beta\Delta T)\), with \(\Delta T=T-T_0\), through bounded exponential factors.
        !!
        !! Guard the exponential before evaluating it: no hard-coded
        !! physical conductivity cap or silent underflow to zero.
        class(fourier_exponential_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `temperature.gradient` and `temperature.absolute` inputs `(ninput,npoints)` at the current time step.
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
        real(wp) :: delta_temperature
        real(wp) :: temperature_exponent
        real(wp) :: exponential_factor
        real(wp) :: exponent_limit
        real(wp) :: conductivity
        real(wp) :: conductivity_mantissa
        real(wp) :: beta
        integer :: conductivity_exponent
        integer :: ipoint
        integer :: iparam
        integer :: i

        exponent_limit = log(huge(1.0_wp))-log(nearest(0.0_wp,1.0_wp))
        beta = self%temperature_coefficient

        point_loop: do ipoint = 1, size(x_new,2)
            point_status(ipoint) = cf_invalid_argument
            if (.not. all(ieee_is_finite(x_new(:,ipoint)))) cycle point_loop
            if (x_new(temperature_component,ipoint) < 0.0_wp) cycle point_loop
            delta_temperature = x_new(temperature_component,ipoint) - self%reference_temperature
            if (abs(beta) > 1.0_wp) then
                if (abs(delta_temperature) > huge(1.0_wp)/abs(beta)) cycle point_loop
            end if
            temperature_exponent = beta*delta_temperature
            if (.not. ieee_is_finite(temperature_exponent)) cycle point_loop
            if (abs(temperature_exponent) > exponent_limit) cycle point_loop
            exponential_factor = exp(temperature_exponent/4.0_wp)
            conductivity_mantissa = fraction(self%reference_conductivity)*fraction(exponential_factor)**4
            conductivity_exponent = exponent(self%reference_conductivity)+4*exponent(exponential_factor)
            conductivity = scale(conductivity_mantissa,conductivity_exponent)
            if (.not. ieee_is_finite(conductivity)) cycle point_loop
            if (conductivity <= 0.0_wp) cycle point_loop
            response(:,ipoint) = scaled_conductivity_product(-conductivity_mantissa,conductivity_exponent,1.0_wp,x_new(1:ndim,ipoint))
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop

            if (present(tangent)) then
                tangent(:,:,ipoint) = 0.0_wp
                do i = 1, ndim
                    tangent(i,i,ipoint) = -conductivity
                end do
                tangent(:,temperature_component,ipoint) = scaled_conductivity_product(-conductivity_mantissa,conductivity_exponent,beta,x_new(1:ndim,ipoint))
                if (.not. all(ieee_is_finite(tangent(:,:,ipoint)))) cycle point_loop
            end if

            if (present(selected_parameters)) then
                parameter_loop: do iparam = 1, size(selected_parameters)
                    select case (selected_parameters(iparam))
                    case (1)
                        dsdp(:,iparam,ipoint) = scaled_reference_derivative(exponential_factor,x_new(1:ndim,ipoint))
                    case (2)
                        dsdp(:,iparam,ipoint) = scaled_conductivity_product(-conductivity_mantissa,conductivity_exponent,delta_temperature,x_new(1:ndim,ipoint))
                    case (3)
                        dsdp(:,iparam,ipoint) = scaled_conductivity_product(conductivity_mantissa,conductivity_exponent,beta,x_new(1:ndim,ipoint))
                    case default
                        exit parameter_loop
                    end select
                end do parameter_loop
                if (iparam <= size(selected_parameters)) cycle point_loop
                if (.not. all(ieee_is_finite(dsdp(:,:,ipoint)))) cycle point_loop
            end if
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_points

    pure elemental function scaled_conductivity_product(conductivity_mantissa, conductivity_exponent, factor, gradient) result(derivative)
        !! Combine the unrounded conductivity representation with two factors
        !! before rounding the final flux or derivative to the working kind.
        real(wp), intent(in) :: conductivity_mantissa
            !! Signed conductivity mantissa before restoring its magnitude.
        integer, intent(in) :: conductivity_exponent
            !! Radix exponent paired with the conductivity mantissa.
        real(wp), intent(in) :: factor
            !! Unit factor, temperature coefficient or prescribed temperature difference.
        real(wp), intent(in) :: gradient
            !! One component of the prescribed temperature gradient.
        real(wp) :: derivative

        derivative = scale((conductivity_mantissa*fraction(factor))*fraction(gradient), conductivity_exponent+exponent(factor)+exponent(gradient))
    end function scaled_conductivity_product

    pure elemental function scaled_reference_derivative(exponential_factor, gradient) result(derivative)
        !! Form \(-\exp[\beta(T-T_0)]g_i\) before rounding the conductivity.
        real(wp), intent(in) :: exponential_factor
            !! Bounded quarter-exponential \(\exp[\beta(T-T_0)/4]\).
        real(wp), intent(in) :: gradient
            !! One component of the prescribed temperature gradient.
        real(wp) :: derivative

        derivative = scale(-(fraction(exponential_factor)**4)*fraction(gradient), 4*exponent(exponential_factor)+exponent(gradient))
    end function scaled_reference_derivative

end module constituf_fourier_exponential

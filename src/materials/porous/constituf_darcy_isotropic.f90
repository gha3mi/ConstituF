! SPDX-License-Identifier: Apache-2.0

module constituf_darcy_isotropic
    !! Three-dimensional isotropic saturated single-phase Darcy flow.
    !!
    !! For a prescribed spatial pressure gradient
    !! \(\boldsymbol{g}=\nabla p\) (pressure per length), Darcy's superficial (specific
    !! discharge) velocity \(\boldsymbol{v}\) (length per time) is
    !! \[\boldsymbol{v}=-\frac{k}{\mu}\boldsymbol{g}.\]
    !! The intrinsic permeability \(k>0\) has units of length squared and the
    !! dynamic viscosity \(\mu>0\) has units of stress times time. The mobility \(k/\mu\)
    !! must be finite and strictly positive in the working precision at
    !! model initialization.
    !! The exact local constitutive Jacobian and optional derivatives at
    !! fixed prescribed pressure gradient are
    !! \[\frac{\partial v_i}{\partial g_j}=-\frac{k}{\mu}\delta_{ij},\]
    !! \[\frac{\partial\boldsymbol{v}}{\partial k}
    !! =-\frac{\boldsymbol{g}}{\mu},\qquad
    !! \frac{\partial\boldsymbol{v}}{\partial\mu}
    !! =\frac{k}{\mu^2}\boldsymbol{g}.\]
    !! The pressure-work dissipation rate per unit bulk volume (power per volume) is
    !! \[-\boldsymbol{v}\cdot\boldsymbol{g}
    !! =\frac{k}{\mu}\|\boldsymbol{g}\|^2\geq0.\]
    !! This derived rate is not a returned response field.
    !!
    !! Input `pressure.gradient` is the spatial pressure gradient (pressure per length).
    !! Output `fluid.darcy_velocity` is the spatial specific discharge
    !! (length per time), *not* the interstitial pore velocity, for a stationary porous
    !! skeleton. Both have three Cartesian components in the same fixed
    !! frame. Gravity is omitted for the stated pressure-gradient input.
    !! A caller including a body acceleration \(\boldsymbol b\) may supply
    !! \(\boldsymbol g_{\mathrm{eff}}=\nabla p-\rho\boldsymbol b\)
    !! externally in the same field and units. Both the flux law and
    !! dissipation then use this effective driving gradient. Capillary
    !! pressure, relative permeability, saturation changes, thermal
    !! variation, porosity, compressibility and fluid inertia are not
    !! modeled. This is only a local constitutive law, not a mass
    !! conservation or groundwater solver.
    !!
    !! The material has no integration-point history, `dt` dependence, or
    !! parameter-dependent prescribed inputs. Failed nonfinite inputs,
    !! fluxes or requested sensitivities invalidate the complete trial.
    !! No heap allocations, I/O or name lookups occur in the batch loop.
    !!
    !! Original experimental source (not code): H. Darcy,
    !! *Les fontaines publiques de la ville de Dijon*, V. Dalmont,
    !! Paris, 1856, Appendix D. Historical discussion:
    !! G. O. Brown, "Henry Darcy and the making of a law",
    !! Water Resources Research 38(7) (2002), DOI 10.1029/2001WR000727.
    !! The 3D pressure-gradient law is a modern vector generalization.
    !! Its pressure form, gravity correction and specific-discharge
    !! convention are stated in the COMSOL Porous Media Flow Module
    !! User's Guide 6.1 (2022), "Darcy's Law--Equation Formulation",
    !! Eqs. (4-1) and (4-3), pp. 114-116:
    !! https://doc.comsol.com/6.1/doc/com.comsol.help.porous/PorousMediaFlowModuleUsersGuide.pdf .
    !! No third-party software code is reused in this revision.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: ndim = 3

    public :: darcy_isotropic_t

    type, extends(material_model_t) :: darcy_isotropic_t
        !! Construct as `darcy_isotropic_t(permeability,dynamic_viscosity)`.
        !! The caller's `material_t` clones both constant scalar parameters.
        real(wp) :: permeability = 0.0_wp
            !! Positive intrinsic permeability \(k\) in squared-length units.
        real(wp) :: dynamic_viscosity = 0.0_wp
            !! Positive fluid viscosity \(\mu\) in stress times time.
    contains
        procedure :: describe => describe_darcy
            !! Declare material fields, parameters and supported operations.
        procedure :: eval_batch => evaluate_darcy
            !! Evaluate the batched material response, optional tangent and trial state.
        procedure :: eval_sensitivity_batch => sensitivity_darcy
            !! Evaluate analytical parameter derivatives through the discrete material update.
    end type darcy_isotropic_t

contains

    subroutine describe_darcy(self, spec, status)
        !! Declare pressure-gradient input, superficial velocity and parameters.
        !! The tangent and both named parameter sensitivities are exact.
        class(darcy_isotropic_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema populated with fields, physical metadata and model capabilities.
        integer, intent(out) :: status
            !! Operation status code. `cf_ok` indicates success.
        real(wp) :: mobility

        status = cf_invalid_argument
        if (.not. all(ieee_is_finite([self%permeability,self%dynamic_viscosity]))) return
        if (self%permeability <= 0.0_wp) return
        if (self%dynamic_viscosity <= 0.0_wp) return
        mobility = self%permeability/self%dynamic_viscosity
        if (.not. ieee_is_finite(mobility)) return
        if (mobility <= 0.0_wp) return

        call spec%add_input("pressure.gradient", ndim, "vector", "pressure_gradient", "spatial", [ndim])
        call spec%add_output("fluid.darcy_velocity", ndim, "vector", "specific_discharge", "spatial", [ndim])
        call spec%add_parameter("permeability")
        call spec%add_parameter("dynamic_viscosity")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_darcy

    subroutine evaluate_darcy(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Evaluate the stateless 3D velocity response and optional Jacobian.
        class(darcy_isotropic_t), intent(in) :: self
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
    end subroutine evaluate_darcy

    subroutine sensitivity_darcy(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Local parameter derivatives in selected order. No derivative history.
        !! Parameter IDs are `permeability=1` and `dynamic_viscosity=2`.
        class(darcy_isotropic_t), intent(in) :: self
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
    end subroutine sensitivity_darcy

    subroutine evaluate_points(self, gradient, response, point_status, tangent, selected_parameters, dsdp)
        !! One batch dispatch, finite values and no point-local allocation.
        class(darcy_isotropic_t), intent(in) :: self
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
        real(wp) :: mobility
        real(wp) :: mobility_fraction
        integer :: mobility_exponent
        integer :: ipoint
        integer :: i
        integer :: iparam

        point_status = cf_invalid_argument
        if (present(selected_parameters)) then
            if (any(selected_parameters < 1) .or. any(selected_parameters > 2)) return
        end if
        mobility = self%permeability/self%dynamic_viscosity
        mobility_fraction = fraction(self%permeability)/fraction(self%dynamic_viscosity)
        mobility_exponent = exponent(self%permeability)-exponent(self%dynamic_viscosity)

        point_loop: do ipoint = 1, size(gradient,2)
            if (.not. all(ieee_is_finite(gradient(:,ipoint)))) cycle point_loop
            if (mobility < tiny(1.0_wp)) then
                ! Retain the unrounded mobility until the final flux product.
                response(:,ipoint) = -scale(mobility_fraction*fraction(gradient(:,ipoint)),mobility_exponent+exponent(gradient(:,ipoint)))
            else
                response(:,ipoint) = -mobility*gradient(:,ipoint)
            end if
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop

            if (present(tangent)) then
                tangent(:,:,ipoint) = 0.0_wp
                do i = 1, ndim
                    tangent(i,i,ipoint) = -mobility
                end do
            end if

            if (present(selected_parameters)) then
                do iparam = 1, size(selected_parameters)
                    select case (selected_parameters(iparam))
                    case (1)
                        dsdp(:,iparam,ipoint) = -gradient(:,ipoint)/self%dynamic_viscosity
                    case (2)
                        dsdp(:,iparam,ipoint) = scaled_viscosity_derivative(self%permeability,self%dynamic_viscosity,gradient(:,ipoint))
                    case default
                        cycle point_loop
                    end select
                end do
                if (.not. all(ieee_is_finite(dsdp(:,:,ipoint)))) cycle point_loop
            end if
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_points

    pure elemental function scaled_viscosity_derivative(permeability, dynamic_viscosity, gradient) result(value)
        !! Evaluate \(k g_i/\mu^2\) without losing a finite derivative when velocity underflows.
        real(wp), intent(in) :: permeability
            !! Positive intrinsic permeability in squared-length units.
        real(wp), intent(in) :: dynamic_viscosity
            !! Positive dynamic viscosity in stress times time.
        real(wp), intent(in) :: gradient
            !! One prescribed pressure-gradient component in pressure per length.
        real(wp) :: value

        value = scale((fraction(permeability)*fraction(gradient))/fraction(dynamic_viscosity)**2, exponent(permeability)+exponent(gradient)-2*exponent(dynamic_viscosity))
    end function scaled_viscosity_derivative

end module constituf_darcy_isotropic

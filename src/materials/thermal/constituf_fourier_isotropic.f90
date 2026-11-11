! SPDX-License-Identifier: Apache-2.0

module constituf_fourier_isotropic
    !! Homogeneous, isotropic Fourier heat conduction in three spatial dimensions.
    !!
    !! For a prescribed spatial temperature gradient \(\boldsymbol{g}=\nabla T\),
    !! the spatial heat-flux vector \(\boldsymbol{q}\) is
    !! \[\boldsymbol{q}=-k\boldsymbol{g},\qquad k>0.\]
    !! The exact constitutive Jacobian is
    !! \[\frac{\partial\boldsymbol{q}}{\partial\boldsymbol{g}}
    !! =-k\boldsymbol{I},\]
    !! and the optional local analytical parameter derivative, at fixed
    !! temperature gradient, is
    !! \[\frac{\partial\boldsymbol{q}}{\partial k}=-\boldsymbol{g}.\]
    !! The down-gradient condition
    !! \[-\boldsymbol{q}\cdot\nabla T=k\|\nabla T\|^2\geq0\]
    !! holds for every gradient. At absolute \(T>0\), conductive entropy
    !! production is this expression divided by \(T^2\). This model does
    !! not take absolute temperature as an input.
    !!
    !! Input `temperature.gradient` is a spatial vector `(3,npoints)` in temperature interval per length.
    !! Output `heat.flux` is a spatial vector `(3,npoints)` in power per area.
    !! The conductivity \(k\) has units of power per length per temperature interval.
    !! There is no internal history,
    !! time-step dependence, heat capacity, energy balance or temperature
    !! dependence. This law alone does not solve the heat equation. Use a
    !! single coherent unit system. Positive flux denotes heat transported
    !! towards the positive coordinate direction.
    !!
    !! Nonfinite inputs, nonfinite responses and invalid conductivity
    !! report failure status. No allocations occur in the integration loop.
    !!
    !! Original reference: J. Fourier, *Théorie analytique de la chaleur*,
    !! Firmin Didot, Paris, 1822. English translation: *The Analytical
    !! Theory of Heat*, translated by Alexander Freeman, Cambridge
    !! University Press, 1878:
    !! https://www.cambridge.org/core/books/analytical-theory-of-heat/F6D4802336FABD1116DDA4AA3FE6EFAA .
    !! Modern formulation and conductive entropy production: G. P. Beretta,
    !! MIT 2.43 Advanced Thermodynamics (2024), Lecture 23, slides 23.14
    !! and 23.07:
    !! https://ocw.mit.edu/courses/2-43-advanced-thermodynamics-spring-2024/mit2_43_s24_lec23.pdf .
    !! Code provenance: the flux and its analytical derivatives directly
    !! implement the stated Fourier law. No third-party software code is
    !! reused in this revision.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: ndim = 3

    public :: fourier_isotropic_t

    type, extends(material_model_t) :: fourier_isotropic_t
        !! Spatial isotropic conductive material with constant conductivity.
        !!
        !! Construct as `fourier_isotropic_t(conductivity)` and pass to
        !! `material_t%init`. The material handle clones the coefficient.
        !! No thermomechanical coupling, finite-rate or gradient regularization.
        real(wp) :: conductivity = 0.0_wp
            !! Conductivity \(k\) in power per length per temperature interval, finite and strictly positive.
    contains
        procedure :: describe => describe_fourier
            !! Declare material fields, parameters and supported operations.
        procedure :: eval_batch => evaluate_fourier
            !! Evaluate the batched material response, optional tangent and trial state.
        procedure :: eval_sensitivity_batch => sensitivity_fourier
            !! Evaluate analytical parameter derivatives through the discrete material update.
    end type fourier_isotropic_t

contains

    subroutine describe_fourier(self, spec, status)
        !! Declare a 3D temperature-gradient input and heat-flux response.
        !!
        !! Parameter `conductivity` is scalar, power per length per temperature interval. Stateless, without
        !! required `dt`. Exact tangent and parameter sensitivity supported.
        class(fourier_isotropic_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema populated with fields, physical metadata and model capabilities.
        integer, intent(out) :: status
            !! Operation status code. `cf_ok` indicates success.

        status = cf_invalid_argument
        if (.not. ieee_is_finite(self%conductivity)) return
        if (self%conductivity <= 0.0_wp) return

        call spec%add_input("temperature.gradient", ndim, "vector", "temperature_gradient", "spatial", [ndim])
        call spec%add_output("heat.flux", ndim, "vector", "heat_flux", "spatial", [ndim])
        call spec%add_parameter("conductivity")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_fourier

    subroutine evaluate_fourier(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Evaluate \(\boldsymbol{q}=-k\nabla T\) without persistent state.
        !!
        !! The caller owns `x_new(ndim,npoints)`, `response(ndim,npoints)`
        !! and optional `tangent(ndim,ndim,npoints)`. A failed point makes the
        !! entire trial invalid. Response buffers are unspecified on failure.
        class(fourier_isotropic_t), intent(in) :: self
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
        !! Exact derivative with respect to the selected conductivity parameter.
        !!
        !! `dsdp(ndim,1,npoints)` has temperature-gradient units and the
        !! zero-component sensitivity state remains unchanged.
        class(fourier_isotropic_t), intent(in) :: self
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

    subroutine evaluate_points(self, x_new, response, point_status, tangent, selected_parameters, dsdp)
        !! Fused flux, exact input Jacobian and optional selected derivatives.
        !! Each point is independent. Nonfinite inputs or outputs fail locally.
        class(fourier_isotropic_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `temperature.gradient` inputs `(ninput,npoints)` at the current time step.
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
        integer :: i
        integer :: j
        integer :: ipoint

        point_loop: do ipoint = 1, size(x_new,2)
            point_status(ipoint) = cf_invalid_argument
            if (.not. all(ieee_is_finite(x_new(:,ipoint)))) cycle point_loop
            response(:,ipoint) = -self%conductivity*x_new(:,ipoint)
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop

            if (present(tangent)) then
                tangent(:,:,ipoint) = 0.0_wp
                do i = 1, ndim
                    tangent(i,i,ipoint) = -self%conductivity
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

end module constituf_fourier_isotropic

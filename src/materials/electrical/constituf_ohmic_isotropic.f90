! SPDX-License-Identifier: Apache-2.0

module constituf_ohmic_isotropic
    !! Homogeneous isotropic Ohmic electrical conduction in 3D.
    !!
    !! For a prescribed spatial electric field \(\boldsymbol{E}\) (electric potential per length),
    !! the conventional conduction current density \(\boldsymbol{J}\) (current per area)
    !! is \[\boldsymbol{J}=\sigma\boldsymbol{E},\quad
    !! \sigma>0.\]
    !! Conductivity has units of current density divided by electric field.
    !! The exact input Jacobian and local scalar parameter sensitivity are
    !! \[\frac{\partial J_i}{\partial E_j}=\sigma\delta_{ij},\quad
    !! \frac{\partial\boldsymbol{J}}{\partial\sigma}=\boldsymbol{E}.\]
    !! Instantaneous Joule dissipation per volume (not an output field) satisfies
    !! \[\boldsymbol{E}\cdot\boldsymbol{J}
    !! =\sigma\|\boldsymbol{E}\|^2\geq0.\]
    !! The external \(\boldsymbol{E}\) and \(\boldsymbol{J}\) are spatial Cartesian vectors ordered \([1,2,3]\).
    !! The medium is stationary in this frame. Motional electromotive
    !! terms and convective charge transport are excluded.
    !! \(\boldsymbol{E}\) is the physical field, not a potential gradient needing a sign
    !! change. This model is stateless, time-step-independent and isothermal.
    !! Local sensitivities hold \(\boldsymbol{E}\) fixed and do not include field-solver effects.
    !! No displacement current, dielectric polarization, Hall current,
    !! temperature feedback, nonlinear conduction or global conservation
    !! equation is included. Nonfinite field/response invalidates a trial.
    !! The point loop performs no allocation, file I/O or field-name lookup.
    !!
    !! Scientific origin: G. S. Ohm, *Die galvanische Kette, mathematisch
    !! bearbeitet*, T. H. Riemann, Berlin (1827). Digitized edition:
    !! DOI 10.5479/sil.354716.39088005838644.
    !! Direct continuum form: S. Hughes, MIT 8.022 (2005), Lecture 7,
    !! section 7.4, pp. 66--67:
    !! https://web.mit.edu/sahughes/www/8.022/lec07.pdf .
    !! Joule power density: MIT, *Electromagnetic Fields and Energy*,
    !! section 11.3:
    !! https://web.mit.edu/6.013_book/www/chapter11/11.3.html .
    !! Code provenance: this continuum specialization and its
    !! analytical derivatives directly implement the stated law.
    !! No third-party software code is reused in this revision.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: ndim = 3

    public :: ohmic_isotropic_t

    type, extends(material_model_t) :: ohmic_isotropic_t
        !! Construct as `ohmic_isotropic_t(conductivity)` in current density per electric field.
        real(wp) :: conductivity = 0.0_wp
            !! Finite, positive absolute conductivity \(\sigma\), current density per electric field.
    contains
        procedure :: describe => describe_ohmic
            !! Declare input and output fields, parameter names and constitutive capabilities.
        procedure :: eval_batch => evaluate_ohmic
            !! Evaluate the batched constitutive response and optional consistent tangent.
        procedure :: eval_sensitivity_batch => sensitivity_ohmic
            !! Evaluate exact scalar-parameter derivatives of the constitutive response.
    end type ohmic_isotropic_t

contains

    subroutine describe_ohmic(self, spec, status)
        !! Declare electric-field/current-density spatial 3-vectors.
        class(ohmic_isotropic_t), intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        type(material_spec_t), intent(out) :: spec
            !! Packed field schema, material parameter metadata and declared capabilities.
        integer, intent(out) :: status
            !! Operation status; `cf_ok` indicates success.

        status = cf_invalid_argument
        if (.not. ieee_is_finite(self%conductivity)) return
        if (self%conductivity <= 0.0_wp) return

        call spec%add_input("electric.field", ndim, "vector", "electric_field", "spatial", [ndim])
        call spec%add_output("electric.current_density", ndim, "vector", "electric_current_density", "spatial", [ndim])
        call spec%add_parameter("conductivity")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_ohmic

    subroutine evaluate_ohmic(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Compute local \(\boldsymbol{J}(\boldsymbol{E})\) with no history or time-step dependence.
        class(ohmic_isotropic_t), intent(in) :: self
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
            !! Packed `electric.current_density` output `(noutput,npoints)`.
        real(wp), intent(in) :: dt
            !! Prescribed time increment; required to be positive only for time-dependent models.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic response Jacobian `(noutput,ninput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Status of each point; a successful update must set its code to `cf_ok`.

        state_trial = state_old
        call evaluate_points(self, x_new, response, point_status, tangent)
    end subroutine evaluate_ohmic

    subroutine sensitivity_ohmic(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Local derivative of \(\boldsymbol{J}\) with respect to selected absolute conductivity.
        !! Parameter ID is `1`. `dsdp` has shape `(3,nselected,npoints)`.
        class(ohmic_isotropic_t), intent(in) :: self
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
            !! Packed `electric.current_density` output `(noutput,npoints)`.
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
        call evaluate_points(self, x_new, response, point_status, tangent, selected_parameters, dsdp)
    end subroutine sensitivity_ohmic

    subroutine evaluate_points(self, x_new, response, point_status, tangent, selected_parameters, dsdp)
        !! Fused allocation-free batch path. Unrequested derivatives are skipped.
        class(ohmic_isotropic_t), intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `electric.field` input `(ninput,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `electric.current_density` output `(noutput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Status of each point; a successful update must set its code to `cf_ok`.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic response Jacobian `(noutput,ninput,npoints)`.
        integer, intent(in), optional :: selected_parameters(:)
            !! One-based IDs of scalar parameters in the requested sensitivity-axis order.
        real(wp), intent(out), optional :: dsdp(:,:,:)
            !! Exact discrete response derivatives `(noutput,nselected,npoints)` at fixed input histories.
        integer :: ipoint
        integer :: i
        integer :: iparam

        point_loop: do ipoint = 1, size(x_new,2)
            point_status(ipoint) = cf_invalid_argument
            if (.not. all(ieee_is_finite(x_new(:,ipoint)))) cycle point_loop
            response(:,ipoint) = self%conductivity*x_new(:,ipoint)
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop

            if (present(tangent)) then
                tangent(:,:,ipoint) = 0.0_wp
                do i = 1, ndim
                    tangent(i,i,ipoint) = self%conductivity
                end do
            end if

            if (present(selected_parameters)) then
                parameter_loop: do iparam = 1, size(selected_parameters)
                    select case (selected_parameters(iparam))
                    case (1)
                        dsdp(:,iparam,ipoint) = x_new(:,ipoint)
                    case default
                        exit parameter_loop
                    end select
                end do parameter_loop
                if (iparam <= size(selected_parameters)) cycle point_loop
            end if
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_points

end module constituf_ohmic_isotropic

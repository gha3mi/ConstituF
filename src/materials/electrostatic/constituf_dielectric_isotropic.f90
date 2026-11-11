! SPDX-License-Identifier: Apache-2.0

module constituf_dielectric_isotropic
    !! Linear isotropic, homogeneous electrostatic dielectric in 3D.
    !!
    !! The prescribed spatial electric field \(\boldsymbol{E}\) (electric potential per length)
    !! produces electric displacement \(\boldsymbol{D}\) (charge per area) via
    !! \[\boldsymbol{D}=\varepsilon\boldsymbol{E},\qquad
    !! \varepsilon>0.\]
    !! Permittivity has units of electric displacement divided by electric field.
    !! The parameter is the *absolute* permittivity, including the vacuum
    !! contribution: \(\boldsymbol D=\varepsilon_0\boldsymbol E+\boldsymbol P\)
    !! and \(\varepsilon=\varepsilon_0\varepsilon_r\), where \(\boldsymbol{P}\) is the
    !! material polarization density.
    !! It is neither relative permittivity nor dielectric susceptibility.
    !! The supplied absolute coefficient is used directly, with no extra
    !! vacuum-permittivity multiplier applied by the implementation.
    !! The exact local \(3\times3\) tangent and selected scalar sensitivity are
    !! \[\frac{\partial D_i}{\partial E_j}=
    !! \varepsilon\delta_{ij},\qquad
    !! \frac{\partial\boldsymbol{D}}{\partial\varepsilon}=
    !! \boldsymbol{E}.\]
    !! The stored electrostatic energy density is
    !! \(w=\tfrac12\boldsymbol{E}\cdot\boldsymbol{D}\geq0\) (energy per volume).
    !! The model does not output that energy. The sign convention for
    !! \(\boldsymbol{D}\) is aligned with \(\boldsymbol{E}\). No minus sign is applied to the field itself.
    !! Both field and displacement use the spatial Cartesian frame.
    !!
    !! This is a stateless, rate-independent, reversible constitutive
    !! relation. It does not solve electrostatic field equations, provide
    !! nonlinear polarization, hysteresis, conductivity, Maxwell stress,
    !! temperature dependence or electromechanical coupling.
    !! Invalid coefficients are rejected at initialization. Nonfinite input,
    !! response or requested derivatives invalidate the trial by point status.
    !! No dynamic allocations occur in the hot loop.
    !!
    !! Scientific basis (not implementation code): J. D. Jackson,
    !! *Classical Electrodynamics*, 3rd ed., Wiley, 1998,
    !! ISBN 978-0-471-30932-1, macroscopic dielectrics in Ch. 4.
    !! Direct permittivity convention: MIT, *Electromagnetic Fields and
    !! Energy*, section 6.4, equations (3)--(4):
    !! https://web.mit.edu/6.013_book/www/chapter6/6.4.html .
    !! Reversible field energy: section 11.4, "Energy Densities" and
    !! Example 11.4.1:
    !! https://web.mit.edu/6.013_book/www/chapter11/11.4.html .
    !! Code provenance: the displacement and its analytical derivatives
    !! directly implement the stated constitutive relation. No third-party
    !! software code is reused in this revision.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: ndim = 3

    public :: dielectric_isotropic_t

    type, extends(material_model_t) :: dielectric_isotropic_t
        !! Construct as `dielectric_isotropic_t(permittivity)` in electric displacement per electric field.
        real(wp) :: permittivity = 0.0_wp
            !! Finite, positive absolute permittivity \(\varepsilon\), electric displacement per electric field.
    contains
        procedure :: describe => describe_dielectric
            !! Declare input and output fields, parameter names and constitutive capabilities.
        procedure :: eval_batch => evaluate_dielectric
            !! Evaluate the batched constitutive response and optional consistent tangent.
        procedure :: eval_sensitivity_batch => sensitivity_dielectric
            !! Evaluate exact scalar-parameter derivatives of the constitutive response.
    end type dielectric_isotropic_t

contains

    subroutine describe_dielectric(self, spec, status)
        !! Declare electric-field/displacement spatial 3-vectors.
        class(dielectric_isotropic_t), intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        type(material_spec_t), intent(out) :: spec
            !! Packed field schema, material parameter metadata and declared capabilities.
        integer, intent(out) :: status
            !! Operation status; `cf_ok` indicates success.

        status = cf_invalid_argument
        if (.not. ieee_is_finite(self%permittivity)) return
        if (self%permittivity <= 0.0_wp) return

        call spec%add_input("electric.field", ndim, "vector", "electric_field", "spatial", [ndim])
        call spec%add_output("electric.displacement", ndim, "vector", "electric_displacement", "spatial", [ndim])
        call spec%add_parameter("permittivity")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_dielectric

    subroutine evaluate_dielectric(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Compute local \(\boldsymbol{D}(\boldsymbol{E})\) with no history or time-step dependence.
        class(dielectric_isotropic_t), intent(in) :: self
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
        call evaluate_points(self, x_new, response, point_status, tangent)
    end subroutine evaluate_dielectric

    subroutine sensitivity_dielectric(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Local derivative of \(\boldsymbol{D}\) with respect to selected absolute permittivity.
        !! Parameter ID is `1`. `dsdp` has shape `(3,nselected,npoints)`.
        class(dielectric_isotropic_t), intent(in) :: self
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
        call evaluate_points(self, x_new, response, point_status, tangent, selected_parameters, dsdp)
    end subroutine sensitivity_dielectric

    subroutine evaluate_points(self, x_new, response, point_status, tangent, selected_parameters, dsdp)
        !! Fused allocation-free batch path. Unrequested derivatives are skipped.
        class(dielectric_isotropic_t), intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `electric.field` input `(ninput,npoints)`.
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
        integer :: ipoint
        integer :: i
        integer :: iparam

        point_loop: do ipoint = 1, size(x_new,2)
            point_status(ipoint) = cf_invalid_argument
            if (.not. all(ieee_is_finite(x_new(:,ipoint)))) cycle point_loop
            response(:,ipoint) = self%permittivity*x_new(:,ipoint)
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop

            if (present(tangent)) then
                tangent(:,:,ipoint) = 0.0_wp
                do i = 1, ndim
                    tangent(i,i,ipoint) = self%permittivity
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

end module constituf_dielectric_isotropic

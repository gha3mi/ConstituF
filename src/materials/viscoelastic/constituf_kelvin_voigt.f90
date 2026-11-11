! SPDX-License-Identifier: Apache-2.0

module constituf_kelvin_voigt
    !! Three-dimensional, infinitesimal-strain, isotropic Kelvin--Voigt solid.
    !!
    !! A linear spring and a Newtonian dashpot act in parallel. Cauchy stress
    !! at prescribed symmetric infinitesimal strain and strain rate is
    !! \[\boldsymbol{\sigma}
    !! =K\mathrm{tr}(\boldsymbol{\varepsilon})\boldsymbol{I}
    !! +2G\mathrm{dev}\boldsymbol{\varepsilon}
    !! +\zeta\mathrm{tr}(\dot{\boldsymbol{\varepsilon}})\boldsymbol{I}
    !! +2\eta\mathrm{dev}\dot{\boldsymbol{\varepsilon}}.\]
    !! \(K,G>0\) are bulk and shear moduli in the chosen stress unit.
    !! The bulk and shear viscosities \(\zeta,\eta\geq0\) have dimensions
    !! of stress times time, with the same time unit as `dt`. The viscous stress is
    !! dissipative because its local power is
    !! \[\boldsymbol{\sigma}_{\mathrm{visc}}:
    !! \dot{\boldsymbol{\varepsilon}}
    !! =\zeta[\mathrm{tr}(\dot{\boldsymbol{\varepsilon}})]^2
    !! +2\eta\|\mathrm{dev}\dot{\boldsymbol{\varepsilon}}\|^2\geq0.\]
    !! Zero viscosities recover linear isotropic elasticity.
    !!
    !! ConstituF owns the *previous committed input strain*, not a separate
    !! Kelvin--Voigt internal state. The initial strain is explicitly
    !! required at `init` and a positive increment \(\Delta t\) is mandatory
    !! at each `eval`. The discrete rate is the backward difference:
    !! \[\dot{\boldsymbol{\varepsilon}}_{n+1}
    !! =(\boldsymbol{\varepsilon}_{n+1}
    !! -\boldsymbol{\varepsilon}_{n})/\Delta t.\]
    !! This determines a discrete constitutive response and an exact
    !! algorithmic tangent at fixed committed previous strain:
    !! \[\mathbb{C}_{\mathrm{alg}}=
    !! (K+\zeta/\Delta t)\boldsymbol{I}\otimes\boldsymbol{I}
    !! +2(G+\eta/\Delta t)\mathbb{I}_{\mathrm{dev}}.\]
    !! For fixed prescribed strain and time-increment histories, exact
    !! local analytical parameter sensitivities are
    !! \[\frac{\partial\sigma}{\partial K}
    !! =\mathrm{tr}(\varepsilon_{n+1})I,\quad
    !! \frac{\partial\sigma}{\partial G}=2\mathrm{dev}\varepsilon_{n+1},\]
    !! \[\frac{\partial\sigma}{\partial\zeta}
    !! =\mathrm{tr}(\Delta\varepsilon)I/\Delta t,\quad
    !! \frac{\partial\sigma}{\partial\eta}
    !! =2\mathrm{dev}(\Delta\varepsilon)/\Delta t.\]
    !! The prescribed input history is parameter-independent, so no
    !! additional derivative history is needed. This is not a global FEM
    !! equilibrium sensitivity.
    !!
    !! Input `strain.small` and output `stress.cauchy` are symmetric
    !! 3D spatial tensors. Each such tensor \(\boldsymbol{A}\) is packed as
    !! \[\operatorname{mandel}(\boldsymbol{A})=
    !! [A_{11},A_{22},A_{33},\sqrt{2}A_{23},
    !! \sqrt{2}A_{13},\sqrt{2}A_{12}]^{\mathsf T}.\]
    !! Strain is dimensionless. Cauchy stress uses the same unit as the
    !! elastic moduli. Tension is positive.
    !! The material has no thermal dependence, aging, volumetric creep
    !! internal variables, or finite-strain objective stress rates.
    !! Holding strain fixed after an accepted ramp removes viscous stress
    !! in the following time step: Kelvin--Voigt is **not** a stress-
    !! relaxation model with exponentially decaying internal overstress.
    !! Use the Maxwell model for that different constitutive behavior.
    !!
    !! References: R. M. Christensen, *Theory of Viscoelasticity:
    !! An Introduction*, second edition, Academic Press, 1982,
    !! ISBN 978-0-12-174252-2, Section 1.5.
    !! J. Linn, H. Lang and A. Tuganov, "Derivation of a Viscoelastic
    !! Constitutive Model of Kelvin--Voigt Type for Cosserat Rods",
    !! Fraunhofer ITWM Report 225 (2013), Section 4, Eq. (29),
    !! gives the underlying three-dimensional isotropic viscosity tensor:
    !! https://publica-rest.fraunhofer.de/server/api/core/bitstreams/22d99f07-6ca4-4770-a5e1-89dcb9dffb07/content.
    !! The parallel spring/dashpot model here is the classical single
    !! Kelvin--Voigt element, not a generalized series arrangement.
    !! The backward-difference response and analytical derivatives are
    !! derived directly from the stated equations. This revision reuses
    !! no third-party source code.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: ndim = 3
    integer, parameter :: nmandel = 6

    public :: kelvin_voigt_t

    type, extends(material_model_t) :: kelvin_voigt_t
        !! Parallel isotropic spring/dashpot with fixed prescribed history.
        !!
        !! Construct `kelvin_voigt_t(K,G,zeta,eta)`, pass to `material_t%init`
        !! with `initial_input(nmandel,npoints)`, then `eval` with positive `dt`.
        !! The handle owns and commits the previous strain for the rate.
        real(wp) :: bulk_modulus = 0.0_wp
            !! Bulk elastic modulus \(K>0\), in the chosen stress unit.
        real(wp) :: shear_modulus = 0.0_wp
            !! Shear elastic modulus \(G>0\), in the chosen stress unit.
        real(wp) :: bulk_viscosity = 0.0_wp
            !! Nonnegative bulk viscosity \(\zeta\), stress times time.
        real(wp) :: shear_viscosity = 0.0_wp
            !! Nonnegative shear viscosity \(\eta\), stress times time.
    contains
        procedure :: describe => describe_kelvin_voigt
            !! Declare material fields, parameters and supported operations.
        procedure :: eval_batch => evaluate_kelvin_voigt
            !! Evaluate the batched material response, optional tangent and trial state.
        procedure :: eval_sensitivity_batch => sensitivity_kelvin_voigt
            !! Evaluate analytical parameter derivatives through the discrete material update.
    end type kelvin_voigt_t

contains

    subroutine describe_kelvin_voigt(self, spec, status)
        !! Validate passive elastic/viscous coefficients and field metadata.
        !!
        !! Previous input and initial input are both required. The material
        !! declares no separate state components or zero-initial-state hook.
        class(kelvin_voigt_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema populated with fields, physical metadata and model capabilities.
        integer, intent(out) :: status
            !! Operation status code. `cf_ok` indicates success.

        status = cf_invalid_argument
        if (.not. all(ieee_is_finite([self%bulk_modulus,self%shear_modulus, self%bulk_viscosity,self%shear_viscosity]))) return
        if (self%bulk_modulus <= 0.0_wp .or. self%shear_modulus <= 0.0_wp) return
        if (self%bulk_viscosity < 0.0_wp .or. self%shear_viscosity < 0.0_wp) return
        if (.not. ieee_is_finite(2.0_wp*self%shear_modulus)) return
        if (.not. ieee_is_finite(self%bulk_modulus + (4.0_wp/3.0_wp)*self%shear_modulus)) return

        call spec%add_input("strain.small", nmandel, "mandel", "small_strain", "spatial", [ndim,ndim])
        call spec%add_output("stress.cauchy", nmandel, "mandel", "cauchy_stress", "spatial", [ndim,ndim])
        call spec%add_parameter("bulk_modulus")
        call spec%add_parameter("shear_modulus")
        call spec%add_parameter("bulk_viscosity")
        call spec%add_parameter("shear_viscosity")
        call spec%set_requires_dt(.true.)
        call spec%set_requires_initial_input(.true.)
        call spec%set_requires_previous_input(.true.)
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_kelvin_voigt

    subroutine evaluate_kelvin_voigt(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Calculate an uncommitted strain-rate stress at each material point.
        !!
        !! The previous strain `x_old` is the last accepted input. `state_old` and
        !! `state_trial` are empty. No per-point allocation or name lookup.
        class(kelvin_voigt_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted strain inputs `(ninput,npoints)` used to form the strain increment or rate.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `strain.small` inputs `(ninput,npoints)` at the current time step.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued history `(nstate,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Updated trial real-valued history `(nstate,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `stress.cauchy` outputs `(noutput,npoints)`.
        real(wp), intent(in) :: dt
            !! Prescribed time increment, strictly positive for time-dependent materials.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact discrete input Jacobian `(noutput,ninput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Per-point status codes. successful points must be assigned `cf_ok`.

        state_trial = state_old
        call evaluate_points(self, x_old, x_new, response, dt, point_status, tangent)
    end subroutine evaluate_kelvin_voigt

    subroutine sensitivity_kelvin_voigt(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Exact local derivatives of stress with respect to four coefficients.
        !!
        !! dsdp(6,nselected,npoints) follows setup-selected scalar parameters.
        !! The zero-sized sensitivity state is unchanged on every trial.
        class(kelvin_voigt_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted strain inputs `(ninput,npoints)` used to form the strain increment or rate.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `strain.small` inputs `(ninput,npoints)` at the current time step.
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
            !! Packed `stress.cauchy` outputs `(noutput,npoints)`.
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
        call evaluate_points(self, x_old, x_new, response, dt, point_status, tangent, selected_parameters, dsdp)
    end subroutine sensitivity_kelvin_voigt

    subroutine evaluate_points(self, x_old, x_new, response, dt, point_status, tangent, selected_parameters, dsdp)
        !! Evaluate elastic and viscous stresses for the discrete strain rate.
        class(kelvin_voigt_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted strain inputs `(ninput,npoints)` used to form the strain increment or rate.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `strain.small` inputs `(ninput,npoints)` at the current time step.
        real(wp), intent(out) :: response(:,:)
            !! Packed `stress.cauchy` outputs `(noutput,npoints)`.
        real(wp), intent(in) :: dt
            !! Prescribed time increment, strictly positive for time-dependent materials.
        integer, intent(inout) :: point_status(:)
            !! Per-point status codes. successful points must be assigned `cf_ok`.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact discrete input Jacobian `(noutput,ninput,npoints)`.
        integer, intent(in), optional :: selected_parameters(:)
            !! One-based IDs of scalar parameters in the requested derivative-axis order.
        real(wp), intent(out), optional :: dsdp(:,:,:)
            !! Exact response derivatives `(noutput,nselected,npoints)` at fixed prescribed input histories.
        real(wp) :: tangent_matrix(nmandel,nmandel)
        real(wp) :: dev_strain(nmandel)
        real(wp) :: dev_increment(nmandel)
        real(wp) :: increment(nmandel)
        real(wp) :: effective_bulk
        real(wp) :: effective_shear
        real(wp) :: viscous_bulk
        real(wp) :: viscous_shear
        real(wp) :: trace_strain
        real(wp) :: trace_increment
        real(wp) :: spherical
        integer :: ipoint
        integer :: i
        integer :: j

        point_status = cf_invalid_argument
        viscous_bulk = self%bulk_viscosity/dt
        viscous_shear = self%shear_viscosity/dt
        effective_bulk = self%bulk_modulus + viscous_bulk
        effective_shear = self%shear_modulus + viscous_shear
        if (.not. ieee_is_finite(effective_bulk)) return
        if (.not. ieee_is_finite(effective_shear)) return
        if (present(tangent)) then
            tangent_matrix = 0.0_wp
            tangent_matrix(1:ndim,1:ndim) = effective_bulk-(2.0_wp/3.0_wp)*effective_shear
            do i = 1, ndim
                tangent_matrix(i,i) = effective_bulk+(4.0_wp/3.0_wp)*effective_shear
            end do
            do i = ndim+1, nmandel
                tangent_matrix(i,i) = 2.0_wp*effective_shear
            end do
            if (.not. all(ieee_is_finite(tangent_matrix))) return
        end if

        point_loop: do ipoint = 1, size(x_new,2)
            if (.not. all(ieee_is_finite(x_new(:,ipoint)))) cycle point_loop
            if (.not. all(ieee_is_finite(x_old(:,ipoint)))) cycle point_loop
            trace_strain = sum(x_new(1:ndim,ipoint))
            increment = x_new(:,ipoint)-x_old(:,ipoint)
            trace_increment = sum(increment(1:ndim))
            dev_strain = x_new(:,ipoint)
            dev_strain(1:ndim) = dev_strain(1:ndim)-trace_strain/3.0_wp
            dev_increment = increment
            dev_increment(1:ndim) = dev_increment(1:ndim)-trace_increment/3.0_wp
            if (.not. ieee_is_finite(trace_strain)) cycle point_loop
            if (.not. ieee_is_finite(trace_increment)) cycle point_loop
            if (.not. all(ieee_is_finite(dev_strain))) cycle point_loop
            if (.not. all(ieee_is_finite(dev_increment))) cycle point_loop

            spherical = self%bulk_modulus*trace_strain + viscous_bulk*trace_increment
            response(:,ipoint) = 2.0_wp*self%shear_modulus*dev_strain + 2.0_wp*(viscous_shear*dev_increment)
            response(1:ndim,ipoint) = response(1:ndim,ipoint)+spherical
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop

            if (present(selected_parameters)) then
                parameter_loop: do j = 1, size(selected_parameters)
                    dsdp(:,j,ipoint) = 0.0_wp
                    select case (selected_parameters(j))
                    case (1)
                        dsdp(1:ndim,j,ipoint) = trace_strain
                    case (2)
                        dsdp(:,j,ipoint) = 2.0_wp*dev_strain
                    case (3)
                        dsdp(1:ndim,j,ipoint) = trace_increment/dt
                    case (4)
                        dsdp(:,j,ipoint) = 2.0_wp*(dev_increment/dt)
                    case default
                        exit parameter_loop
                    end select
                end do parameter_loop
                if (j <= size(selected_parameters)) cycle point_loop
                if (.not. all(ieee_is_finite(dsdp(:,:,ipoint)))) cycle point_loop
            end if

            if (present(tangent)) tangent(:,:,ipoint) = tangent_matrix
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_points

end module constituf_kelvin_voigt

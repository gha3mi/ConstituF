! SPDX-License-Identifier: Apache-2.0

module constituf_isotropic_thermoelastic
    !! Three-dimensional, small-strain, isotropic linear thermoelasticity.
    !!
    !! The prescribed fields are total infinitesimal strain
    !! \(\boldsymbol{\varepsilon}\) and absolute temperature \(T\), in a
    !! single shared spatial Cartesian frame. For temperature-independent
    !! moduli \(K>0\), \(G>0\), expansion coefficient \(\alpha\) and
    !! stress-free reference temperature \(T_0\), thermal strain is
    !! \[\boldsymbol{\varepsilon}_{\mathrm{th}}
    !! =\alpha(T-T_0)\boldsymbol{I},\]
    !! and the Cauchy stress is
    !! \[\boldsymbol{\sigma}
    !! =K[\mathrm{tr}(\boldsymbol{\varepsilon})-3\alpha(T-T_0)]
    !! \boldsymbol{I}+2G\mathrm{dev}(\boldsymbol{\varepsilon}).\]
    !! Here \(K,G\) share the chosen stress unit. The coefficient
    !! \(\alpha\) has inverse-temperature units and \(T,T_0\) use the
    !! same absolute-temperature scale.
    !! This is the Duhamel--Neumann constitutive relation with prescribed
    !! temperature, not a coupled heat-equation or entropy-balance solver.
    !!
    !! The strain and stress use the spatial symmetric 3D Mandel map
    !! \[\boldsymbol{A}_M=
    !! [A_{11},A_{22},A_{33},\sqrt{2}A_{23},\sqrt{2}A_{13},\sqrt{2}A_{12}]^T,\]
    !! where \(\boldsymbol{A}\) denotes either strain or Cauchy stress.
    !! Inputs are packed as `strain.small(1:6)`, then
    !! `temperature.absolute(7)`, both at the current time level.
    !! Output `stress.cauchy(1:6)` is positive in tension.
    !!
    !! The exact local \(6\times7\) Jacobian at fixed parameters is
    !! \[\frac{\partial\boldsymbol{\sigma}}{\partial
    !! \boldsymbol{\varepsilon}}=
    !! K\boldsymbol{I}\otimes\boldsymbol{I}
    !! +2G\mathbb{I}_{\mathrm{dev}},\qquad
    !! \frac{\partial\boldsymbol{\sigma}}{\partial T}
    !! =-3K\alpha\boldsymbol{I}.\]
    !! At fixed prescribed strain and temperature, the exact analytical
    !! parameter derivatives, in declaration order \(K,G,\alpha,T_0\), are
    !! \[\frac{\partial\boldsymbol{\sigma}}{\partial K}
    !! =[\mathrm{tr}(\varepsilon)-3\alpha(T-T_0)]\boldsymbol{I},\]
    !! \[\frac{\partial\boldsymbol{\sigma}}{\partial G}
    !! =2\mathrm{dev}(\varepsilon),\qquad
    !! \frac{\partial\boldsymbol{\sigma}}{\partial\alpha}
    !! =-3K(T-T_0)\boldsymbol{I},\qquad
    !! \frac{\partial\boldsymbol{\sigma}}{\partial T_0}
    !! =3K\alpha\boldsymbol{I}.\]
    !! These are local derivatives, not global equilibrium sensitivities.
    !!
    !! The moduli are temperature-independent, and the constitutive response
    !! is quasistatic, rate-independent and stateless. No heat flux, heat capacity,
    !! temperature-dependent moduli or finite thermal deformation is
    !! modeled. The solver owns temperature and thermal balance.
    !! Both temperatures must be finite, nonnegative and measured from
    !! absolute zero.
    !! Negative finite \(\alpha\) (thermal contraction upon heating) is
    !! allowed. Excessive inputs causing nonfinite stresses/derivatives
    !! fail the entire trial via per-point status.
    !!
    !! References: R. Abeyaratne, *Continuum Mechanics*, Volume II,
    !! section 11.1, linearized thermoelasticity, equations (11.19)--(11.23).
    !! https://web.mit.edu/abeyaratne/Volumes/RCA_Vol_II.pdf.
    !! MIT 16.20, *Structural Mechanics*, section 3.7
    !! "Thermoelastic Effects", equations (3.27)--(3.28).
    !! https://web.mit.edu/16.20/homepage/3_Constitutive/Constitutive_files/module_3_no_solutions.pdf.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: ndim = 3
    integer, parameter :: nmandel = 6
    integer, parameter :: ninput = nmandel + 1
    integer, parameter :: temperature_component = ninput

    public :: isotropic_thermoelastic_t

    type, extends(material_model_t) :: isotropic_thermoelastic_t
        !! Material parameters for isotropic thermoelastic Cauchy stress.
        !!
        !! Construct as `isotropic_thermoelastic_t(K,G,alpha,T0)`, then pass
        !! to `material_t%init`. The handle clones the model parameters.
        real(wp) :: bulk_modulus = 0.0_wp
            !! Bulk modulus \(K>0\), in the chosen stress unit.
        real(wp) :: shear_modulus = 0.0_wp
            !! Shear modulus \(G>0\), in the chosen stress unit.
        real(wp) :: thermal_expansion = 0.0_wp
            !! Linear expansion coefficient \(\alpha\), inverse temperature interval. May be negative.
        real(wp) :: reference_temperature = 0.0_wp
            !! Stress-free absolute temperature \(T_0\geq0\), in the chosen absolute-temperature unit.
    contains
        procedure :: describe => describe_thermoelastic
            !! Declare material fields, parameters and supported operations.
        procedure :: eval_batch => evaluate_thermoelastic
            !! Evaluate the batched material response, optional tangent and trial state.
        procedure :: eval_sensitivity_batch => sensitivity_thermoelastic
            !! Evaluate analytical parameter derivatives through the discrete material update.
    end type isotropic_thermoelastic_t

contains

    subroutine describe_thermoelastic(self, spec, status)
        !! Declare two independent input fields and one Mandel response.
        !!
        !! No history or `dt` is required. Invalid material parameters
        !! return `cf_invalid_argument` before the material is constructed.
        class(isotropic_thermoelastic_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema populated with fields, physical metadata and model capabilities.
        integer, intent(out) :: status
            !! Operation status code. `cf_ok` indicates success.

        status = cf_invalid_argument
        if (.not. ieee_is_finite(self%bulk_modulus)) return
        if (.not. ieee_is_finite(self%shear_modulus)) return
        if (.not. ieee_is_finite(self%thermal_expansion)) return
        if (.not. ieee_is_finite(self%reference_temperature)) return
        if (self%bulk_modulus <= 0.0_wp .or. self%shear_modulus <= 0.0_wp) return
        if (self%reference_temperature < 0.0_wp) return
        if (.not. ieee_is_finite(2.0_wp*self%shear_modulus)) return
        if (.not. ieee_is_finite(self%bulk_modulus + (4.0_wp/3.0_wp)*self%shear_modulus)) return
        if (.not. ieee_is_finite(3.0_wp*(self%bulk_modulus*self%thermal_expansion))) return

        call spec%add_input("strain.small", nmandel, "mandel", "small_strain", "spatial", [ndim,ndim])
        call spec%add_input("temperature.absolute", 1, "scalar", "absolute_temperature", "material")
        call spec%add_output("stress.cauchy", nmandel, "mandel", "cauchy_stress", "spatial", [ndim,ndim])
        call spec%add_parameter("bulk_modulus")
        call spec%add_parameter("shear_modulus")
        call spec%add_parameter("thermal_expansion")
        call spec%add_parameter("reference_temperature")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_thermoelastic

    subroutine evaluate_thermoelastic(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Compute stress at fixed current total strain and temperature.
        !!
        !! `x_new(ninput,npoints)`, `response(nmandel,npoints)` and optional
        !! `tangent(nmandel,ninput,npoints)` are caller-owned. No allocation
        !! or persistent state changes occur inside the batch loop.
        class(isotropic_thermoelastic_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previous packed inputs `(ninput,npoints)` when required by the model, otherwise a zero-sized array.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `strain.small` and `temperature.absolute` inputs `(ninput,npoints)` at the current time step.
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
        call evaluate_points(self, x_new, response, point_status, tangent)
    end subroutine evaluate_thermoelastic

    subroutine sensitivity_thermoelastic(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Fused stress, exact Jacobian and analytical parameter derivatives.
        !!
        !! Selected indices \(1,\ldots,4\) are resolved during initialization. `dsdp` has
        !! shape `(nmandel,nselected,npoints)` in the requested name order.
        !! Stateless derivative history is empty and requires no storage.
        class(isotropic_thermoelastic_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previous packed inputs `(ninput,npoints)` when required by the model, otherwise a zero-sized array.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `strain.small` and `temperature.absolute` inputs `(ninput,npoints)` at the current time step.
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
        call evaluate_points(self, x_new, response, point_status, tangent, selected_parameters, dsdp)
    end subroutine sensitivity_thermoelastic

    subroutine evaluate_points(self, x_new, response, point_status, tangent, selected_parameters, dsdp)
        !! One polymorphic dispatch per batch. Scalar parameter logic only
        !! executes when opt-in sensitivities have been requested.
        class(isotropic_thermoelastic_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `strain.small` and `temperature.absolute` inputs `(ninput,npoints)` at the current time step.
        real(wp), intent(out) :: response(:,:)
            !! Packed `stress.cauchy` outputs `(noutput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Per-point status codes. successful points must be assigned `cf_ok`.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact discrete input Jacobian `(noutput,ninput,npoints)`.
        integer, intent(in), optional :: selected_parameters(:)
            !! One-based IDs of scalar parameters in the requested derivative-axis order.
        real(wp), intent(out), optional :: dsdp(:,:,:)
            !! Exact response derivatives `(noutput,nselected,npoints)` at fixed prescribed input histories.
        real(wp) :: tangent_matrix(nmandel,ninput)
        real(wp) :: trace_strain
        real(wp) :: temperature_change
        real(wp) :: spherical
        real(wp) :: deviator(nmandel)
        integer :: ipoint
        integer :: i
        integer :: j

        if (present(tangent)) then
            tangent_matrix = 0.0_wp
            tangent_matrix(1:ndim,1:ndim) = self%bulk_modulus-(2.0_wp/3.0_wp)*self%shear_modulus
            do i = 1, ndim
                tangent_matrix(i,i) = self%bulk_modulus+(4.0_wp/3.0_wp)*self%shear_modulus
            end do
            do i = ndim+1, nmandel
                tangent_matrix(i,i) = 2.0_wp*self%shear_modulus
            end do
            tangent_matrix(1:ndim,temperature_component) = -3.0_wp*(self%bulk_modulus*self%thermal_expansion)
            if (.not. all(ieee_is_finite(tangent_matrix))) then
                point_status = cf_invalid_argument
                return
            end if
        end if

        point_loop: do ipoint = 1, size(x_new,2)
            point_status(ipoint) = cf_invalid_argument
            if (.not. all(ieee_is_finite(x_new(:,ipoint)))) cycle point_loop
            if (x_new(temperature_component,ipoint) < 0.0_wp) cycle point_loop

            trace_strain = sum(x_new(1:ndim,ipoint))
            temperature_change = x_new(temperature_component,ipoint) - self%reference_temperature
            deviator = x_new(1:nmandel,ipoint)
            deviator(1:ndim) = deviator(1:ndim)-trace_strain/3.0_wp
            spherical = self%bulk_modulus*(trace_strain- 3.0_wp*self%thermal_expansion*temperature_change)
            if (.not. ieee_is_finite(spherical)) cycle point_loop
            response(:,ipoint) = 2.0_wp*self%shear_modulus*deviator
            response(1:ndim,ipoint) = response(1:ndim,ipoint)+spherical
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop

            if (present(selected_parameters)) then
                parameter_loop: do j = 1, size(selected_parameters)
                    dsdp(:,j,ipoint) = 0.0_wp
                    select case (selected_parameters(j))
                    case (1)
                        dsdp(1:ndim,j,ipoint) = trace_strain-3.0_wp*self%thermal_expansion*temperature_change
                    case (2)
                        dsdp(:,j,ipoint) = 2.0_wp*deviator
                    case (3)
                        dsdp(1:ndim,j,ipoint) = -3.0_wp*(self%bulk_modulus*temperature_change)
                    case (4)
                        dsdp(1:ndim,j,ipoint) = 3.0_wp*(self%bulk_modulus*self%thermal_expansion)
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

end module constituf_isotropic_thermoelastic

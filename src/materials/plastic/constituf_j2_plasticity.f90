! SPDX-License-Identifier: Apache-2.0

module constituf_j2_plasticity
    !! Three-dimensional, infinitesimal-strain, associative \(J_2\) plasticity.
    !!
    !! This model uses an additive elastic/plastic strain split and linear
    !! isotropic hardening. At any prescribed Mandel strain
    !! \(\boldsymbol{\varepsilon}\), with committed plastic strain
    !! \(\boldsymbol{\varepsilon}^{p}_n\) and accumulated equivalent plastic
    !! strain \(\alpha_n\), the elastic predictor is
    !! \[\boldsymbol{s}_{\mathrm{tr}}=
    !! 2G(\mathrm{dev}\boldsymbol{\varepsilon}
    !! -\boldsymbol{\varepsilon}^{p}_n),\qquad
    !! q_{\mathrm{tr}}=\sqrt{3/2}\|\boldsymbol{s}_{\mathrm{tr}}\|.\]
    !! Define
    !! \[f_{\mathrm{tr}}=q_{\mathrm{tr}}-(\sigma_y+H\alpha_n).\]
    !! For \(f_{\mathrm{tr}}\leq0\), the update is elastic. Otherwise the
    !! backward-Euler radial-return solution with \(G>0\), \(H\geq0\) is
    !! \[\Delta\gamma=\frac{f_{\mathrm{tr}}}{3G+H},\qquad
    !! \boldsymbol{n}=\frac{3\boldsymbol{s}_{\mathrm{tr}}}{2q_{\mathrm{tr}}},\]
    !! \[\boldsymbol{\varepsilon}^{p}_{n+1}
    !! =\boldsymbol{\varepsilon}^{p}_n+\Delta\gamma\boldsymbol{n},\qquad
    !! \alpha_{n+1}=\alpha_n+\Delta\gamma,\]
    !! \[\boldsymbol{\sigma}_{n+1}=
    !! K\,\mathrm{tr}(\boldsymbol{\varepsilon})\boldsymbol{I}
    !! +\boldsymbol{s}_{\mathrm{tr}}-2G\Delta\gamma\boldsymbol{n}.\]
    !! Here \(\boldsymbol{n}:\boldsymbol{n}=3/2\), so \(\Delta\gamma\)
    !! is the equivalent plastic-strain increment. There is no local Newton
    !! iteration for linear hardening.
    !! The exact *algorithmic*, not continuum elastic, Mandel tangent on
    !! the plastic branch is
    !! \[\mathbb{C}_{\mathrm{alg}}=K\boldsymbol{I}\otimes\boldsymbol{I}
    !! +2G\left(1-\frac{3G\Delta\gamma}{q_{\mathrm{tr}}}\right)
    !! \mathbb{I}_{\mathrm{dev}}
    !! +4G^2\left(\frac{\Delta\gamma}{q_{\mathrm{tr}}}
    !! -\frac{1}{3G+H}\right)\boldsymbol{n}\otimes\boldsymbol{n}.\]
    !! On the elastic branch, \(\mathbb{C}_{\mathrm{alg}}=
    !! K\boldsymbol{I}\otimes\boldsymbol{I}+2G\mathbb{I}_{\mathrm{dev}}\).
    !!
    !! Positive stress/strain denotes tension. Input `strain.small` is
    !! dimensionless. Response `stress.cauchy` is in the chosen stress unit. Both are symmetric
    !! \(3\times3\) spatial tensors. Each such tensor \(\boldsymbol{A}\)
    !! is packed as
    !! \[\operatorname{mandel}(\boldsymbol{A})=
    !! [A_{11},A_{22},A_{33},\sqrt{2}A_{23},
    !! \sqrt{2}A_{13},\sqrt{2}A_{12}]^{\mathsf T}.\]
    !! History contains the six Mandel plastic-strain components and scalar
    !! \(\alpha\geq0\). Initial state and all initial path derivatives are
    !! zero (virgin, initially unplastified solid). No `dt` is required.
    !! Parameters: positive bulk modulus \(K\), positive shear modulus \(G\),
    !! positive yield stress \(\sigma_y\) and nonnegative hardening
    !! modulus \(H\), all in the same stress unit. \(H=0\) is perfect plasticity.
    !!
    !! Optional analytical parameter sensitivities differentiate the exact
    !! *discrete* return algorithm through all accepted history steps.
    !! For any selected \(p\in\{K,G,\sigma_y,H\}\), with prescribed strains
    !! held fixed and \(\dot{(\ )}=\mathrm{d}(\ )/\mathrm{d}p\):
    !! \[\dot{\boldsymbol{s}}_{\mathrm{tr}}
    !! =2\dot G(\mathrm{dev}\varepsilon-\varepsilon_n^p)
    !! -2G\dot{\varepsilon}_n^p,\quad
    !! \dot q_{\mathrm{tr}}=\boldsymbol{n}:
    !! \dot{\boldsymbol{s}}_{\mathrm{tr}},\]
    !! \[\dot{\Delta\gamma}=
    !! \frac{\dot q_{\mathrm{tr}}-\dot{\sigma}_y-\dot H\alpha_n
    !! -H\dot\alpha_n-\Delta\gamma(3\dot G+\dot H)}{3G+H},\]
    !! \[\dot{\boldsymbol n}=\frac{3}{2q_{\mathrm{tr}}}
    !! \left(\dot{\boldsymbol s}_{\mathrm{tr}}-
    !! \frac{\boldsymbol s_{\mathrm{tr}}}{q_{\mathrm{tr}}}
    !! \dot q_{\mathrm{tr}}\right).\]
    !! Updated state/stress derivatives follow by the product rule.
    !! These derivatives are local: they do NOT include global equilibrium
    !! or parameter-dependent prescribed-strain derivatives.
    !! Exactly at yield activation the classical derivative is not unique.
    !! Sensitivity-enabled evaluations return `cf_nondifferentiable` for
    !! \(|f_{\mathrm{tr}}|\) within a scale-aware roundoff band. Ordinary
    !! evaluation uses the elastic side when \(f_{\mathrm{tr}}=0\).
    !!
    !! Restrictions: isothermal, small strain/rotation, rate-independent,
    !! associative von Mises yield with isotropic linear hardening only.
    !! No plane-stress reduction, kinematic hardening, finite strain or
    !! damage is included. Invalid physical inputs invalidate the trial.
    !!
    !! References: J. C. Simo and R. L. Taylor, *Consistent Tangent
    !! Operators for Rate-Independent Elastoplasticity*, Computer Methods
    !! in Applied Mechanics and Engineering 48 (1985), 101--118,
    !! DOI: 10.1016/0045-7825(85)90070-2. Sections 3 and 4,
    !! Eqs. (3.8) and (4.12)--(4.13), specialized to zero back stress and
    !! linear isotropic hardening. The paper uses a unit flow normal.
    !! Its plastic multiplier increment is \(\sqrt{3/2}\Delta\gamma\)
    !! in the normalization used here.
    !! J. C. Simo and T. J. R. Hughes, *Computational Inelasticity*,
    !! Springer, 1998, DOI: 10.1007/b98904, ISBN 978-0-387-97520-7.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument, cf_nondifferentiable
    implicit none (type, external)
    private

    integer, parameter :: nmandel = 6
    integer, parameter :: nstate = nmandel + 1

    public :: j2_plasticity_t

    type, extends(material_model_t) :: j2_plasticity_t
        !! Virgin three-dimensional \(J_2\) solid with linear isotropic hardening.
        !!
        !! Construct as `j2_plasticity_t(K,G,yield_stress,hardening_modulus)`
        !! and call `material_t%init`. All history is owned by `material_t`.
        !! Every `eval` starts from the last committed state.
        real(wp) :: bulk_modulus = 0.0_wp
            !! \(K>0\), bulk modulus in the chosen stress unit.
        real(wp) :: shear_modulus = 0.0_wp
            !! \(G>0\), shear modulus in the chosen stress unit.
        real(wp) :: yield_stress = 0.0_wp
            !! \(\sigma_y>0\), initial yield stress in the chosen stress unit.
        real(wp) :: hardening_modulus = 0.0_wp
            !! \(H\geq0\), isotropic linear hardening in the chosen stress unit.
    contains
        procedure :: describe => describe_j2
            !! Describe the material input, response, state and scalar parameter schema.
        procedure :: eval_batch => evaluate_j2
            !! Evaluate the batched material response and trial history.
        procedure :: initialize_sensitivity => initialize_j2_sensitivity
            !! Initialize analytical derivatives of the reference real history.
        procedure :: eval_sensitivity_batch => sensitivity_j2
            !! Evaluate exact selected-parameter sensitivities through the discrete material update.
    end type j2_plasticity_t

contains

    subroutine describe_j2(self, spec, status)
        !! Specify spatial Mandel fields and explicit real-valued history.
        class(j2_plasticity_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema receiving the packed field descriptions and supported capabilities.
        integer, intent(out) :: status
            !! Status code; `cf_ok` indicates success.
        real(wp) :: denominator

        status = cf_invalid_argument
        if (.not. ieee_is_finite(self%bulk_modulus)) return
        if (.not. ieee_is_finite(self%shear_modulus)) return
        if (.not. ieee_is_finite(self%yield_stress)) return
        if (.not. ieee_is_finite(self%hardening_modulus)) return
        if (self%bulk_modulus <= 0.0_wp .or. self%shear_modulus <= 0.0_wp) return
        if (self%yield_stress <= 0.0_wp .or. self%hardening_modulus < 0.0_wp) return
        denominator = 3.0_wp*self%shear_modulus + self%hardening_modulus
        if (.not. ieee_is_finite(denominator)) return
        if (.not. ieee_is_finite(self%bulk_modulus + 4.0_wp*self%shear_modulus/3.0_wp)) return

        call spec%add_input("strain.small", nmandel, "mandel", "small_strain", "spatial", [3,3])
        call spec%add_output("stress.cauchy", nmandel, "mandel", "cauchy_stress", "spatial", [3,3])
        call spec%add_state("plastic.strain", nmandel, "mandel", "plastic_strain", "spatial", [3,3])
        call spec%add_state("plastic.accumulated", 1, "scalar", "equivalent_plastic_strain", "material")
        call spec%add_parameter("bulk_modulus")
        call spec%add_parameter("shear_modulus")
        call spec%add_parameter("yield_stress")
        call spec%add_parameter("hardening_modulus")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_j2

    subroutine initialize_j2_sensitivity(self, initial_input, initial_state, selected_parameters, sensitivity, status)
        !! Virgin zero plastic strain and \(\alpha\) are parameter-independent.
        class(j2_plasticity_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        real(wp), intent(in) :: initial_input(:,:)
            !! Prescribed input state used for material initialization.
        real(wp), intent(in) :: initial_state(:,:)
            !! Initial real-valued history `(nstate,npoints)`.
        integer, intent(in) :: selected_parameters(:)
            !! One-based parameter IDs in the caller's selected derivative-axis order.
        real(wp), intent(out) :: sensitivity(:,:,:)
            !! Initial exact history derivatives `(nstate,nselected,npoints)`.
        integer, intent(out) :: status
            !! Status code; `cf_ok` indicates success.

        sensitivity = 0.0_wp
        status = cf_ok
    end subroutine initialize_j2_sensitivity

    subroutine evaluate_j2(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Discrete radial-return update. No allocations in the batch kernel.
        !!
        !! Optional tangent(6,6,npoints) is the exact derivative of the
        !! updated Cauchy stress with respect to current Mandel strain
        !! at fixed committed history and parameters.
        class(j2_plasticity_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed inputs `(ninput,npoints)` when required, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `strain.small` inputs `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued history `(nstate,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Computed uncommitted real-valued history `(nstate,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `stress.cauchy` responses `(noutput,npoints)`.
        real(wp), intent(in) :: dt
            !! Prescribed time increment; a positive value is required for time-dependent laws.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic tangent `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(inout) :: point_status(:)
            !! Per-point outcome; write `cf_ok` only for a successful constitutive update.

        call evaluate_points(self, x_new, state_old, state_trial, response, point_status, tangent)
    end subroutine evaluate_j2

    subroutine sensitivity_j2(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Exact local parameter-path derivatives of each return-map branch.
        !!
        !! dsdp(6,nselected,npoints) preserves the requested parameter
        !! order. sensitivity_trial(7,nselected,npoints) is retained and
        !! rolled back/committed with the plastic history by `material_t`.
        class(j2_plasticity_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed inputs `(ninput,npoints)` when required, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `strain.small` inputs `(ninput,npoints)`.
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
            !! Packed `stress.cauchy` responses `(noutput,npoints)`.
        real(wp), intent(out) :: dsdp(:,:,:)
            !! Exact total response derivatives `(noutput,nselected,npoints)` at fixed prescribed driving histories.
        real(wp), intent(in) :: dt
            !! Prescribed time increment; a positive value is required for time-dependent laws.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic tangent `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(inout) :: point_status(:)
            !! Per-point outcome; write `cf_ok` only for a successful constitutive update.

        call evaluate_points(self, x_new, state_old, state_trial, response, point_status, tangent, selected_parameters, sensitivity_old, sensitivity_trial, dsdp)
    end subroutine sensitivity_j2

    subroutine evaluate_points(self, x_new, state_old, state_trial, response, point_status, tangent, selected_parameters, sensitivity_old, sensitivity_trial, dsdp)
        !! Shared path computes the same stress and tangent with/without dsdp.
        class(j2_plasticity_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `strain.small` inputs `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued history `(nstate,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Computed uncommitted real-valued history `(nstate,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `stress.cauchy` responses `(noutput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Per-point outcome; write `cf_ok` only for a successful constitutive update.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic tangent `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(in), optional :: selected_parameters(:)
            !! One-based parameter IDs in the caller's selected derivative-axis order.
        real(wp), intent(in), optional :: sensitivity_old(:,:,:)
            !! Committed real-history parameter derivatives `(nstate,nselected,npoints)`.
        real(wp), intent(out), optional :: sensitivity_trial(:,:,:)
            !! Updated real-history parameter derivatives `(nstate,nselected,npoints)`.
        real(wp), intent(out), optional :: dsdp(:,:,:)
            !! Exact total response derivatives `(noutput,nselected,npoints)` at fixed prescribed driving histories.
        real(wp) :: strain_dev(nmandel)
        real(wp) :: trial_stress(nmandel)
        real(wp) :: direction(nmandel)
        real(wp) :: elastic_dev(nmandel)
        real(wp) :: plastic_stress(nmandel)
        real(wp) :: derivative_trial(nmandel)
        real(wp) :: derivative_direction(nmandel)
        real(wp) :: projector(nmandel,nmandel)
        real(wp) :: tangent_matrix(nmandel,nmandel)
        real(wp) :: trace_strain
        real(wp) :: equivalent_trial
        real(wp) :: yield_function
        real(wp) :: denominator
        real(wp) :: increment
        real(wp) :: effective_yield
        real(wp) :: tangent_factor
        real(wp) :: d_bulk
        real(wp) :: d_shear
        real(wp) :: d_yield
        real(wp) :: d_hardening
        real(wp) :: derivative_equivalent
        real(wp) :: derivative_increment
        real(wp) :: derivative_denominator
        integer :: ipoint
        integer :: iparam
        integer :: i
        integer :: j
        logical :: plastic_branch

        denominator = 3.0_wp*self%shear_modulus + self%hardening_modulus
        if (present(tangent)) call deviator_projector(projector)

        point_loop: do ipoint = 1, size(x_new,2)
            point_status(ipoint) = cf_invalid_argument
            if (.not. all(ieee_is_finite(x_new(:,ipoint)))) cycle point_loop
            if (.not. all(ieee_is_finite(state_old(:,ipoint)))) cycle point_loop
            if (present(selected_parameters)) then
                if (.not. all(ieee_is_finite(sensitivity_old(:,:,ipoint)))) then
                    cycle point_loop
                end if
            end if

            trace_strain = sum(x_new(1:3,ipoint))
            strain_dev = x_new(:,ipoint)
            strain_dev(1:3) = strain_dev(1:3) - trace_strain/3.0_wp
            elastic_dev = strain_dev - state_old(1:nmandel,ipoint)
            trial_stress = 2.0_wp*self%shear_modulus*elastic_dev
            if (.not. all(ieee_is_finite(trial_stress))) cycle point_loop
            equivalent_trial = sqrt(1.5_wp*dot_product(trial_stress,trial_stress))
            if (.not. ieee_is_finite(equivalent_trial)) cycle point_loop
            effective_yield = self%yield_stress + self%hardening_modulus*state_old(nstate,ipoint)
            if (.not. ieee_is_finite(effective_yield)) cycle point_loop
            yield_function = equivalent_trial - effective_yield
            if (.not. ieee_is_finite(yield_function)) cycle point_loop

            plastic_branch = yield_function > 0.0_wp
            increment = 0.0_wp
            direction = 0.0_wp
            plastic_stress = trial_stress
            if (plastic_branch) then
                increment = yield_function/denominator
                direction = (1.5_wp/equivalent_trial)*trial_stress
                plastic_stress = trial_stress - 2.0_wp*self%shear_modulus* increment*direction
                if (.not. ieee_is_finite(increment)) cycle point_loop
                if (.not. all(ieee_is_finite(direction))) cycle point_loop
                if (.not. all(ieee_is_finite(plastic_stress))) cycle point_loop
            end if

            response(:,ipoint) = plastic_stress
            response(1:3,ipoint) = response(1:3,ipoint) + self%bulk_modulus*trace_strain
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop

            if (present(tangent)) then
                tangent_factor = 1.0_wp
                if (plastic_branch) then
                    tangent_factor = 1.0_wp - 3.0_wp*self%shear_modulus*increment/equivalent_trial
                end if
                tangent_matrix = 2.0_wp*self%shear_modulus*tangent_factor*projector
                tangent_matrix(1:3,1:3) = tangent_matrix(1:3,1:3) + self%bulk_modulus
                if (plastic_branch) then
                    tangent_factor = 4.0_wp*self%shear_modulus**2* (increment/equivalent_trial - 1.0_wp/denominator)
                    do j = 1, nmandel
                        do i = 1, nmandel
                            tangent_matrix(i,j) = tangent_matrix(i,j) + tangent_factor*direction(i)*direction(j)
                        end do
                    end do
                end if
                if (.not. all(ieee_is_finite(tangent_matrix))) cycle point_loop
            end if

            if (present(selected_parameters)) then
                if (abs(yield_function)/max(self%yield_stress,equivalent_trial,abs(effective_yield)) <= 32.0_wp*epsilon(1.0_wp)) then
                    point_status(ipoint) = cf_nondifferentiable
                    cycle point_loop
                end if
                parameter_loop: do iparam = 1, size(selected_parameters)
                    d_bulk = 0.0_wp
                    d_shear = 0.0_wp
                    d_yield = 0.0_wp
                    d_hardening = 0.0_wp
                    select case (selected_parameters(iparam))
                    case (1)
                        d_bulk = 1.0_wp
                    case (2)
                        d_shear = 1.0_wp
                    case (3)
                        d_yield = 1.0_wp
                    case (4)
                        d_hardening = 1.0_wp
                    case default
                        exit parameter_loop
                    end select
                    derivative_trial = 2.0_wp*d_shear*elastic_dev - 2.0_wp*self%shear_modulus* sensitivity_old(1:nmandel,iparam,ipoint)
                    sensitivity_trial(:,iparam,ipoint) = sensitivity_old(:,iparam,ipoint)
                    dsdp(:,iparam,ipoint) = derivative_trial
                    if (plastic_branch) then
                        derivative_equivalent = dot_product(direction,derivative_trial)
                        derivative_denominator = 3.0_wp*d_shear+d_hardening
                        derivative_increment = (derivative_equivalent-d_yield- d_hardening*state_old(nstate,ipoint) - self%hardening_modulus*sensitivity_old(nstate,iparam,ipoint) - increment*derivative_denominator)/denominator
                        derivative_direction = (1.5_wp/equivalent_trial)* (derivative_trial- (trial_stress/equivalent_trial)*derivative_equivalent)
                        sensitivity_trial(1:nmandel,iparam,ipoint) = sensitivity_old(1:nmandel,iparam,ipoint) + derivative_increment*direction + increment*derivative_direction
                        sensitivity_trial(nstate,iparam,ipoint) = sensitivity_old(nstate,iparam,ipoint)+derivative_increment
                        dsdp(:,iparam,ipoint) = derivative_trial - 2.0_wp*d_shear*increment*direction - 2.0_wp*self%shear_modulus*(derivative_increment*direction + increment*derivative_direction)
                    end if
                    dsdp(1:3,iparam,ipoint) = dsdp(1:3,iparam,ipoint) + d_bulk*trace_strain
                end do parameter_loop
                if (iparam <= size(selected_parameters)) cycle point_loop
                if (.not. all(ieee_is_finite(dsdp(:,:,ipoint)))) cycle point_loop
                if (.not. all(ieee_is_finite(sensitivity_trial(:,:,ipoint)))) then
                    cycle point_loop
                end if
            end if

            state_trial(:,ipoint) = state_old(:,ipoint)
            if (plastic_branch) then
                state_trial(1:nmandel,ipoint) = state_old(1:nmandel,ipoint) + increment*direction
                state_trial(nstate,ipoint) = state_old(nstate,ipoint) + increment
            end if
            if (.not. all(ieee_is_finite(state_trial(:,ipoint)))) cycle point_loop
            if (present(tangent)) tangent(:,:,ipoint) = tangent_matrix
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_points

    pure subroutine deviator_projector(projector)
        !! Orthogonal symmetric-deviatoric projector in 3D Mandel notation.
        real(wp), intent(out) :: projector(nmandel,nmandel)
            !! Six-by-six Mandel projector onto deviatoric symmetric tensors.
        integer :: i

        projector = 0.0_wp
        projector(1:3,1:3) = -1.0_wp/3.0_wp
        do i = 1, nmandel
            projector(i,i) = projector(i,i) + 1.0_wp
        end do
    end subroutine deviator_projector

end module constituf_j2_plasticity

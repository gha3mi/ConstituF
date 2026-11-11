! SPDX-License-Identifier: Apache-2.0

module constituf_ductile_damage
    !! Small-strain \(J_2\) plasticity with irreversible, scalar ductile damage.
    !!
    !! This is an independently implemented **Lemaitre-inspired**
    !! strain-threshold damage specialization, not a reproduction of the
    !! fully coupled, stress-triaxiality-dependent Lemaitre law. The
    !! undamaged effective response follows the existing
    !! [[constituf_j2_plasticity(module)]] radial-return algorithm with
    !! bulk modulus \(K>0\), shear modulus \(G>0\), yield stress \(Y_0>0\) and
    !! linear isotropic hardening \(H\geq0\), all in the chosen stress unit. Strain, plastic strain
    !! and Cauchy stress are symmetric three-dimensional spatial tensors.
    !! Each such tensor \(\boldsymbol{A}\) is packed as
    !! \[\operatorname{mandel}(\boldsymbol{A})=
    !! [A_{11},A_{22},A_{33},\sqrt{2}A_{23},
    !! \sqrt{2}A_{13},\sqrt{2}A_{12}]^{\mathsf T}.\]
    !! Tension is positive.
    !!
    !! After the effective-stress backward-Euler update, with plastic
    !! strain \(\boldsymbol{\varepsilon}^p\) and accumulated equivalent plastic strain \(\alpha\), define
    !! the nonnegative damage variable
    !! \[D(\alpha)=D_{\max}
    !! \bigl[1-\exp(-b\langle\alpha-\alpha_0\rangle)\bigr],
    !! \quad 0\leq D_{\max}<1,\quad b>0,\quad\alpha_0\geq0.\]
    !! The positive part is zero before damage initiation. There is no
    !! separate damage return equation or stress-triaxiality multiplier.
    !! Its derivative on the active branch is
    !! \[D_\alpha=bD_{\max}\exp[-b(\alpha-\alpha_0)].\]
    !! Let \(\boldsymbol{\varepsilon}_e=\boldsymbol{\varepsilon}-\boldsymbol{\varepsilon}^p\) and
    !! \[W_{\rm e}=\frac K2[\mathrm{tr}(\varepsilon_e)]^2+
    !! G\,\mathrm{dev}(\varepsilon_e):
    !! \mathrm{dev}(\varepsilon_e),\quad
    !! Y=W_{\rm e}+\frac H2\alpha^2,\]
    !! Here \(Y=-\partial\Psi/\partial D\) has units of stress and is the
    !! positive elastic-plus-hardening energy *release rate* conjugate
    !! to \(D\) (not the yield stress \(Y_0\)).
    !! Degrading both elastic and hardening energy is a choice of this
    !! adaptation. Lemaitre's cited coupled law degrades the elastic term.
    !! The Helmholtz free-energy density and nominal degraded Cauchy
    !! stress are
    !! \[\Psi=(1-D)Y,\qquad
    !! \boldsymbol\sigma=(1-D)\boldsymbol\sigma_{\rm eff},\quad
    !! \boldsymbol\sigma_{\rm eff}=
    !! \mathbb C_{\rm e}:(\boldsymbol\varepsilon-\boldsymbol\varepsilon^p).\]
    !! Effective-stress plasticity remains \(J_2\) with effective yield
    !! \(q_{\rm eff}=Y_0+H\alpha\) on plastic loading.
    !!
    !! Continuous isothermal dissipation along an admissible plastic
    !! trajectory is
    !! \[\mathcal D=(1-D)Y_0\dot\alpha+Y\dot D\ge0.\]
    !! The *reported discrete increment* is the explicit end-point
    !! quadrature diagnostic
    !! \[\Delta\mathcal D_{\rm diag}=
    !! (1-D_{n+1})Y_0(\alpha_{n+1}-\alpha_n)
    !! +Y_{n+1}(D_{n+1}-D_n)\ge0.\]
    !! This nonnegative local quantity is a quadrature approximation of
    !! the continuous dissipation. It must NOT be interpreted as an
    !! exact work-energy residual for a finite increment. The solver
    !! is responsible for any mesh-objective fracture regularization.
    !!
    !! The inelastic model returns three distinct outputs:
    !! `stress.cauchy` (6 Mandel components, stress units),
    !! `energy.helmholtz` (1 scalar, energy per volume) and
    !! `energy.dissipation_increment` (1 scalar, energy per volume).
    !! History contains plastic.strain (6), plastic.accumulated (1) and
    !! damage.scalar (1), all dimensionless. It is committed/reverted
    !! atomically by the existing `material_t` handle.
    !!
    !! The \(8\times6\) exact **algorithmic tangent** at fixed committed history
    !! is obtained by differentiating the discrete \(J_2\) radial return and
    !! the damage degradation simultaneously. On a smooth branch:
    !! \[\partial_j\sigma_i=(1-D)C^{\rm J2}_{ij}
    !! -\sigma_{{\rm eff},i}D_\alpha\,\partial_j\alpha,\]
    !! \[\partial_j\Psi=(1-D)\partial_jY-YD_\alpha\partial_j\alpha,\]
    !! \[\partial_j\Delta\mathcal D_{\rm diag}
    !! =(1-D)Y_0\partial_j\alpha
    !! -Y_0\Delta\alpha D_\alpha\partial_j\alpha
    !! +\Delta D\partial_jY+YD_\alpha\partial_j\alpha.\]
    !! The local derivatives \(\partial\alpha/\partial\varepsilon\)
    !! and \(\partial\varepsilon^p/\partial\varepsilon\) follow by
    !! differentiation of the *same* \(J_2\) return. No finite differences
    !! are computed in material evaluation. At initial yielding or
    !! \(\alpha=\alpha_0\) a classical branch derivative need not exist. The
    !! inactive-side convention is used for ordinary tangent evaluations
    !! at exact equality. Sensitivity-enabled trials instead report
    !! `cf_nondifferentiable` in roundoff bands around these switches.
    !!
    !! No rate, temperature, length-scale, crack-normal, gradient,
    !! nonlocal interaction, finite-strain kinematics, or fully coupled
    !! Lemaitre triaxial damage force is implemented. Analytical parameter
    !! sensitivities propagate \(J_2\) plastic history and explicitly
    !! differentiate the threshold damage, Helmholtz energy and diagnostic
    !! dissipation increment. No parameter perturbation is used.
    !! Local softening is generally mesh-dependent in FEM. Do not claim
    !! fracture-energy convergence without a host-solver regularization.
    !!
    !! Scientific background: J. Lemaitre, "A Continuous Damage Mechanics
    !! Model for Ductile Fracture", J. Eng. Mater. Technol. 107 (1985),
    !! 83--89, DOI 10.1115/1.3225775. J. Lemaitre, "Coupled
    !! elasto-plasticity and damage constitutive equations", CMAME 51
    !! (1985), 31--49, DOI 10.1016/0045-7825(85)90026-X.
    !! Section 2.2, Eqs. (2.3) and (2.5)--(2.8), gives the effective
    !! stress and elastic-only damage-energy coupling, with a nonpositive
    !! damage force \(Y\). The present positive release rate also includes
    !! hardening energy.
    !! The \(\alpha\)-threshold exponential law and its explicit effective
    !! return-map coupling here are original adaptations, not claimed
    !! to be the original Lemaitre triaxiality-dependent damage law.
    !! Code provenance: independently authored from the stated equations
    !! and the project's own \(J_2\) solver. No third-party code was copied.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument, cf_evaluation_failed, cf_nondifferentiable
    use constituf_j2_plasticity, only: j2_plasticity_t
    implicit none (type, external)
    private

    integer, parameter :: nmandel = 6
    integer, parameter :: nbase_state = 7
    integer, parameter :: nstate = 8
    integer, parameter :: noutput = 8

    public :: ductile_damage_t

    type, extends(material_model_t) :: ductile_damage_t
        !! Construct with `ductile_damage_t(K,G,Y0,H,damage_threshold,damage_rate,maximum_damage)`.
        !! No external state is supplied. `npoints` is set in `material_t%init`.
        real(wp) :: bulk_modulus = 0.0_wp
            !! \(K>0\) in the chosen stress unit.
        real(wp) :: shear_modulus = 0.0_wp
            !! \(G>0\) in the chosen stress unit.
        real(wp) :: yield_stress = 0.0_wp
            !! Effective \(J_2\) yield stress \(Y_0>0\) in the chosen stress unit.
        real(wp) :: hardening_modulus = 0.0_wp
            !! Effective isotropic hardening \(H\geq0\) in the chosen stress unit.
        real(wp) :: damage_threshold = 0.0_wp
            !! Accumulated equivalent strain threshold \(\alpha_0\geq0\).
        real(wp) :: damage_rate = 1.0_wp
            !! Dimensionless positive post-threshold exponential rate \(b\).
        real(wp) :: maximum_damage = 0.0_wp
            !! Saturation \(D_{\max}\in[0,1)\). Zero exactly reduces to \(J_2\).
    contains
        procedure :: describe => describe_ductile_damage
            !! Declare the constitutive fields, parameter metadata and supported capabilities.
        procedure :: eval_batch => evaluate_ductile_damage
            !! Evaluate the batched constitutive response, optional tangent and trial history.
        procedure :: initialize_sensitivity => initialize_damage_sensitivity
            !! Initialize analytical sensitivities of the reference material history.
        procedure :: eval_sensitivity_batch => sensitivity_damage
            !! Evaluate exact parameter-path derivatives and corresponding trial history.
    end type ductile_damage_t

contains

    subroutine describe_ductile_damage(self,spec,status)
        !! Validate the parameters and declare `strain.small` input and `stress.cauchy` and `energy.helmholtz` and `energy.dissipation_increment` output with their physical metadata.
        class(ductile_damage_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema receiving field layouts, parameters and supported operations.
        integer, intent(out) :: status
            !! Operation status. `cf_ok` indicates a valid model or successful initialization.
        type(j2_plasticity_t) :: base
        type(material_spec_t) :: base_spec

        status = cf_invalid_argument
        if (.not. all(ieee_is_finite([self%bulk_modulus,self%shear_modulus,self%yield_stress, self%hardening_modulus,self%damage_threshold,self%damage_rate,self%maximum_damage]))) return
        if (self%damage_threshold < 0.0_wp .or. self%damage_rate <= 0.0_wp) return
        if (self%maximum_damage < 0.0_wp .or. self%maximum_damage >= 1.0_wp) return
        base = j2_plasticity_t(self%bulk_modulus,self%shear_modulus,self%yield_stress,self%hardening_modulus)
        call base%describe(base_spec,status)
        if (status /= cf_ok) return

        call spec%add_input("strain.small", nmandel, "mandel", "small_strain", "spatial", [3,3])
        call spec%add_output("stress.cauchy", nmandel, "mandel", "cauchy_stress", "spatial", [3,3])
        call spec%add_output("energy.helmholtz", 1, "scalar", "helmholtz_energy", "material")
        call spec%add_output("energy.dissipation_increment", 1, "scalar", "dissipation_increment", "material")
        call spec%add_state("plastic.strain", nmandel, "mandel", "plastic_strain", "spatial", [3,3])
        call spec%add_state("plastic.accumulated", 1, "scalar", "equivalent_plastic_strain", "material")
        call spec%add_state("damage.scalar", 1, "scalar", "isotropic_damage", "material")
        call spec%add_parameter("bulk_modulus")
        call spec%add_parameter("shear_modulus")
        call spec%add_parameter("yield_stress")
        call spec%add_parameter("hardening_modulus")
        call spec%add_parameter("damage_threshold")
        call spec%add_parameter("damage_rate")
        call spec%add_parameter("maximum_damage")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_ductile_damage

    subroutine evaluate_ductile_damage(self,x_old,x_new,state_old,state_trial,response,dt,tangent,point_status)
        !! \(J_2\) effective return followed by damage degradation and energies.
        !! All per-point scratch arrays have fixed bounds. Production
        !! evaluation has no heap allocation or I/O. The underlying
        !! effective \(J_2\) call uses exactly its original return algorithm.
        class(ductile_damage_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed input `(ninput,npoints)` if required, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `strain.small` input `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued state `(nstate,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Trial real-valued state `(nstate,npoints)` returned by this evaluation.
        real(wp), intent(out) :: response(:,:)
            !! Packed `stress.cauchy` and `energy.helmholtz` and `energy.dissipation_increment` response `(noutput,npoints)` for the evaluated material points.
        real(wp), intent(in) :: dt
            !! Prescribed time increment. It is positive when the material requires time integration.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact discrete response Jacobian `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(inout) :: point_status(:)
            !! Per-point status values; successful evaluations must set each entry to `cf_ok`.
        type(j2_plasticity_t) :: base
        real(wp) :: dummy_input(0,0)
        real(wp) :: local_state(nbase_state,1)
        real(wp) :: effective(nmandel,1)
        real(wp) :: effective_tangent(nmandel,nmandel,1)
        real(wp) :: elastic_strain(nmandel)
        real(wp) :: elastic_deviator(nmandel)
        real(wp) :: trial_deviator(nmandel)
        real(wp) :: flow(nmandel)
        real(wp) :: projector(nmandel,nmandel)
        real(wp) :: identity(nmandel)
        real(wp) :: derivative_trial(nmandel)
        real(wp) :: derivative_flow(nmandel)
        real(wp) :: derivative_plastic(nmandel)
        real(wp) :: depsilon(nmandel)
        real(wp) :: effective_stress(nmandel)
        real(wp) :: alpha_old
        real(wp) :: alpha
        real(wp) :: damage_old
        real(wp) :: damage
        real(wp) :: ddamage_dalpha
        real(wp) :: expected_old_damage
        real(wp) :: ignored_derivative
        real(wp) :: gamma
        real(wp) :: delta_damage
        real(wp) :: trace_strain
        real(wp) :: elastic_energy
        real(wp) :: available_energy
        real(wp) :: qtrial
        real(wp) :: denominator
        real(wp) :: dq
        real(wp) :: dgamma
        real(wp) :: ddamage
        real(wp) :: denergy
        real(wp) :: derivative_dissipation
        integer :: ipoint
        integer :: j
        integer :: i
        integer :: local_status(1)
        logical :: ok

        base = j2_plasticity_t(self%bulk_modulus,self%shear_modulus,self%yield_stress,self%hardening_modulus)
        identity = [1.0_wp,1.0_wp,1.0_wp,0.0_wp,0.0_wp,0.0_wp]
        if (present(tangent)) then
            projector = 0.0_wp
            do i = 1,nmandel
                projector(i,i) = 1.0_wp
            end do
            projector(1:3,1:3) = projector(1:3,1:3)-1.0_wp/3.0_wp
        end if

        point_status = cf_invalid_argument
        point_loop: do ipoint = 1,size(x_new,2)
            if (.not. all(ieee_is_finite(x_new(:,ipoint)))) cycle point_loop
            if (.not. all(ieee_is_finite(state_old(:,ipoint)))) cycle point_loop
            alpha_old = state_old(nbase_state,ipoint)
            damage_old = state_old(nstate,ipoint)
            if (alpha_old < 0.0_wp .or. damage_old < 0.0_wp) cycle point_loop
            if (damage_old > self%maximum_damage) cycle point_loop
            call evaluate_damage(self,alpha_old,expected_old_damage,ignored_derivative,ok)
            if (.not. ok) cycle point_loop
            if (abs(damage_old-expected_old_damage) > 256.0_wp*epsilon(1.0_wp)) cycle point_loop

            if (present(tangent)) then
                call base%eval_batch(dummy_input,x_new(:,ipoint:ipoint),state_old(1:nbase_state,ipoint:ipoint), local_state,effective,0.0_wp,effective_tangent,local_status)
            else
                call base%eval_batch(dummy_input,x_new(:,ipoint:ipoint),state_old(1:nbase_state,ipoint:ipoint), local_state,effective,0.0_wp,point_status=local_status)
            end if
            if (local_status(1) /= cf_ok) then
                point_status(ipoint) = local_status(1)
                cycle point_loop
            end if

            alpha = local_state(nbase_state,1)
            gamma = alpha-alpha_old
            if (gamma < 0.0_wp .or. .not. ieee_is_finite(gamma)) cycle point_loop
            call evaluate_damage(self,alpha,damage,ddamage_dalpha,ok)
            if (.not. ok) cycle point_loop
            delta_damage = damage-damage_old
            if (delta_damage < 0.0_wp) cycle point_loop

            elastic_strain = x_new(:,ipoint)-local_state(1:nmandel,1)
            trace_strain = sum(elastic_strain(1:3))
            elastic_deviator = elastic_strain-trace_strain*identity/3.0_wp
            elastic_energy = 0.5_wp*self%bulk_modulus*trace_strain**2 + self%shear_modulus*dot_product(elastic_deviator,elastic_deviator)
            available_energy = elastic_energy + 0.5_wp*self%hardening_modulus*alpha**2
            effective_stress = effective(:,1)
            response(1:nmandel,ipoint) = (1.0_wp-damage)*effective_stress
            response(nbase_state,ipoint) = (1.0_wp-damage)*available_energy
            response(noutput,ipoint) = (1.0_wp-damage)*self%yield_stress*gamma + available_energy*delta_damage
            if (.not. all(ieee_is_finite(response(:,ipoint)))) then
                point_status(ipoint) = cf_evaluation_failed
                cycle point_loop
            end if
            if (response(nbase_state,ipoint) < 0.0_wp .or. response(noutput,ipoint) < 0.0_wp) then
                point_status(ipoint) = cf_evaluation_failed
                cycle point_loop
            end if

            if (present(tangent)) then
                trial_deviator = x_new(:,ipoint)-state_old(1:nmandel,ipoint)
                trial_deviator(1:3) = trial_deviator(1:3)-sum(x_new(1:3,ipoint))/3.0_wp
                trial_deviator = 2.0_wp*self%shear_modulus*trial_deviator
                qtrial = sqrt(1.5_wp*dot_product(trial_deviator,trial_deviator))
                denominator = 3.0_wp*self%shear_modulus+self%hardening_modulus
                flow = 0.0_wp
                if (gamma > 0.0_wp) then
                    if (qtrial <= 0.0_wp .or. .not. ieee_is_finite(qtrial)) then
                        point_status(ipoint) = cf_evaluation_failed
                        cycle point_loop
                    end if
                    flow = 1.5_wp*trial_deviator/qtrial
                end if
                do j = 1,nmandel
                    dgamma = 0.0_wp
                    derivative_plastic = 0.0_wp
                    if (gamma > 0.0_wp) then
                        derivative_trial = 2.0_wp*self%shear_modulus*projector(:,j)
                        dq = dot_product(flow,derivative_trial)
                        dgamma = dq/denominator
                        derivative_flow = (1.5_wp/qtrial)*(derivative_trial-trial_deviator*dq/qtrial)
                        derivative_plastic = dgamma*flow+gamma*derivative_flow
                    end if
                    ddamage = ddamage_dalpha*dgamma
                    depsilon = -derivative_plastic
                    depsilon(j) = depsilon(j)+1.0_wp
                    denergy = dot_product(effective_stress,depsilon)+self%hardening_modulus*alpha*dgamma
                    tangent(1:nmandel,j,ipoint) = (1.0_wp-damage)*effective_tangent(:,j,1)-effective_stress*ddamage
                    tangent(nbase_state,j,ipoint) = (1.0_wp-damage)*denergy-available_energy*ddamage
                    derivative_dissipation = (1.0_wp-damage)*self%yield_stress*dgamma - self%yield_stress*gamma*ddamage+delta_damage*denergy+available_energy*ddamage
                    tangent(noutput,j,ipoint) = derivative_dissipation
                end do
                if (.not. all(ieee_is_finite(tangent(:,:,ipoint)))) then
                    point_status(ipoint) = cf_evaluation_failed
                    cycle point_loop
                end if
            end if

            state_trial(1:nbase_state,ipoint) = local_state(:,1)
            state_trial(nstate,ipoint) = damage
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_ductile_damage

    pure subroutine evaluate_damage(self,alpha,damage,derivative,ok,fraction,decay)
        !! Explicit monotone saturation. Guard multiplication overflow.
        class(ductile_damage_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        real(wp), intent(in) :: alpha
            !! Accumulated equivalent plastic strain governing damage evolution.
        real(wp), intent(out) :: damage
            !! Scalar damage predicted by the hardening-dependent evolution law.
        real(wp), intent(out) :: derivative
            !! Local damage derivative \(\mathrm{d}D/\mathrm{d}\alpha\).
        logical, intent(out) :: ok
            !! True if the local computation satisfies its admissibility conditions.
        real(wp), intent(out), optional :: fraction
            !! Cancellation-free damage fraction before multiplication by the damage cap.
        real(wp), intent(out), optional :: decay
            !! Exponential decay factor shared with analytical parameter sensitivities.
        real(wp) :: elapsed
        real(wp) :: argument
        real(wp) :: exponential
        real(wp) :: activation

        damage = 0.0_wp
        derivative = 0.0_wp
        if (present(fraction)) fraction = 0.0_wp
        if (present(decay)) decay = 1.0_wp
        ok = ieee_is_finite(alpha)
        if (.not. ok) return
        if (alpha <= self%damage_threshold) return
        elapsed = alpha-self%damage_threshold
        if (.not. ieee_is_finite(elapsed)) then
            ok = .false.
            return
        end if
        argument = huge(1.0_wp)
        if (elapsed > 1.0_wp) then
            if (self%damage_rate <= huge(1.0_wp)/elapsed) argument = self%damage_rate*elapsed
        else
            argument = self%damage_rate*elapsed
        end if
        exponential = exp(-argument)
        if (argument <= 1.0_wp) then
            activation = 2.0_wp*exp(-0.5_wp*argument)*sinh(0.5_wp*argument)
        else
            activation = 1.0_wp-exponential
        end if
        damage = self%maximum_damage*activation
        derivative = self%maximum_damage*self%damage_rate*exponential
        if (present(fraction)) fraction = activation
        if (present(decay)) decay = exponential
        ok = ieee_is_finite(damage) .and. ieee_is_finite(derivative)
    end subroutine evaluate_damage


    subroutine initialize_damage_sensitivity(self,initial_input,initial_state,selected_parameters,sensitivity,status)
        !! Initialize exact parameter sensitivities of the reference material history.
        class(ductile_damage_t),intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        real(wp),intent(in) :: initial_input(:,:)
            !! Initial packed material inputs used to establish the reference history.
        real(wp),intent(in) :: initial_state(:,:)
            !! Initial real-valued history `(nstate,npoints)`.
        integer,intent(in) :: selected_parameters(:)
            !! One-based IDs of selected model parameters, in sensitivity-axis order.
        real(wp),intent(out) :: sensitivity(:,:,:)
            !! Reference-state derivatives `(nstate,nselected,npoints)` initialized consistently with the model.
        integer,intent(out) :: status
            !! Operation status. `cf_ok` indicates a valid model or successful initialization.

        sensitivity=0.0_wp
        status=cf_ok
    end subroutine initialize_damage_sensitivity

    subroutine sensitivity_damage(self,x_old,x_new,state_old,state_trial,selected_parameters, sensitivity_old,sensitivity_trial,response,dsdp,dt,tangent,point_status)
        !! Combine exact \(J_2\) path sensitivities with analytic energy and
        !! damage derivatives. The base \(J_2\) model is reused, never
        !! differentiated by parameter perturbation.
        class(ductile_damage_t),intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        real(wp),intent(in) :: x_old(:,:)
            !! Previously accepted packed input `(ninput,npoints)` if required, otherwise an empty array.
        real(wp),intent(in) :: x_new(:,:)
            !! Prescribed `strain.small` input `(ninput,npoints)`.
        real(wp),intent(in) :: state_old(:,:)
            !! Committed real-valued state `(nstate,npoints)`.
        integer,intent(in) :: selected_parameters(:)
            !! One-based IDs of selected model parameters, in sensitivity-axis order.
        real(wp),intent(in) :: sensitivity_old(:,:,:)
            !! Committed real-state derivatives `(nstate,nselected,npoints)`.
        real(wp),intent(out) :: state_trial(:,:)
            !! Trial real-valued state `(nstate,npoints)` returned by this evaluation.
        real(wp),intent(out) :: sensitivity_trial(:,:,:)
            !! Trial real-state derivatives `(nstate,nselected,npoints)`.
        real(wp),intent(out) :: response(:,:)
            !! Packed `stress.cauchy` and `energy.helmholtz` and `energy.dissipation_increment` response `(noutput,npoints)` for the evaluated material points.
        real(wp),intent(out) :: dsdp(:,:,:)
            !! Exact total response derivatives `(noutput,nselected,npoints)` at fixed prescribed input history.
        real(wp),intent(in) :: dt
            !! Prescribed time increment. It is positive when the material requires time integration.
        real(wp),intent(out),optional :: tangent(:,:,:)
            !! Optional exact discrete response Jacobian `(noutput,ninput,npoints)` at fixed committed history.
        integer,intent(inout) :: point_status(:)
            !! Per-point status values; successful evaluations must set each entry to `cf_ok`.
        type(j2_plasticity_t) :: base
        integer :: base_selected(size(selected_parameters))
        integer :: local_status(1)
        real(wp) :: dummy_input(0,0)
        real(wp) :: base_state(nbase_state,1)
        real(wp) :: base_stress(nmandel,1)
        real(wp) :: base_dsdp(nmandel,size(selected_parameters),1)
        real(wp) :: base_sensitivity_old(nbase_state,size(selected_parameters),1)
        real(wp) :: base_sensitivity_trial(nbase_state,size(selected_parameters),1)
        real(wp) :: identity(nmandel)
        real(wp) :: elastic_strain(nmandel)
        real(wp) :: elastic_dev(nmandel)
        real(wp) :: dep(nmandel)
        real(wp) :: dev_dep(nmandel)
        real(wp) :: dsigma(nmandel)
        real(wp) :: alpha
        real(wp) :: alpha_old
        real(wp) :: damage
        real(wp) :: damage_old
        real(wp) :: delta_alpha
        real(wp) :: delta_damage
        real(wp) :: trace_e
        real(wp) :: energy
        real(wp) :: damage_alpha
        real(wp) :: exp_term
        real(wp) :: damage_fraction
        real(wp) :: offset
        real(wp) :: strain_scale
        real(wp) :: d_bulk
        real(wp) :: d_shear
        real(wp) :: d_yield
        real(wp) :: d_hardening
        real(wp) :: ddamage
        real(wp) :: dalphap
        real(wp) :: dalpha_old
        real(wp) :: dtrace
        real(wp) :: denergy
        real(wp) :: ddissipation
        integer :: ipoint
        integer :: iparam
        integer :: id
        logical :: ok

        call evaluate_ductile_damage(self,x_old,x_new,state_old,state_trial,response,dt,tangent,point_status)
        base=j2_plasticity_t(self%bulk_modulus,self%shear_modulus,self%yield_stress,self%hardening_modulus)
        identity=[1.0_wp,1.0_wp,1.0_wp,0.0_wp,0.0_wp,0.0_wp]
        sensitivity_trial=sensitivity_old
        dsdp=0.0_wp
        do iparam=1,size(selected_parameters)
            id=selected_parameters(iparam)
            if (id>=1.and.id<=4) then
                base_selected(iparam)=id
            else if (id>=5.and.id<=7) then
                ! A dummy bulk-modulus column is discarded for damage-only
                ! parameters. Their effective J2 state derivative is zero.
                base_selected(iparam)=1
            else
                point_status=cf_invalid_argument
                return
            end if
        end do

        point_loop: do ipoint=1,size(x_new,2)
            if (point_status(ipoint)/=cf_ok) cycle point_loop
            alpha=state_trial(nbase_state,ipoint)
            alpha_old=state_old(nbase_state,ipoint)
            damage=state_trial(nstate,ipoint)
            damage_old=state_old(nstate,ipoint)
            delta_alpha=alpha-alpha_old
            delta_damage=damage-damage_old
            if (self%maximum_damage>0.0_wp) then
                strain_scale=max(alpha,self%damage_threshold)
                if (strain_scale<=0.0_wp) then
                    point_status(ipoint)=cf_nondifferentiable
                    cycle point_loop
                end if
                if (abs(alpha-self%damage_threshold)/strain_scale<=64.0_wp*epsilon(1.0_wp)) then
                    point_status(ipoint)=cf_nondifferentiable
                    cycle point_loop
                end if
            end if

            base_sensitivity_old=sensitivity_old(1:nbase_state,:,ipoint:ipoint)
            do iparam=1,size(selected_parameters)
                if (selected_parameters(iparam)>4) base_sensitivity_old(:,iparam,:)=0.0_wp
            end do
            local_status=cf_evaluation_failed
            call base%eval_sensitivity_batch(dummy_input,x_new(:,ipoint:ipoint), state_old(1:nbase_state,ipoint:ipoint),base_state,base_selected, base_sensitivity_old,base_sensitivity_trial,base_stress,base_dsdp,0.0_wp, point_status=local_status)
            if (local_status(1)/=cf_ok) then
                point_status(ipoint)=local_status(1)
                cycle point_loop
            end if

            elastic_strain=x_new(:,ipoint)-state_trial(1:nmandel,ipoint)
            trace_e=sum(elastic_strain(1:3))
            elastic_dev=elastic_strain-trace_e*identity/3.0_wp
            energy=0.5_wp*self%bulk_modulus*trace_e**2+ self%shear_modulus*dot_product(elastic_dev,elastic_dev)+ 0.5_wp*self%hardening_modulus*alpha**2
            call evaluate_damage(self,alpha,damage,damage_alpha,ok,fraction=damage_fraction,decay=exp_term)
            if (.not. ok) then
                point_status(ipoint)=cf_evaluation_failed
                cycle point_loop
            end if
            offset=max(0.0_wp,alpha-self%damage_threshold)
            parameter_loop: do iparam=1,size(selected_parameters)
                d_bulk=0.0_wp
                d_shear=0.0_wp
                d_yield=0.0_wp
                d_hardening=0.0_wp
                id=selected_parameters(iparam)
                select case(id)
                case(1)
                    d_bulk=1.0_wp
                case(2)
                    d_shear=1.0_wp
                case(3)
                    d_yield=1.0_wp
                case(4)
                    d_hardening=1.0_wp
                case(5:7)
                case default
                    point_status(ipoint)=cf_invalid_argument
                    exit parameter_loop
                end select
                dalphap=0.0_wp
                dep=0.0_wp
                dsigma=0.0_wp
                if (id<=4) then
                    sensitivity_trial(1:nbase_state,iparam,ipoint)=base_sensitivity_trial(:,iparam,1)
                    dalphap=base_sensitivity_trial(nbase_state,iparam,1)
                    dep=base_sensitivity_trial(1:nmandel,iparam,1)
                    dsigma=base_dsdp(:,iparam,1)
                else
                    sensitivity_trial(1:nbase_state,iparam,ipoint)=0.0_wp
                end if
                dalpha_old=sensitivity_old(nbase_state,iparam,ipoint)
                ddamage=damage_alpha*dalphap
                if (alpha>self%damage_threshold) then
                    select case(id)
                    case(5)
                        ddamage=ddamage-self%maximum_damage*self%damage_rate*exp_term
                    case(6)
                        ddamage=ddamage+self%maximum_damage*offset*exp_term
                    case(7)
                        ddamage=ddamage+damage_fraction
                    case default
                    end select
                end if
                sensitivity_trial(nstate,iparam,ipoint)=ddamage
                dtrace=-sum(dep(1:3))
                dev_dep=dep-sum(dep(1:3))*identity/3.0_wp
                denergy=0.5_wp*d_bulk*trace_e**2+ self%bulk_modulus*trace_e*dtrace+ d_shear*dot_product(elastic_dev,elastic_dev)- 2.0_wp*self%shear_modulus*dot_product(elastic_dev,dev_dep)+ 0.5_wp*d_hardening*alpha**2+ self%hardening_modulus*alpha*dalphap
                dsdp(1:nmandel,iparam,ipoint)=(1.0_wp-damage)*dsigma-base_stress(:,1)*ddamage
                dsdp(nbase_state,iparam,ipoint)=(1.0_wp-damage)*denergy-energy*ddamage
                ddissipation=(1.0_wp-damage)*(d_yield*delta_alpha+ self%yield_stress*(dalphap-dalpha_old)) - ddamage*self%yield_stress*delta_alpha+ denergy*delta_damage+ energy*(ddamage-sensitivity_old(nstate,iparam,ipoint))
                dsdp(noutput,iparam,ipoint)=ddissipation
            end do parameter_loop
            if (point_status(ipoint)/=cf_ok) cycle point_loop
            if (.not. all(ieee_is_finite(dsdp(:,:,ipoint))) .or. .not. all(ieee_is_finite(sensitivity_trial(:,:,ipoint)))) then
                point_status(ipoint)=cf_evaluation_failed
            end if
        end do point_loop
    end subroutine sensitivity_damage

end module constituf_ductile_damage

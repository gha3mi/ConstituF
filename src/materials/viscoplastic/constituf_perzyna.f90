! SPDX-License-Identifier: Apache-2.0

module constituf_perzyna
    !! Three-dimensional isothermal small-strain \(J_2\) Perzyna viscoplasticity.
    !!
    !! This independently derived *linear-overstress specialization* of
    !! Perzyna flow has additive strain split and an isotropic elastic
    !! stiffness \(K>0,G>0\). With current strain \(\varepsilon\),
    !! committed plastic strain \(\varepsilon^p_n\) and accumulated
    !! equivalent plastic strain \(\alpha_n\), define
    !! \[\boldsymbol{s}_{\rm tr}=2G\,\mathrm{dev}
    !! (\boldsymbol\varepsilon-\boldsymbol\varepsilon^p_n),\quad
    !! q_{\rm tr}=\sqrt{3/2}\|\boldsymbol{s}_{\rm tr}\|,\quad
    !! f_{\rm tr}=q_{\rm tr}-(Y_0+H\alpha_n).\]
    !! The positive overstress law uses viscosity \(\eta>0\) with dimensions
    !! of stress times time, using the same time unit as `dt`:
    !! \[\dot{\gamma}=\frac{\langle
    !! q-(Y_0+H\alpha)\rangle}{\eta},\qquad
    !! \dot{\boldsymbol\varepsilon}^p=\dot\gamma
    !! \frac{3\boldsymbol{s}}{2q},\qquad \dot\alpha=\dot\gamma.\]
    !! Backward Euler over \(\delta t>0\) gives an *exact scalar* solve:
    !! \[\Delta\gamma =
    !! \frac{\delta t\,\langle f_{\rm tr}\rangle}
    !! {\eta+\delta t(3G+H)},\quad
    !! \boldsymbol n=\frac{3\boldsymbol{s}_{\rm tr}}{2q_{\rm tr}},\]
    !! \[\boldsymbol\varepsilon^p_{n+1}
    !! =\boldsymbol\varepsilon^p_n+\Delta\gamma\boldsymbol n,\quad
    !! \alpha_{n+1}=\alpha_n+\Delta\gamma,\quad
    !! \boldsymbol\sigma_{n+1}=K\,\mathrm{tr}
    !! (\boldsymbol\varepsilon_{n+1})\boldsymbol I
    !! +\boldsymbol{s}_{\rm tr}-2G\Delta\gamma\boldsymbol n.\]
    !! The yield condition may have *positive* overstress after a time step.
    !! Enforcing \(f=0\) would instead describe rate-independent plasticity.
    !! For \(\Delta\gamma>0\), viscous dissipation is
    !! \(\Delta\gamma\,(q_{n+1}-Y_0-H\alpha_{n+1})
    !! =\eta\,\Delta\gamma^2/\delta t\geq0\).
    !!
    !! The model also offers strain-increment-controlled subdivision of an
    !! accepted global step. Given previous committed strain
    !! \(\varepsilon_n\), split the prescribed line segment to
    !! \(\varepsilon_{n+1}\) into \(N\) equal backward-Euler steps:
    !! \[\varepsilon_k=\varepsilon_n+\frac{k}{N}
    !! (\varepsilon_{n+1}-\varepsilon_n),\quad \delta t=\Delta t/N.\]
    !! A positive `maximum_substep_strain` caps the maximum packed Mandel
    !! strain-component increment in each substep, up to `max_substeps`.
    !! Zero disables subdivision. This is *increment-controlled*
    !! substepping, not a local-error-estimated adaptive integrator.
    !! Exceeding `max_substeps` returns `cf_evaluation_failed`. The committed
    !! initial strain must be supplied via `material_t%init(initial_input=...)`.
    !!
    !! With a fixed active-set sequence and substep count, the exact tangent
    !! with respect to the **global endpoint** strain is propagated through
    !! each substep by the chain rule. Let \(P_{\rm dev}\) be the Mandel
    !! deviatoric projector, \(a_k=k/N\) and \(D\) denote
    !! differentiation with respect to \(\varepsilon_{n+1}\). Then
    !! \[D\boldsymbol{s}_{\rm tr}
    !! =2G(a_k P_{\rm dev}-D\boldsymbol\varepsilon^p_{k-1}),
    !! \qquad Dq_{\rm tr}=\boldsymbol n^T
    !! D\boldsymbol{s}_{\rm tr},\]
    !! \[D\Delta\gamma=\frac{\delta t}{
    !! \eta+\delta t(3G+H)}
    !! (Dq_{\rm tr}-H D\alpha_{k-1}),\quad
    !! D\boldsymbol n=\frac{3}{2q_{\rm tr}}
    !! \left(D\boldsymbol{s}_{\rm tr}-
    !! \frac{\boldsymbol{s}_{\rm tr}}{q_{\rm tr}}Dq_{\rm tr}\right).\]
    !! Update plastic-history derivatives before proceeding to the next
    !! substep. The final tangent is \(D\boldsymbol\sigma_{n+1}\).
    !! No finite differences occur in this model's production tangent.
    !! At active-set or subdivision-count changes the map is piecewise
    !! differentiable. A classical derivative may not exist.
    !!
    !! Input `strain.small` and output `stress.cauchy` are symmetric
    !! spatial tensors. Each such tensor \(\boldsymbol{A}\) is packed as
    !! \[\operatorname{mandel}(\boldsymbol{A})=
    !! [A_{11},A_{22},A_{33},\sqrt{2}A_{23},
    !! \sqrt{2}A_{13},\sqrt{2}A_{12}]^{\mathsf T}.\]
    !! Strain is dimensionless and stress uses the chosen stress unit.
    !! The model owns seven real history components:
    !! six plastic-strain components and accumulated \(\alpha\). \(K,G,Y_0,H\) are
    !! in the chosen stress unit. Viscosity \(\eta\) has dimensions of stress
    !! times time. \(H\geq0\), while \(K,G,Y_0,\eta>0\).
    !! Requires positive `dt` and an initial input. Exact parameter-path
    !! sensitivities differentiate each discrete backward-Euler substep,
    !! including the entire accepted plastic-history path. The temperature
    !! and material
    !! orientation are fixed. No finite-strain or kinematic hardening.
    !!
    !! References: P. Perzyna, "The Constitutive Equations for Rate
    !! Sensitive Plastic Materials", Quarterly of Applied Mathematics
    !! **20** (1963), 321--332, DOI: 10.1090/qam/144536.
    !! Equations (8), (46) and the \(S=1\) specialization in Section 5
    !! give the normalized linear-overstress law. For perfect plasticity,
    !! \(k=Y_0/\sqrt3\) and the original rate coefficient
    !! \(\gamma_{\rm ref}\) maps to
    !! \(\eta=\sqrt3\,Y_0/(2\gamma_{\rm ref})\).
    !!
    !! P. Perzyna, "Fundamental Problems in Viscoplasticity",
    !! Advances in Applied Mechanics **9** (1966), pp. 243--377,
    !! DOI: 10.1016/S0065-2156(08)70009-7. The exponent-one overstress
    !! law with constant dimensional viscosity, linear hardening and
    !! strain-limited substepping are the stated specialization and
    !! numerical adaptation, not a claim of identical published algorithm.
    !! Code provenance: original from these equations. No third-party code.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument, cf_evaluation_failed, cf_nondifferentiable
    implicit none (type, external)
    private

    integer, parameter :: nmandel = 6
    integer, parameter :: nstate = 7
    public :: perzyna_t

    type, extends(material_model_t) :: perzyna_t
        !! Linear-overstress \(J_2\) viscoplasticity with optional subdivision.
        !! Construct as `perzyna_t(K,G,Y0,H,eta,max_substep_strain,max_substeps)`.
        real(wp) :: bulk_modulus = 0.0_wp
            !! Isotropic bulk modulus \(K>0\) in the chosen stress unit.
        real(wp) :: shear_modulus = 0.0_wp
            !! Elastic shear modulus \(G>0\) in the chosen stress unit.
        real(wp) :: yield_stress = 0.0_wp
            !! Initial equivalent yield stress \(Y_0>0\) in the chosen stress unit.
        real(wp) :: hardening_modulus = 0.0_wp
            !! Linear isotropic hardening \(H\geq0\) in the chosen stress unit.
        real(wp) :: viscosity = 0.0_wp
            !! Linear overstress viscosity \(\eta>0\) in stress times time.
        real(wp) :: maximum_substep_strain = 0.0_wp
            !! Maximum absolute Mandel strain component per substep,
            !! zero disables subdivision.
        integer :: max_substeps = 128
            !! Upper bound \(N_{\max}\geq1\) on incremental subdivision.
    contains
        procedure :: describe => describe_perzyna
            !! Declare material fields, parameters and supported operations.
        procedure :: eval_batch => evaluate_perzyna
            !! Evaluate the batched material response, optional tangent and trial state.
        procedure :: initialize_sensitivity => initialize_perzyna_sensitivity
            !! Initialize exact derivatives of the material's reference state.
        procedure :: eval_sensitivity_batch => sensitivity_perzyna
            !! Evaluate analytical parameter derivatives through the discrete material update.
    end type perzyna_t

contains

    subroutine describe_perzyna(self, spec, status)
        !! Validate material parameters and register `strain.small` input and `stress.cauchy` response fields.
        class(perzyna_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema populated with fields, physical metadata and model capabilities.
        integer, intent(out) :: status
            !! Operation status code. `cf_ok` indicates success.

        status = cf_invalid_argument
        if (.not. all(ieee_is_finite([self%bulk_modulus,self%shear_modulus,self%yield_stress, self%hardening_modulus,self%viscosity,self%maximum_substep_strain]))) return
        if (self%bulk_modulus <= 0.0_wp .or. self%shear_modulus <= 0.0_wp) return
        if (self%yield_stress <= 0.0_wp .or. self%hardening_modulus < 0.0_wp .or. self%viscosity <= 0.0_wp) return
        if (self%maximum_substep_strain < 0.0_wp .or. self%max_substeps < 1) return
        if (.not. ieee_is_finite(3.0_wp*self%shear_modulus+self%hardening_modulus)) return
        if (.not. ieee_is_finite(self%bulk_modulus+(4.0_wp/3.0_wp)*self%shear_modulus)) return

        call spec%add_input("strain.small", nmandel, "mandel", "small_strain", "spatial", [3,3])
        call spec%add_output("stress.cauchy", nmandel, "mandel", "cauchy_stress", "spatial", [3,3])
        call spec%add_state("plastic.strain", nmandel, "mandel", "plastic_strain", "spatial", [3,3])
        call spec%add_state("plastic.accumulated", 1, "scalar", "equivalent_plastic_strain", "material")
        call spec%add_parameter("bulk_modulus")
        call spec%add_parameter("shear_modulus")
        call spec%add_parameter("yield_stress")
        call spec%add_parameter("hardening_modulus")
        call spec%add_parameter("viscosity")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        call spec%set_requires_dt(.true.)
        call spec%set_requires_initial_input(.true.)
        call spec%set_requires_previous_input(.true.)
        status = cf_ok
    end subroutine describe_perzyna

    subroutine evaluate_perzyna(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Piecewise-exact analytical tangent for every accepted substep.
        !!
        !! Each integration point starts from its committed old strain and
        !! history. All local arrays are fixed-size and there is no I/O,
        !! heap allocation or dynamic field lookup in the batch kernel.
        class(perzyna_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previous packed inputs `(ninput,npoints)` when required by the model, otherwise a zero-sized array.
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
        real(wp) :: projector(nmandel,nmandel)
        real(wp) :: strain(nmandel)
        real(wp) :: deviation(nmandel)
        real(wp) :: trial_stress(nmandel)
        real(wp) :: normal(nmandel)
        real(wp) :: plastic_strain(nmandel)
        real(wp) :: current_stress(nmandel)
        real(wp) :: delta(nmandel)
        real(wp) :: derivative_plastic(nmandel,nmandel)
        real(wp) :: derivative_alpha(nmandel)
        real(wp) :: derivative_stress(nmandel,nmandel)
        real(wp) :: derivative_trial(nmandel)
        real(wp) :: derivative_normal(nmandel)
        real(wp) :: plastic_increment(nmandel)
        real(wp) :: alpha
        real(wp) :: trace_strain
        real(wp) :: equivalent_trial
        real(wp) :: threshold
        real(wp) :: yield_function
        real(wp) :: sub_dt
        real(wp) :: denominator
        real(wp) :: coefficient
        real(wp) :: gamma
        real(wp) :: fraction
        real(wp) :: maximum_change
        real(wp) :: derivative_equivalent
        real(wp) :: derivative_gamma
        integer :: ipoint
        integer :: i
        integer :: j
        integer :: istep
        integer :: nsub
        logical :: valid

        projector = 0.0_wp
        do i = 1, nmandel
            projector(i,i) = 1.0_wp
        end do
        projector(1:3,1:3) = projector(1:3,1:3)-1.0_wp/3.0_wp

        point_status = cf_invalid_argument
        point_loop: do ipoint = 1, size(x_new,2)
            if (.not. all(ieee_is_finite(x_new(:,ipoint)))) cycle point_loop
            if (.not. all(ieee_is_finite(x_old(:,ipoint)))) cycle point_loop
            if (.not. all(ieee_is_finite(state_old(:,ipoint)))) cycle point_loop
            alpha = state_old(nstate,ipoint)
            if (alpha < 0.0_wp) cycle point_loop
            delta = x_new(:,ipoint)-x_old(:,ipoint)
            if (.not. all(ieee_is_finite(delta))) cycle point_loop
            nsub = 1
            if (self%maximum_substep_strain > 0.0_wp) then
                maximum_change = maxval(abs(delta))
                if (maximum_change > real(self%max_substeps,wp)*self%maximum_substep_strain) then
                    point_status(ipoint) = cf_evaluation_failed
                    cycle point_loop
                end if
                nsub = max(1,ceiling(maximum_change/self%maximum_substep_strain))
            end if

            sub_dt = dt/real(nsub,wp)
            ! Use equivalent forms that avoid multiplying a long time step by stiffness.
            if (sub_dt >= 1.0_wp) then
                denominator = self%viscosity/sub_dt+3.0_wp*self%shear_modulus+self%hardening_modulus
            else
                denominator = self%viscosity+sub_dt*(3.0_wp*self%shear_modulus+self%hardening_modulus)
            end if
            if (.not. ieee_is_finite(denominator)) cycle point_loop
            if (denominator <= 0.0_wp) cycle point_loop
            if (sub_dt >= 1.0_wp) then
                coefficient = 1.0_wp/denominator
            else
                coefficient = sub_dt/denominator
            end if
            if (.not. ieee_is_finite(coefficient)) cycle point_loop

            plastic_strain = state_old(1:nmandel,ipoint)
            if (present(tangent)) then
                derivative_plastic = 0.0_wp
                derivative_alpha = 0.0_wp
            end if
            valid = .true.
            substep_loop: do istep = 1, nsub
                fraction = real(istep,wp)/real(nsub,wp)
                strain = x_old(:,ipoint)+fraction*delta
                trace_strain = sum(strain(1:3))
                deviation = strain
                deviation(1:3) = deviation(1:3)-trace_strain/3.0_wp
                trial_stress = 2.0_wp*self%shear_modulus*(deviation-plastic_strain)
                equivalent_trial = sqrt(1.5_wp*dot_product(trial_stress,trial_stress))
                threshold = self%yield_stress+self%hardening_modulus*alpha
                if (.not. all(ieee_is_finite(trial_stress)) .or. .not. ieee_is_finite(equivalent_trial) .or. .not. ieee_is_finite(threshold)) then
                    valid = .false.
                    exit substep_loop
                end if
                if (threshold <= 0.0_wp) then
                    valid = .false.
                    exit substep_loop
                end if
                yield_function = equivalent_trial-threshold
                gamma = 0.0_wp
                normal = 0.0_wp
                if (yield_function > 0.0_wp) then
                    gamma = coefficient*yield_function
                    normal = (1.5_wp/equivalent_trial)*trial_stress
                end if
                if (.not. ieee_is_finite(gamma) .or. .not. all(ieee_is_finite(normal))) then
                    valid = .false.
                    exit substep_loop
                end if

                current_stress = trial_stress-2.0_wp*self%shear_modulus*gamma*normal
                current_stress(1:3) = current_stress(1:3)+self%bulk_modulus*trace_strain
                if (.not. all(ieee_is_finite(current_stress))) then
                    valid = .false.
                    exit substep_loop
                end if
                if (present(tangent)) then
                    do j = 1, nmandel
                        derivative_trial = 2.0_wp*self%shear_modulus*(fraction*projector(:,j)-derivative_plastic(:,j))
                        derivative_stress(:,j) = derivative_trial
                        if (j <= 3) derivative_stress(1:3,j) = derivative_stress(1:3,j) + self%bulk_modulus*fraction
                        if (yield_function > 0.0_wp) then
                            derivative_equivalent = dot_product(normal,derivative_trial)
                            derivative_gamma = coefficient*(derivative_equivalent-self%hardening_modulus*derivative_alpha(j))
                            derivative_normal = (1.5_wp/equivalent_trial)* (derivative_trial-trial_stress*derivative_equivalent/equivalent_trial)
                            derivative_stress(:,j) = derivative_stress(:,j)-2.0_wp*self%shear_modulus* (derivative_gamma*normal+gamma*derivative_normal)
                            derivative_plastic(:,j) = derivative_plastic(:,j)+derivative_gamma*normal+gamma*derivative_normal
                            derivative_alpha(j) = derivative_alpha(j)+derivative_gamma
                        end if
                    end do
                end if

                plastic_increment = gamma*normal
                plastic_strain = plastic_strain+plastic_increment
                alpha = alpha+gamma
                if (.not. all(ieee_is_finite(plastic_strain)) .or. .not. ieee_is_finite(alpha)) then
                    valid = .false.
                    exit substep_loop
                end if
            end do substep_loop
            if (.not. valid) then
                point_status(ipoint) = cf_evaluation_failed
                cycle point_loop
            end if

            state_trial(1:nmandel,ipoint) = plastic_strain
            state_trial(nstate,ipoint) = alpha
            response(:,ipoint) = current_stress
            if (present(tangent)) then
                if (.not. all(ieee_is_finite(derivative_stress))) then
                    point_status(ipoint) = cf_evaluation_failed
                    cycle point_loop
                end if
                tangent(:,:,ipoint) = derivative_stress
            end if
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_perzyna


    subroutine initialize_perzyna_sensitivity(self,initial_input,initial_state,selected_parameters,sensitivity,status)
        !! Virgin state and its derivatives are independent of material parameters.
        class(perzyna_t),intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        real(wp),intent(in) :: initial_input(:,:)
            !! Prescribed initial input used to initialize the model's reference state.
        real(wp),intent(in) :: initial_state(:,:)
            !! Initial real-valued history `(nstate,npoints)`.
        integer,intent(in) :: selected_parameters(:)
            !! One-based IDs of scalar parameters in the requested derivative-axis order.
        real(wp),intent(out) :: sensitivity(:,:,:)
            !! Exact initial real-history derivatives `(nstate,nselected,npoints)`.
        integer,intent(out) :: status
            !! Operation status code. `cf_ok` indicates success.

        sensitivity=0.0_wp
        status=cf_ok
    end subroutine initialize_perzyna_sensitivity

    subroutine sensitivity_perzyna(self,x_old,x_new,state_old,state_trial,selected_parameters, sensitivity_old,sensitivity_trial,response,dsdp,dt,tangent,point_status)
        !! Direct differentiation of every backward-Euler substep. No
        !! parameter perturbation, finite-difference or numerical Jacobian.
        class(perzyna_t),intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        real(wp),intent(in) :: x_old(:,:)
            !! Previous packed inputs `(ninput,npoints)` when required by the model, otherwise a zero-sized array.
        real(wp),intent(in) :: x_new(:,:)
            !! Prescribed `strain.small` inputs `(ninput,npoints)` at the current time step.
        real(wp),intent(in) :: state_old(:,:)
            !! Committed real-valued history `(nstate,npoints)`.
        integer,intent(in) :: selected_parameters(:)
            !! One-based IDs of scalar parameters in the requested derivative-axis order.
        real(wp),intent(in) :: sensitivity_old(:,:,:)
            !! Committed history derivatives `(nstate,nselected,npoints)`.
        real(wp),intent(out) :: state_trial(:,:)
            !! Updated trial real-valued history `(nstate,npoints)`.
        real(wp),intent(out) :: sensitivity_trial(:,:,:)
            !! Updated trial history derivatives `(nstate,nselected,npoints)`.
        real(wp),intent(out) :: response(:,:)
            !! Packed `stress.cauchy` outputs `(noutput,npoints)`.
        real(wp),intent(out) :: dsdp(:,:,:)
            !! Exact response derivatives `(noutput,nselected,npoints)` at fixed prescribed input histories.
        real(wp),intent(in) :: dt
            !! Prescribed time increment, strictly positive for time-dependent materials.
        real(wp),intent(out),optional :: tangent(:,:,:)
            !! Optional exact discrete input Jacobian `(noutput,ninput,npoints)`.
        integer,intent(inout) :: point_status(:)
            !! Per-point status codes. successful points must be assigned `cf_ok`.
        real(wp) :: identity(nmandel)
        real(wp) :: strain(nmandel)
        real(wp) :: eps_dev(nmandel)
        real(wp) :: plastic_strain(nmandel)
        real(wp) :: dplastic(nmandel)
        real(wp) :: stress_trial(nmandel)
        real(wp) :: dstress_trial(nmandel)
        real(wp) :: normal(nmandel)
        real(wp) :: dnormal(nmandel)
        real(wp) :: dstress(nmandel)
        real(wp) :: d_bulk
        real(wp) :: d_shear
        real(wp) :: d_yield
        real(wp) :: d_hardening
        real(wp) :: d_viscosity
        real(wp) :: alpha
        real(wp) :: dalpha
        real(wp) :: trace_strain
        real(wp) :: threshold
        real(wp) :: dthreshold
        real(wp) :: q
        real(wp) :: dq
        real(wp) :: overstress
        real(wp) :: fraction
        real(wp) :: sub_dt
        real(wp) :: denominator
        real(wp) :: coefficient
        real(wp) :: dcoefficient
        real(wp) :: gamma
        real(wp) :: dgamma
        real(wp) :: maximum_change
        real(wp) :: stress_scale
        real(wp) :: delta(nmandel)
        integer :: ipoint
        integer :: iparam
        integer :: istep
        integer :: nsub
        logical :: smooth

        call evaluate_perzyna(self,x_old,x_new,state_old,state_trial,response,dt,tangent,point_status)
        identity=[1.0_wp,1.0_wp,1.0_wp,0.0_wp,0.0_wp,0.0_wp]
        sensitivity_trial=sensitivity_old
        dsdp=0.0_wp

        point_loop: do ipoint=1,size(x_new,2)
            if (point_status(ipoint)/=cf_ok) cycle point_loop
            delta=x_new(:,ipoint)-x_old(:,ipoint)
            nsub=1
            if (self%maximum_substep_strain>0.0_wp) then
                maximum_change=maxval(abs(delta))
                nsub=max(1,ceiling(maximum_change/self%maximum_substep_strain))
            end if
            sub_dt=dt/real(nsub,wp)
            if (sub_dt>=1.0_wp) then
                denominator=self%viscosity/sub_dt+3.0_wp*self%shear_modulus+self%hardening_modulus
                coefficient=1.0_wp/denominator
            else
                denominator=self%viscosity+sub_dt*(3.0_wp*self%shear_modulus+self%hardening_modulus)
                coefficient=sub_dt/denominator
            end if
            smooth=.true.
            parameter_loop: do iparam=1,size(selected_parameters)
                d_bulk=0.0_wp
                d_shear=0.0_wp
                d_yield=0.0_wp
                d_hardening=0.0_wp
                d_viscosity=0.0_wp
                select case(selected_parameters(iparam))
                case(1)
                    d_bulk=1.0_wp
                case(2)
                    d_shear=1.0_wp
                case(3)
                    d_yield=1.0_wp
                case(4)
                    d_hardening=1.0_wp
                case(5)
                    d_viscosity=1.0_wp
                case default
                    point_status(ipoint)=cf_invalid_argument
                    exit parameter_loop
                end select
                if (sub_dt>=1.0_wp) then
                    dcoefficient=-coefficient*((d_viscosity/sub_dt+3.0_wp*d_shear+d_hardening)/denominator)
                else
                    dcoefficient=-coefficient*((d_viscosity+sub_dt*(3.0_wp*d_shear+d_hardening))/denominator)
                end if
                plastic_strain=state_old(1:nmandel,ipoint)
                dplastic=sensitivity_old(1:nmandel,iparam,ipoint)
                alpha=state_old(nstate,ipoint)
                dalpha=sensitivity_old(nstate,iparam,ipoint)

                substep_loop: do istep=1,nsub
                    fraction=real(istep,wp)/real(nsub,wp)
                    strain=x_old(:,ipoint)+fraction*delta
                    trace_strain=sum(strain(1:3))
                    eps_dev=strain-trace_strain*identity/3.0_wp
                    stress_trial=2.0_wp*self%shear_modulus*(eps_dev-plastic_strain)
                    dstress_trial=2.0_wp*d_shear*(eps_dev-plastic_strain)- 2.0_wp*self%shear_modulus*dplastic
                    q=sqrt(1.5_wp*dot_product(stress_trial,stress_trial))
                    threshold=self%yield_stress+self%hardening_modulus*alpha
                    dthreshold=d_yield+d_hardening*alpha+self%hardening_modulus*dalpha
                    overstress=q-threshold
                    stress_scale=max(q,threshold)
                    if (abs(overstress)/stress_scale<=64.0_wp*epsilon(1.0_wp)) then
                        smooth=.false.
                        exit substep_loop
                    end if
                    gamma=0.0_wp
                    dgamma=0.0_wp
                    normal=0.0_wp
                    dnormal=0.0_wp
                    if (overstress>0.0_wp) then
                        normal=1.5_wp*stress_trial/q
                        dq=dot_product(normal,dstress_trial)
                        gamma=coefficient*overstress
                        dgamma=dcoefficient*overstress+coefficient*(dq-dthreshold)
                        dnormal=(1.5_wp/q)*(dstress_trial-stress_trial*dq/q)
                    end if
                    dstress=dstress_trial+d_bulk*trace_strain*identity- 2.0_wp*d_shear*gamma*normal- 2.0_wp*self%shear_modulus*(dgamma*normal+gamma*dnormal)
                    dplastic=dplastic+dgamma*normal+gamma*dnormal
                    dalpha=dalpha+dgamma
                    plastic_strain=plastic_strain+gamma*normal
                    alpha=alpha+gamma
                end do substep_loop
                if (.not. smooth) exit parameter_loop
                sensitivity_trial(1:nmandel,iparam,ipoint)=dplastic
                sensitivity_trial(nstate,iparam,ipoint)=dalpha
                dsdp(:,iparam,ipoint)=dstress
            end do parameter_loop
            if (.not. smooth) then
                point_status(ipoint)=cf_nondifferentiable
                cycle point_loop
            end if
            if (point_status(ipoint)/=cf_ok) cycle point_loop
            if (.not. all(ieee_is_finite(dsdp(:,:,ipoint))) .or. .not. all(ieee_is_finite(sensitivity_trial(:,:,ipoint)))) then
                point_status(ipoint)=cf_evaluation_failed
            end if
        end do point_loop
    end subroutine sensitivity_perzyna

end module constituf_perzyna

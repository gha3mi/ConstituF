! SPDX-License-Identifier: Apache-2.0

module constituf_maxwell_shear
    !! Three-dimensional infinitesimal-strain generalized Maxwell material.
    !!
    !! The volumetric response and equilibrium deviatoric response are elastic,
    !! with one additional relaxing deviatoric Maxwell branch:
    !! \[\boldsymbol{\sigma}(t)=
    !! K\,\mathrm{tr}(\boldsymbol{\varepsilon}(t))\boldsymbol{I}
    !! +2G_\infty\boldsymbol{e}(t)+\boldsymbol{q}(t),\]
    !! \[\dot{\boldsymbol{q}}+\boldsymbol{q}/\tau
    !! =2G_1\dot{\boldsymbol{e}},\qquad
    !! \boldsymbol{e}=\mathrm{dev}\boldsymbol{\varepsilon}.\]
    !! Here \(K>0\) is the bulk modulus, \(G_\infty\geq 0\) the long-time
    !! shear modulus, \(G_1\geq 0\) the branch shear modulus (with
    !! \(G_\infty+G_1>0\)) and \(\tau>0\) its relaxation time.
    !! Positive \(G_\infty\) gives a viscoelastic solid. Its zero limit
    !! is a Maxwell shear branch with an elastic volumetric response.
    !! Moduli and overstress share the chosen stress unit. Relaxation
    !! time and `dt` share the chosen time unit. Initial branch overstress
    !! is zero at the explicitly required initial strain. This defines an
    !! initially relaxed material, not an arbitrary prescribed prestress.
    !!
    !! Each accepted backward-Euler time step of strictly positive
    !! \(\Delta t\) has the discrete update
    !! \[D=1+\Delta t/\tau,\qquad
    !! \boldsymbol{q}_{n+1}=
    !! \frac{\boldsymbol{q}_n+
    !! 2G_1(\boldsymbol{e}_{n+1}-\boldsymbol{e}_n)}{D}.\]
    !! The exact consistent algorithmic tangent at fixed committed strain
    !! and overstress is
    !! \[\mathbb{C}_{\mathrm{alg}}=
    !! K\boldsymbol{I}\otimes\boldsymbol{I}
    !! +2(G_\infty+G_1/D)\mathbb{I}_{\mathrm{dev}}.\]
    !!
    !! With the model parameters \(p\in\{K,G_\infty,G_1,\tau\}\), the exact
    !! discrete total path sensitivity from accepted history is
    !! \[\frac{\mathrm{d}\boldsymbol{q}_{n+1}}{\mathrm{d}p}
    !! =\frac{\mathrm{d}\boldsymbol{q}_n/\mathrm{d}p+
    !! 2\delta_{pG_1}\Delta\boldsymbol{e}
    !! -\boldsymbol{q}_{n+1}\,\mathrm{d}D/\mathrm{d}p}{D},\]
    !! where \(\mathrm{d}D/\mathrm{d}\tau=-\Delta t/\tau^2\). The stress
    !! sensitivity adds \(\delta_{pK}\mathrm{tr}(\varepsilon)I\)
    !! and \(2\delta_{pG_\infty}e\) to the branch sensitivity.
    !! External strain histories and step lengths are held fixed. No global
    !! solver sensitivities are included.
    !!
    !! Input `strain.small`, response `stress.cauchy` and state
    !! `overstress.deviatoric` are symmetric 3D spatial tensors.
    !! Each such tensor \(\boldsymbol{A}\) is packed as
    !! \[\operatorname{mandel}(\boldsymbol{A})=
    !! [A_{11},A_{22},A_{33},\sqrt{2}A_{23},
    !! \sqrt{2}A_{13},\sqrt{2}A_{12}]^{\mathsf T}.\]
    !! Engineering Voigt shear strains must not be passed as Mandel strains.
    !! Tensile stress/strain is positive. Material is linear viscoelastic,
    !! isothermal, small-rotation and small-strain. No finite-strain,
    !! temperature, creep-plasticity, or volumetric relaxation is modeled.
    !!
    !! References: R. M. Christensen, *Theory of Viscoelasticity: An
    !! Introduction*, 2nd ed., Academic Press, 1982,
    !! ISBN 978-0-12-174252-2, DOI: 10.1016/B978-0-12-174252-2.X5001-7.
    !! FrontISTR documentation, "Generalized Maxwell Model", deviatoric
    !! relaxation formulation:
    !! https://docs.frontistr.com/en/theory/constitutive/04_viscoelastic.html.
    !! Its one-branch form maps as \(G_\infty=G_0\mu_\infty\),
    !! \(G_1=G_0\mu_1\), \(\tau=\lambda_1\) and
    !! \(\boldsymbol q=2G_1\boldsymbol h^{(1)}\).
    !! The backward-Euler update and analytical derivatives above are
    !! independently derived from that continuous law. This revision
    !! reuses no third-party source code.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: nmandel = 6

    public :: maxwell_shear_t

    type, extends(material_model_t) :: maxwell_shear_t
        !! One-branch linear-viscoelastic Maxwell material in 3D.
        !!
        !! Construct with `maxwell_shear_t(bulk_modulus, equilibrium_shear_modulus, branch_shear_modulus, relaxation_time)`
        !! and initialize through `material_t%init` with `initial_input(6,npoints)`.
        !! Each `eval` requires `dt` in the chosen time unit and is only a trial until commit.
        real(wp) :: bulk_modulus = 0.0_wp
            !! Positive unrelaxed/relaxed bulk modulus \(K\), in the chosen stress unit.
        real(wp) :: equilibrium_shear_modulus = 0.0_wp
            !! Nonnegative long-time shear modulus \(G_\infty\), in the chosen stress unit.
        real(wp) :: branch_shear_modulus = 0.0_wp
            !! Nonnegative Maxwell branch shear modulus \(G_1\), in the chosen stress unit.
        real(wp) :: relaxation_time = 0.0_wp
            !! Positive Maxwell relaxation time \(\tau\), in the same time unit as `dt`.
    contains
        procedure :: describe => describe_maxwell_shear
            !! Declare material fields, parameters and supported operations.
        procedure :: eval_batch => evaluate_maxwell_shear
            !! Evaluate the batched material response, optional tangent and trial state.
        procedure :: initialize_sensitivity => initialize_maxwell_sensitivity
            !! Initialize exact derivatives of the material's reference state.
        procedure :: eval_sensitivity_batch => sensitivity_maxwell_shear
            !! Evaluate analytical parameter derivatives through the discrete material update.
    end type maxwell_shear_t

contains

    subroutine describe_maxwell_shear(self, spec, status)
        !! Declare 6-component Mandel input, response, real history and scalars.
        !!
        !! Nonfinite, invalid, or unrepresentable modulus combinations fail
        !! setup with `cf_invalid_argument`. This model requires `initial_input`,
        !! previous accepted strain and a strictly positive evaluation `dt`.
        class(maxwell_shear_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema populated with fields, physical metadata and model capabilities.
        integer, intent(out) :: status
            !! Operation status code. `cf_ok` indicates success.
        real(wp) :: total_shear

        status = cf_invalid_argument
        if (.not. ieee_is_finite(self%bulk_modulus)) return
        if (.not. ieee_is_finite(self%equilibrium_shear_modulus)) return
        if (.not. ieee_is_finite(self%branch_shear_modulus)) return
        if (.not. ieee_is_finite(self%relaxation_time)) return
        if (self%bulk_modulus <= 0.0_wp .or. self%relaxation_time <= 0.0_wp) return
        if (self%equilibrium_shear_modulus < 0.0_wp .or. self%branch_shear_modulus < 0.0_wp) return
        total_shear = self%equilibrium_shear_modulus + self%branch_shear_modulus
        if (.not. ieee_is_finite(total_shear)) return
        if (total_shear <= 0.0_wp) return
        if (.not. ieee_is_finite(2.0_wp*total_shear)) return
        if (.not. ieee_is_finite(self%bulk_modulus + (4.0_wp/3.0_wp)*total_shear)) return

        call spec%add_input("strain.small", nmandel, "mandel", "small_strain", "spatial", [3,3])
        call spec%add_output("stress.cauchy", nmandel, "mandel", "cauchy_stress", "spatial", [3,3])
        call spec%add_state("overstress.deviatoric", nmandel, "mandel", "internal_stress", "spatial", [3,3])
        call spec%add_parameter("bulk_modulus")
        call spec%add_parameter("equilibrium_shear_modulus")
        call spec%add_parameter("branch_shear_modulus")
        call spec%add_parameter("relaxation_time")
        call spec%set_requires_dt(.true.)
        call spec%set_requires_initial_input(.true.)
        call spec%set_requires_previous_input(.true.)
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_maxwell_shear

    subroutine initialize_maxwell_sensitivity(self, initial_input, initial_state, selected_parameters, sensitivity, status)
        !! Initially relaxed overstress is identically zero for every parameter.
        !!
        !! This exact initialization holds for any finite initial strain.
        !! Parametric derivatives of all six initial history components are zero.
        class(maxwell_shear_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        real(wp), intent(in) :: initial_input(:,:)
            !! Prescribed initial input used to initialize the model's reference state.
        real(wp), intent(in) :: initial_state(:,:)
            !! Initial real-valued history `(nstate,npoints)`.
        integer, intent(in) :: selected_parameters(:)
            !! One-based IDs of scalar parameters in the requested derivative-axis order.
        real(wp), intent(out) :: sensitivity(:,:,:)
            !! Exact initial real-history derivatives `(nstate,nselected,npoints)`.
        integer, intent(out) :: status
            !! Operation status code. `cf_ok` indicates success.

        sensitivity = 0.0_wp
        status = cf_ok
    end subroutine initialize_maxwell_sensitivity

    subroutine evaluate_maxwell_shear(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Trial stress and six-component overstress update. No allocations.
        !!
        !! The previous input/state refer to the last committed time level.
        !! `tangent(6,6,npoints)` is the exact backward-Euler derivative
        !! \(\partial\boldsymbol{\sigma}_{n+1}/\partial\boldsymbol{\varepsilon}_{n+1}\)
        !! at fixed committed history. No history is mutated on failure.
        class(maxwell_shear_t), intent(in) :: self
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

        call evaluate_points(self, x_old, x_new, state_old, state_trial, response, dt, point_status, tangent)
    end subroutine evaluate_maxwell_shear

    subroutine sensitivity_maxwell_shear(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Fused state and stress sensitivities including accepted history.
        !!
        !! `dsdp(6,nselected,npoints)` follows `init`'s parameter name order.
        !! Derivatives of `state_trial` are committed/reverted transactionally.
        class(maxwell_shear_t), intent(in) :: self
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

        call evaluate_points(self, x_old, x_new, state_old, state_trial, response, dt, point_status, tangent, selected_parameters, sensitivity_old, sensitivity_trial, dsdp)
    end subroutine sensitivity_maxwell_shear

    subroutine evaluate_points(self, x_old, x_new, state_old, state_trial, response, dt, point_status, tangent, selected_parameters, sensitivity_old, sensitivity_trial, dsdp)
        !! A preallocated batch updates all points independently, in Mandel order.
        class(maxwell_shear_t), intent(in) :: self
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
        integer, intent(inout) :: point_status(:)
            !! Per-point status codes. successful points must be assigned `cf_ok`.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact discrete input Jacobian `(noutput,ninput,npoints)`.
        integer, intent(in), optional :: selected_parameters(:)
            !! One-based IDs of scalar parameters in the requested derivative-axis order.
        real(wp), intent(in), optional :: sensitivity_old(:,:,:)
            !! Committed history derivatives `(nstate,nselected,npoints)`.
        real(wp), intent(out), optional :: sensitivity_trial(:,:,:)
            !! Updated trial history derivatives `(nstate,nselected,npoints)`.
        real(wp), intent(out), optional :: dsdp(:,:,:)
            !! Exact response derivatives `(noutput,nselected,npoints)` at fixed prescribed input histories.
        real(wp) :: deviation(nmandel)
        real(wp) :: increment(nmandel)
        real(wp) :: delta_deviation(nmandel)
        real(wp) :: q(nmandel)
        real(wp) :: trial_overstress(nmandel)
        real(wp) :: tangent_matrix(nmandel,nmandel)
        real(wp) :: time_ratio
        real(wp) :: attenuation
        real(wp) :: viscous_fraction
        real(wp) :: normal_time_fraction
        real(wp) :: time_derivative
        real(wp) :: effective_shear
        real(wp) :: trace_strain
        real(wp) :: dg_branch
        integer :: ipoint
        integer :: iparam
        integer :: i

        point_status = cf_invalid_argument
        if (dt <= self%relaxation_time) then
            time_ratio = dt/self%relaxation_time
            attenuation = 1.0_wp/(1.0_wp+time_ratio)
            viscous_fraction = time_ratio*attenuation
            normal_time_fraction = attenuation
        else
            time_ratio = self%relaxation_time/dt
            viscous_fraction = 1.0_wp/(1.0_wp+time_ratio)
            attenuation = time_ratio*viscous_fraction
            normal_time_fraction = viscous_fraction
        end if
        if (attenuation >= tiny(1.0_wp)) then
            effective_shear = self%equilibrium_shear_modulus + self%branch_shear_modulus*attenuation
        else
            ! A subnormal factor can lose relative precision before a large
            ! modulus restores scale. Keep the original time factors instead.
            effective_shear = self%equilibrium_shear_modulus + scaled_attenuation_product(self%branch_shear_modulus,self%relaxation_time,dt,viscous_fraction)
        end if
        if (.not. ieee_is_finite(effective_shear)) return

        if (present(tangent)) then
            call elastic_tangent(self%bulk_modulus, effective_shear, tangent_matrix)
            if (.not. all(ieee_is_finite(tangent_matrix))) return
        end if

        point_loop: do ipoint = 1, size(x_new,2)
            if (.not. all(ieee_is_finite(x_new(:,ipoint)))) cycle point_loop
            if (.not. all(ieee_is_finite(x_old(:,ipoint)))) cycle point_loop
            if (.not. all(ieee_is_finite(state_old(:,ipoint)))) cycle point_loop
            trace_strain = sum(x_new(1:3,ipoint))
            deviation = x_new(:,ipoint)
            ! Pairwise differences keep hydrostatic strains and increments
            ! exactly deviator-free, even when their trace rounds.
            deviation(1) = (x_new(1,ipoint)-x_new(2,ipoint))/3.0_wp+(x_new(1,ipoint)-x_new(3,ipoint))/3.0_wp
            deviation(2) = (x_new(2,ipoint)-x_new(1,ipoint))/3.0_wp+(x_new(2,ipoint)-x_new(3,ipoint))/3.0_wp
            deviation(3) = (x_new(3,ipoint)-x_new(1,ipoint))/3.0_wp+(x_new(3,ipoint)-x_new(2,ipoint))/3.0_wp
            increment = x_new(:,ipoint)-x_old(:,ipoint)
            delta_deviation = increment
            delta_deviation(1) = (increment(1)-increment(2))/3.0_wp+(increment(1)-increment(3))/3.0_wp
            delta_deviation(2) = (increment(2)-increment(1))/3.0_wp+(increment(2)-increment(3))/3.0_wp
            delta_deviation(3) = (increment(3)-increment(1))/3.0_wp+(increment(3)-increment(2))/3.0_wp
            if (.not. ieee_is_finite(trace_strain)) cycle point_loop
            if (.not. all(ieee_is_finite(deviation))) cycle point_loop
            if (.not. all(ieee_is_finite(delta_deviation))) cycle point_loop

            trial_overstress = state_old(:,ipoint) + 2.0_wp*self%branch_shear_modulus*delta_deviation
            if (all(ieee_is_finite(trial_overstress))) then
                if (attenuation >= tiny(1.0_wp)) then
                    q = attenuation*trial_overstress
                else
                    q = scaled_attenuation_product(trial_overstress,self%relaxation_time,dt,viscous_fraction)
                end if
            else
                ! Attenuate overflowing forcing from its original factors.
                ! Keep finite trial components on the usual update path.
                q = state_old(:,ipoint)
                where (ieee_is_finite(trial_overstress)) q = trial_overstress
                if (attenuation >= tiny(1.0_wp)) then
                    q = attenuation*q
                else
                    q = scaled_attenuation_product(q,self%relaxation_time,dt,viscous_fraction)
                end if
                do i = 1,nmandel
                    if (.not. ieee_is_finite(trial_overstress(i))) then
                        q(i) = q(i)+scaled_forcing_attenuation(self%branch_shear_modulus,delta_deviation(i),normal_time_fraction,dt,self%relaxation_time)
                    end if
                end do
            end if
            response(:,ipoint) = 2.0_wp*self%equilibrium_shear_modulus*deviation + q
            response(1:3,ipoint) = response(1:3,ipoint) + self%bulk_modulus*trace_strain
            if (.not. all(ieee_is_finite(q))) cycle point_loop
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop

            if (present(selected_parameters)) then
                parameter_loop: do iparam = 1, size(selected_parameters)
                    dg_branch = 0.0_wp
                    select case (selected_parameters(iparam))
                    case (1,2,4)
                        continue
                    case (3)
                        dg_branch = 1.0_wp
                    case default
                        exit parameter_loop
                    end select
                    if (attenuation >= tiny(1.0_wp)) then
                        sensitivity_trial(:,iparam,ipoint) = attenuation*(sensitivity_old(:,iparam,ipoint) + 2.0_wp*dg_branch*delta_deviation)
                    else
                        sensitivity_trial(:,iparam,ipoint) = scaled_attenuation_product(sensitivity_old(:,iparam,ipoint) + 2.0_wp*dg_branch*delta_deviation,self%relaxation_time,dt,viscous_fraction)
                    end if
                    if (selected_parameters(iparam) == 4) then
                        do i = 1,nmandel
                            if (.not. ieee_is_finite(trial_overstress(i)) .or. (self%branch_shear_modulus > 0.0_wp .and. abs(delta_deviation(i)) > 0.0_wp .and. abs(2.0_wp*self%branch_shear_modulus*delta_deviation(i)) < tiny(1.0_wp))) then
                                ! Preserve the original forcing factors when
                                ! rounding or overflow precedes time differentiation.
                                time_derivative = scaled_time_derivative(state_old(i,ipoint),normal_time_fraction,dt,self%relaxation_time) + scaled_forcing_time_derivative(self%branch_shear_modulus,delta_deviation(i),normal_time_fraction,dt,self%relaxation_time)
                            else
                                time_derivative = scaled_time_derivative(trial_overstress(i),normal_time_fraction,dt,self%relaxation_time)
                            end if
                            sensitivity_trial(i,iparam,ipoint) = sensitivity_trial(i,iparam,ipoint) + time_derivative
                        end do
                    end if
                    dsdp(:,iparam,ipoint) = sensitivity_trial(:,iparam,ipoint)
                    select case (selected_parameters(iparam))
                    case (1)
                        dsdp(1:3,iparam,ipoint) = dsdp(1:3,iparam,ipoint) + trace_strain
                    case (2)
                        dsdp(:,iparam,ipoint) = dsdp(:,iparam,ipoint) + 2.0_wp*deviation
                    case default
                        continue
                    end select
                end do parameter_loop
                if (iparam <= size(selected_parameters)) cycle point_loop
                if (.not. all(ieee_is_finite(sensitivity_trial(:,:,ipoint)))) then
                    cycle point_loop
                end if
                if (.not. all(ieee_is_finite(dsdp(:,:,ipoint)))) cycle point_loop
            end if

            state_trial(:,ipoint) = q
            if (present(tangent)) tangent(:,:,ipoint) = tangent_matrix
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_points

    pure elemental function scaled_attenuation_product(payload, relaxation_time, dt, viscous_fraction) result(value)
        !! Evaluate \(\mathrm{payload}\,\tau B/\Delta t\), with \(B=\Delta t/(\tau+\Delta t)\), before rounding a small attenuation.
        real(wp), intent(in) :: payload
            !! Modulus, trial overstress or incoming sensitivity to attenuate.
        real(wp), intent(in) :: relaxation_time
            !! Positive relaxation time \(\tau\).
        real(wp), intent(in) :: dt
            !! Positive time increment, greater than \(\tau\) in this branch.
        real(wp), intent(in) :: viscous_fraction
            !! Normally representable complementary fraction \(B\).
        real(wp) :: value

        value = scale((fraction(payload)*fraction(relaxation_time))*fraction(viscous_fraction)/fraction(dt), exponent(payload)+exponent(relaxation_time)+exponent(viscous_fraction)-exponent(dt))
    end function scaled_attenuation_product

    pure elemental function scaled_forcing_attenuation(branch_modulus, increment, time_fraction, dt, relaxation_time) result(value)
        !! Attenuate branch forcing without first forming an overflowing \(2G_1\Delta e\).
        real(wp), intent(in) :: branch_modulus
            !! Nonnegative branch shear modulus \(G_1\).
        real(wp), intent(in) :: increment
            !! One original deviatoric strain-increment component.
        real(wp), intent(in) :: time_fraction
            !! Normal-side \(A=\tau/(\tau+\Delta t)\) for \(\Delta t\leq\tau\) or \(B=\Delta t/(\tau+\Delta t)\) otherwise.
        real(wp), intent(in) :: dt
            !! Positive time increment.
        real(wp), intent(in) :: relaxation_time
            !! Positive relaxation time \(\tau\).
        real(wp) :: value

        if (dt <= relaxation_time) then
            value = scale(((2.0_wp*fraction(branch_modulus))*fraction(increment))*fraction(time_fraction), exponent(branch_modulus)+exponent(increment)+exponent(time_fraction))
        else
            value = scale(((2.0_wp*fraction(branch_modulus))*fraction(increment))*fraction(relaxation_time)*fraction(time_fraction)/fraction(dt), exponent(branch_modulus)+exponent(increment)+exponent(relaxation_time)+exponent(time_fraction)-exponent(dt))
        end if
    end function scaled_forcing_attenuation

    pure elemental function scaled_time_derivative(payload, time_fraction, dt, relaxation_time) result(value)
        !! Evaluate \(\mathrm{payload}\,\partial A/\partial\tau\) without quantizing a small time ratio or intermediate payload.
        real(wp), intent(in) :: payload
            !! Trial overstress before attenuation.
        real(wp), intent(in) :: time_fraction
            !! Normally representable \(A=\tau/(\tau+\Delta t)\) for \(\Delta t\leq\tau\) or \(B=\Delta t/(\tau+\Delta t)\) otherwise.
        real(wp), intent(in) :: dt
            !! Positive time increment.
        real(wp), intent(in) :: relaxation_time
            !! Positive relaxation time \(\tau\).
        real(wp) :: value

        if (dt <= relaxation_time) then
            value = scale((fraction(payload)*fraction(time_fraction)**2)*fraction(dt)/fraction(relaxation_time)**2, exponent(payload)+2*exponent(time_fraction)+exponent(dt)-2*exponent(relaxation_time))
        else
            value = scale(fraction(payload)*fraction(time_fraction)**2/fraction(dt), exponent(payload)+2*exponent(time_fraction)-exponent(dt))
        end if
    end function scaled_time_derivative

    pure elemental function scaled_forcing_time_derivative(branch_modulus, increment, time_fraction, dt, relaxation_time) result(value)
        !! Differentiate the branch forcing without first rounding \(2G_1\Delta e\).
        real(wp), intent(in) :: branch_modulus
            !! Nonnegative branch shear modulus \(G_1\).
        real(wp), intent(in) :: increment
            !! One original deviatoric strain-increment component.
        real(wp), intent(in) :: time_fraction
            !! Normal-side \(A=\tau/(\tau+\Delta t)\) for \(\Delta t\leq\tau\) or \(B=\Delta t/(\tau+\Delta t)\) otherwise.
        real(wp), intent(in) :: dt
            !! Positive time increment.
        real(wp), intent(in) :: relaxation_time
            !! Positive relaxation time \(\tau\).
        real(wp) :: value

        if (dt <= relaxation_time) then
            value = scale(((2.0_wp*fraction(branch_modulus))*fraction(increment))*fraction(time_fraction)**2*fraction(dt)/fraction(relaxation_time)**2, exponent(branch_modulus)+exponent(increment)+2*exponent(time_fraction)+exponent(dt)-2*exponent(relaxation_time))
        else
            value = scale(((2.0_wp*fraction(branch_modulus))*fraction(increment))*fraction(time_fraction)**2/fraction(dt), exponent(branch_modulus)+exponent(increment)+2*exponent(time_fraction)-exponent(dt))
        end if
    end function scaled_forcing_time_derivative

    pure subroutine elastic_tangent(bulk, shear, tangent)
        !! Isotropic 6 by 6 Mandel tangent, exact for backward Euler.
        real(wp), intent(in) :: bulk
            !! Bulk modulus \(K\) in the chosen stress unit.
        real(wp), intent(in) :: shear
            !! Shear modulus \(G\) in the chosen stress unit.
        real(wp), intent(out) :: tangent(nmandel,nmandel)
            !! Optional exact discrete input Jacobian `(noutput,ninput,npoints)`.
        real(wp) :: normal_coupling
        integer :: i

        normal_coupling = bulk - (2.0_wp/3.0_wp)*shear
        tangent = 0.0_wp
        tangent(1:3,1:3) = normal_coupling
        do i = 1, 3
            tangent(i,i) = bulk + (4.0_wp/3.0_wp)*shear
        end do
        do i = 4, nmandel
            tangent(i,i) = 2.0_wp*shear
        end do
    end subroutine elastic_tangent

end module constituf_maxwell_shear

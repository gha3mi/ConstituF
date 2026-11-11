! SPDX-License-Identifier: Apache-2.0

module constituf_maxwell_multi
    !! Linear three-dimensional generalized Maxwell material with \(N\) branches.
    !!
    !! With deviatoric small strain
    !! \(\boldsymbol{e}=\operatorname{dev}(\boldsymbol{\varepsilon})\),
    !! bulk modulus \(K>0\), equilibrium shear modulus \(G_\infty\geq0\),
    !! branch moduli \(G_j\geq0\) and relaxation times \(\tau_j>0\),
    !! its constitutive law is
    !! \[\boldsymbol{\sigma}
    !! =K\,\operatorname{tr}(\boldsymbol\varepsilon)\boldsymbol I+
    !! 2G_\infty\boldsymbol e+\sum_{j=1}^N\boldsymbol q_j,\qquad
    !! \dot{\boldsymbol q}_j+\boldsymbol q_j/\tau_j
    !! =2G_j\dot{\boldsymbol e}.\]
    !! A material with \(G_\infty+\sum_{j=1}^N G_j>0\) is admissible. \(G_\infty=0\) gives a
    !! generalized Maxwell shear fluid with an elastic volumetric response.
    !! Moduli, branch overstresses and Cauchy stress share the chosen
    !! stress unit. Relaxation times and `dt` share one time unit. Strain is dimensionless.
    !! Each symmetric 3D spatial tensor \(\boldsymbol{A}\) is packed as
    !! \[\operatorname{mandel}(\boldsymbol{A})=
    !! [A_{11},A_{22},A_{33},\sqrt{2}A_{23},
    !! \sqrt{2}A_{13},\sqrt{2}A_{12}]^{\mathsf T}.\]
    !! Strains and stress are work conjugate.
    !!
    !! Given the last *committed* input \(\boldsymbol{\varepsilon}_n\), each backward-Euler
    !! trial with \(\Delta t>0\) evaluates independently, for every branch:
    !! \[D_j=1+\Delta t/\tau_j,\quad
    !! \boldsymbol q_{j,n+1}=
    !! [\boldsymbol q_{j,n}+2G_j(\boldsymbol e_{n+1}-
    !! \boldsymbol e_n)]/D_j.\]
    !! The exact consistent tangent is
    !! \[\boldsymbol C_\mathrm{alg}=
    !! K\boldsymbol I\otimes\boldsymbol I+
    !! 2(G_\infty+\sum_j G_j/D_j)\boldsymbol I_\mathrm{dev}.\]
    !! Initial branch stress and its parameter derivatives are zero:
    !! input at `init` defines an initially relaxed reference state.
    !!
    !! All parameter derivatives below propagate accepted branch histories.
    !! Index-based names `branch_shear_modulus_N` and
    !! `relaxation_time_N` (one-based) are resolved at `init` alongside
    !! `bulk_modulus` and `equilibrium_shear_modulus`. For selected \(p\),
    !! \[\frac{\mathrm d\boldsymbol q_{j,n+1}}{\mathrm dp}
    !! =[\frac{\mathrm d\boldsymbol q_{j,n}}{\mathrm dp}
    !! +2\delta_{p,G_j}\Delta\boldsymbol e-
    !! \boldsymbol q_{j,n+1}\frac{\mathrm dD_j}{\mathrm dp}]/D_j,\quad
    !! \frac{\partial D_j}{\partial \tau_j}
    !! =-\Delta t/\tau_j^2.\]
    !! Parameter sensitivities hold all prescribed strain inputs and time
    !! increments fixed. No structural/global sensitivity is claimed.
    !!
    !! State field `overstress.branches` is packed as a `(6,N)` array for
    !! every point, with six Mandel components per branch and no scalar
    !! prestress. Strictly linear, isothermal, small-strain and deviatoric
    !! relaxation. No nonlinear, volumetric, or finite-strain effects.
    !!
    !! Continuous generalized-Maxwell reference: FrontISTR, Theory Manual,
    !! "Viscoelastic Material", section "Generalized Maxwell Model":
    !! https://docs.frontistr.com/en/theory/constitutive/04_viscoelastic.html .
    !! Its \(G_0=G_\infty+\sum_{j=1}^N G_j\),
    !! \(\mu_\infty=G_\infty/G_0\), \(\mu_j=G_j/G_0\) and
    !! \(\lambda_j=\tau_j\) give this law with
    !! \(\boldsymbol{q}_j=2G_j\boldsymbol{h}_j\). Backward Euler,
    !! as stated above, is the discrete update used here.
    !! Background: R. M. Christensen, Theory of Viscoelasticity: An
    !! Introduction, 2nd ed., Academic Press (1982),
    !! ISBN 978-0-12-174252-2, DOI 10.1016/B978-0-12-174252-2.X5001-7.
    !! Code provenance: this backward-Euler update and its analytical
    !! derivatives are independently derived from the stated continuous
    !! law. No third-party software code is reused in this revision.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: nmandel = 6
    public :: maxwell_multi_t

    type, extends(material_model_t) :: maxwell_multi_t
        !! Construct with `maxwell_multi_t(K,G_inf,[G_1,...],[tau_1,...])`.
        !!
        !! Branch arrays must be allocated, nonempty and equally sized,
        !! with nonnegative moduli and strictly positive relaxation times.
        !! Elements are paired in array order independently of lower bounds.
        !! Parameter suffixes always use one-based branch positions.
        !! The branch arrays are cloned by `material_t%init`. Callers cannot
        !! modify live model parameters afterward. `init` requires `initial_input`.
        real(wp) :: bulk_modulus = 0.0_wp
            !! \(K\) in the chosen stress unit, strictly positive.
        real(wp) :: equilibrium_shear_modulus = 0.0_wp
            !! \(G_\infty\) in the chosen stress unit, nonnegative.
        real(wp), allocatable :: branch_shear_moduli(:)
            !! \(G_j\) in the chosen stress unit, nonnegative.
        real(wp), allocatable :: relaxation_times(:)
            !! \(\tau_j\) in the same time unit as `dt`, strictly positive.
    contains
        procedure :: describe => describe_maxwell_multi
            !! Declare material fields, parameters and supported operations.
        procedure :: eval_batch => evaluate_maxwell_multi
            !! Evaluate the batched material response, optional tangent and trial state.
        procedure :: initialize_sensitivity => initialize_multi_sensitivity
            !! Initialize exact derivatives of the material's reference state.
        procedure :: eval_sensitivity_batch => sensitivity_maxwell_multi
            !! Evaluate analytical parameter derivatives through the discrete material update.
    end type maxwell_multi_t

contains

    subroutine describe_maxwell_multi(self, spec, status)
        !! Validate material parameters and register `strain.small` input and `stress.cauchy` response fields.
        class(maxwell_multi_t), intent(in) :: self
            !! Material model instance containing its fixed constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema populated with fields, physical metadata and model capabilities.
        integer, intent(out) :: status
            !! Operation status code. `cf_ok` indicates success.
        real(wp) :: total_shear
        integer :: j
        integer :: n
        character(48) :: name

        status = cf_invalid_argument
        if (.not. allocated(self%branch_shear_moduli)) return
        if (.not. allocated(self%relaxation_times)) return
        n = size(self%branch_shear_moduli)
        if (n < 1 .or. n > huge(n)/nmandel) return
        if (n /= size(self%relaxation_times)) return
        if (.not. ieee_is_finite(self%bulk_modulus) .or. .not. ieee_is_finite(self%equilibrium_shear_modulus)) return
        if (.not. all(ieee_is_finite(self%branch_shear_moduli))) return
        if (.not. all(ieee_is_finite(self%relaxation_times))) return
        if (self%bulk_modulus <= 0.0_wp .or. self%equilibrium_shear_modulus < 0.0_wp) return
        if (any(self%branch_shear_moduli < 0.0_wp) .or. any(self%relaxation_times <= 0.0_wp)) return
        total_shear = self%equilibrium_shear_modulus + sum(self%branch_shear_moduli)
        if (.not. ieee_is_finite(total_shear)) return
        if (total_shear <= 0.0_wp) return
        if (.not. ieee_is_finite(2.0_wp*total_shear)) return
        if (.not. ieee_is_finite(self%bulk_modulus + (4.0_wp/3.0_wp)*total_shear)) return

        call spec%add_input("strain.small", nmandel, "mandel", "small_strain", "spatial", [3,3])
        call spec%add_output("stress.cauchy", nmandel, "mandel", "cauchy_stress", "spatial", [3,3])
        call spec%add_state("overstress.branches", nmandel*n, "custom.mandel_branches", "internal_stress", "spatial", [nmandel,n])
        call spec%add_parameter("bulk_modulus")
        call spec%add_parameter("equilibrium_shear_modulus")
        do j = 1, n
            write(name, "(a,i0)") "branch_shear_modulus_", j
            call spec%add_parameter(trim(name))
            write(name, "(a,i0)") "relaxation_time_", j
            call spec%add_parameter(trim(name))
        end do
        call spec%set_requires_dt(.true.)
        call spec%set_requires_previous_input(.true.)
        call spec%set_requires_initial_input(.true.)
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_maxwell_multi

    subroutine initialize_multi_sensitivity(self, initial_input, initial_state, selected_parameters, sensitivity, status)
        !! Initialize analytical parameter derivatives of the reference real history.
        class(maxwell_multi_t), intent(in) :: self
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
    end subroutine initialize_multi_sensitivity

    subroutine evaluate_maxwell_multi(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Evaluate `stress.cauchy` for the current prescribed `strain.small` fields and update trial history.
        class(maxwell_multi_t), intent(in) :: self
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
    end subroutine evaluate_maxwell_multi

    subroutine sensitivity_maxwell_multi(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Differentiate the discrete `stress.cauchy` update with respect to selected scalar material parameters.
        class(maxwell_multi_t), intent(in) :: self
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
    end subroutine sensitivity_maxwell_multi

    subroutine evaluate_points(self, x_old, x_new, state_old, state_trial, response, dt, point_status, tangent, selected_parameters, sensitivity_old, sensitivity_trial, dsdp)
        !! Evaluate `stress.cauchy` with optional exact Jacobians and propagated history derivatives at each point.
        class(maxwell_multi_t), intent(in) :: self
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
        real(wp) :: attenuation(size(self%branch_shear_moduli))
        real(wp) :: viscous_fraction(size(self%branch_shear_moduli))
        real(wp) :: tangent_matrix(nmandel,nmandel)
        real(wp) :: dev(nmandel)
        real(wp) :: increment(nmandel)
        real(wp) :: delta(nmandel)
        real(wp) :: q(nmandel)
        real(wp) :: trial_overstress(nmandel)
        real(wp) :: derivative(nmandel)
        real(wp) :: shear
        real(wp) :: lambda
        real(wp) :: tr
        real(wp) :: time_ratio
        real(wp) :: branch_modulus
        real(wp) :: relaxation_time
        real(wp) :: normal_time_fraction
        real(wp) :: time_derivative
        integer :: ipoint
        integer :: i
        integer :: j
        integer :: k
        integer :: first
        integer :: last
        integer :: id
        integer :: shear_lower
        integer :: time_lower

        point_status = cf_invalid_argument
        shear_lower = lbound(self%branch_shear_moduli,1)
        time_lower = lbound(self%relaxation_times,1)
        do j = 1, size(attenuation)
            if (dt <= self%relaxation_times(time_lower+(j-1))) then
                time_ratio = dt/self%relaxation_times(time_lower+(j-1))
                attenuation(j) = 1.0_wp/(1.0_wp+time_ratio)
                viscous_fraction(j) = time_ratio*attenuation(j)
            else
                time_ratio = self%relaxation_times(time_lower+(j-1))/dt
                viscous_fraction(j) = 1.0_wp/(1.0_wp+time_ratio)
                attenuation(j) = time_ratio*viscous_fraction(j)
            end if
        end do
        if (all(attenuation >= tiny(1.0_wp))) then
            shear = self%equilibrium_shear_modulus + sum(self%branch_shear_moduli*attenuation)
        else
            ! Preserve the original time factors before a subnormal attenuation
            ! loses relative precision. Each branch contributes exactly once.
            shear = 0.0_wp
            do j = 1, size(attenuation)
                branch_modulus = self%branch_shear_moduli(shear_lower+(j-1))
                if (attenuation(j) >= tiny(1.0_wp)) then
                    shear = shear + branch_modulus*attenuation(j)
                else
                    relaxation_time = self%relaxation_times(time_lower+(j-1))
                    shear = shear + scaled_attenuation_product(branch_modulus,relaxation_time,dt,viscous_fraction(j))
                end if
            end do
            shear = self%equilibrium_shear_modulus + shear
        end if
        if (.not. ieee_is_finite(shear)) return
        if (present(tangent)) then
            lambda = self%bulk_modulus - (2.0_wp/3.0_wp)*shear
            tangent_matrix = 0.0_wp
            tangent_matrix(1:3,1:3) = lambda
            do j = 1, 3
                tangent_matrix(j,j) = self%bulk_modulus + (4.0_wp/3.0_wp)*shear
            end do
            do j = 4, nmandel
                tangent_matrix(j,j) = 2.0_wp*shear
            end do
            if (.not. all(ieee_is_finite(tangent_matrix))) return
        end if
        point_loop: do ipoint = 1, size(x_new,2)
            if (.not. all(ieee_is_finite(x_new(:,ipoint)))) cycle point_loop
            if (.not. all(ieee_is_finite(x_old(:,ipoint)))) cycle point_loop
            if (.not. all(ieee_is_finite(state_old(:,ipoint)))) cycle point_loop
            if (present(selected_parameters)) then
                if (.not. all(ieee_is_finite(sensitivity_old(:,:,ipoint)))) cycle point_loop
            end if
            tr = sum(x_new(1:3,ipoint))
            dev = x_new(:,ipoint)
            ! Preserve exact zero deviators for hydrostatic strains and increments.
            dev(1) = (x_new(1,ipoint)-x_new(2,ipoint))/3.0_wp+(x_new(1,ipoint)-x_new(3,ipoint))/3.0_wp
            dev(2) = (x_new(2,ipoint)-x_new(1,ipoint))/3.0_wp+(x_new(2,ipoint)-x_new(3,ipoint))/3.0_wp
            dev(3) = (x_new(3,ipoint)-x_new(1,ipoint))/3.0_wp+(x_new(3,ipoint)-x_new(2,ipoint))/3.0_wp
            increment = x_new(:,ipoint)-x_old(:,ipoint)
            delta = increment
            delta(1) = (increment(1)-increment(2))/3.0_wp+(increment(1)-increment(3))/3.0_wp
            delta(2) = (increment(2)-increment(1))/3.0_wp+(increment(2)-increment(3))/3.0_wp
            delta(3) = (increment(3)-increment(1))/3.0_wp+(increment(3)-increment(2))/3.0_wp
            if (.not. all(ieee_is_finite(dev)) .or. .not. all(ieee_is_finite(delta))) cycle point_loop
            response(:,ipoint) = 2.0_wp*self%equilibrium_shear_modulus*dev
            response(1:3,ipoint) = response(1:3,ipoint) + self%bulk_modulus*tr
            if (present(selected_parameters)) then
                sensitivity_trial(:,:,ipoint) = sensitivity_old(:,:,ipoint)
                dsdp(:,:,ipoint) = 0.0_wp
                do k = 1, size(selected_parameters)
                    id = selected_parameters(k)
                    if (id == 1) dsdp(1:3,k,ipoint) = tr
                    if (id == 2) dsdp(:,k,ipoint) = 2.0_wp*dev
                end do
            end if
            do j = 1, size(attenuation)
                first = nmandel*(j-1)+1
                last = first+nmandel-1
                branch_modulus = self%branch_shear_moduli(shear_lower+(j-1))
                relaxation_time = self%relaxation_times(time_lower+(j-1))
                if (dt <= relaxation_time) then
                    normal_time_fraction = attenuation(j)
                else
                    normal_time_fraction = viscous_fraction(j)
                end if
                trial_overstress = state_old(first:last,ipoint) + 2.0_wp*branch_modulus*delta
                if (all(ieee_is_finite(trial_overstress))) then
                    if (attenuation(j) >= tiny(1.0_wp)) then
                        q = attenuation(j)*trial_overstress
                    else
                        q = scaled_attenuation_product(trial_overstress,relaxation_time,dt,viscous_fraction(j))
                    end if
                else
                    ! Attenuate overflowing forcing from its original factors.
                    ! Keep finite trial components on the usual update path.
                    q = state_old(first:last,ipoint)
                    where (ieee_is_finite(trial_overstress)) q = trial_overstress
                    if (attenuation(j) >= tiny(1.0_wp)) then
                        q = attenuation(j)*q
                    else
                        q = scaled_attenuation_product(q,relaxation_time,dt,viscous_fraction(j))
                    end if
                    do i = 1,nmandel
                        if (.not. ieee_is_finite(trial_overstress(i))) then
                            q(i) = q(i)+scaled_forcing_attenuation(branch_modulus,delta(i),normal_time_fraction,dt,relaxation_time)
                        end if
                    end do
                end if
                state_trial(first:last,ipoint) = q
                response(:,ipoint) = response(:,ipoint)+q
                if (present(selected_parameters)) then
                    do k = 1, size(selected_parameters)
                        id = selected_parameters(k)
                        derivative = sensitivity_old(first:last,k,ipoint)
                        if (id == 2*j+1) derivative = derivative + 2.0_wp*delta
                        if (attenuation(j) >= tiny(1.0_wp)) then
                            derivative = attenuation(j)*derivative
                        else
                            derivative = scaled_attenuation_product(derivative,relaxation_time,dt,viscous_fraction(j))
                        end if
                        if (id == 2*j+2) then
                            do i = 1,nmandel
                                if (.not. ieee_is_finite(trial_overstress(i)) .or. (branch_modulus > 0.0_wp .and. abs(delta(i)) > 0.0_wp .and. abs(2.0_wp*branch_modulus*delta(i)) < tiny(1.0_wp))) then
                                    ! Preserve the original forcing factors when
                                    ! rounding or overflow precedes time differentiation.
                                    time_derivative = scaled_time_derivative(state_old(first+i-1,ipoint),normal_time_fraction,dt,relaxation_time) + scaled_forcing_time_derivative(branch_modulus,delta(i),normal_time_fraction,dt,relaxation_time)
                                else
                                    time_derivative = scaled_time_derivative(trial_overstress(i),normal_time_fraction,dt,relaxation_time)
                                end if
                                derivative(i) = derivative(i) + time_derivative
                            end do
                        end if
                        sensitivity_trial(first:last,k,ipoint) = derivative
                        dsdp(:,k,ipoint) = dsdp(:,k,ipoint) + derivative
                    end do
                end if
            end do
            if (.not. all(ieee_is_finite(state_trial(:,ipoint)))) cycle point_loop
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop
            if (present(selected_parameters)) then
                if (.not. all(ieee_is_finite(sensitivity_trial(:,:,ipoint)))) cycle point_loop
                if (.not. all(ieee_is_finite(dsdp(:,:,ipoint)))) cycle point_loop
            end if
            if (present(tangent)) tangent(:,:,ipoint) = tangent_matrix
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_points

    pure elemental function scaled_attenuation_product(payload, relaxation_time, dt, viscous_fraction) result(value)
        !! Evaluate \(\mathrm{payload}\,\tau B/\Delta t\), with \(B=\Delta t/(\tau+\Delta t)\), before rounding a small attenuation.
        real(wp), intent(in) :: payload
            !! Modulus, trial overstress or incoming sensitivity to attenuate.
        real(wp), intent(in) :: relaxation_time
            !! Positive branch relaxation time \(\tau\).
        real(wp), intent(in) :: dt
            !! Positive time increment, greater than \(\tau\) in this branch.
        real(wp), intent(in) :: viscous_fraction
            !! Normally representable complementary fraction \(B\).
        real(wp) :: value

        value = scale((fraction(payload)*fraction(relaxation_time))*fraction(viscous_fraction)/fraction(dt), exponent(payload)+exponent(relaxation_time)+exponent(viscous_fraction)-exponent(dt))
    end function scaled_attenuation_product

    pure elemental function scaled_forcing_attenuation(branch_modulus, increment, time_fraction, dt, relaxation_time) result(value)
        !! Attenuate branch forcing without first forming an overflowing \(2G_j\Delta e\).
        real(wp), intent(in) :: branch_modulus
            !! Nonnegative branch shear modulus \(G_j\).
        real(wp), intent(in) :: increment
            !! One original deviatoric strain-increment component.
        real(wp), intent(in) :: time_fraction
            !! Normal-side \(A=\tau/(\tau+\Delta t)\) for \(\Delta t\leq\tau\) or \(B=\Delta t/(\tau+\Delta t)\) otherwise.
        real(wp), intent(in) :: dt
            !! Positive time increment.
        real(wp), intent(in) :: relaxation_time
            !! Positive branch relaxation time \(\tau\).
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
            !! Positive branch relaxation time \(\tau\).
        real(wp) :: value

        if (dt <= relaxation_time) then
            value = scale((fraction(payload)*fraction(time_fraction)**2)*fraction(dt)/fraction(relaxation_time)**2, exponent(payload)+2*exponent(time_fraction)+exponent(dt)-2*exponent(relaxation_time))
        else
            value = scale(fraction(payload)*fraction(time_fraction)**2/fraction(dt), exponent(payload)+2*exponent(time_fraction)-exponent(dt))
        end if
    end function scaled_time_derivative

    pure elemental function scaled_forcing_time_derivative(branch_modulus, increment, time_fraction, dt, relaxation_time) result(value)
        !! Differentiate the branch forcing without first rounding \(2G_j\Delta e\).
        real(wp), intent(in) :: branch_modulus
            !! Nonnegative branch shear modulus \(G_j\).
        real(wp), intent(in) :: increment
            !! One original deviatoric strain-increment component.
        real(wp), intent(in) :: time_fraction
            !! Normal-side \(A=\tau/(\tau+\Delta t)\) for \(\Delta t\leq\tau\) or \(B=\Delta t/(\tau+\Delta t)\) otherwise.
        real(wp), intent(in) :: dt
            !! Positive time increment.
        real(wp), intent(in) :: relaxation_time
            !! Positive branch relaxation time \(\tau\).
        real(wp) :: value

        if (dt <= relaxation_time) then
            value = scale(((2.0_wp*fraction(branch_modulus))*fraction(increment))*fraction(time_fraction)**2*fraction(dt)/fraction(relaxation_time)**2, exponent(branch_modulus)+exponent(increment)+2*exponent(time_fraction)+exponent(dt)-2*exponent(relaxation_time))
        else
            value = scale(((2.0_wp*fraction(branch_modulus))*fraction(increment))*fraction(time_fraction)**2/fraction(dt), exponent(branch_modulus)+exponent(increment)+2*exponent(time_fraction)-exponent(dt))
        end if
    end function scaled_forcing_time_derivative

end module constituf_maxwell_multi

! SPDX-License-Identifier: Apache-2.0

module constituf_cohesive_mode_i
    !! Rate-independent, tension-only, bilinear Mode I cohesive interface.
    !!
    !! The normal opening \(\delta\geq0\) is a prescribed scalar displacement
    !! jump, not a strain. A scalar history \(\kappa_n\) stores the largest
    !! previously accepted opening. The exact discrete update is
    !! \[\kappa_{n+1}=\max(\kappa_n,\delta_{n+1}).\]
    !! For penalty stiffness \(K>0\), damage-onset opening \(0<\delta_0\)
    !! and complete-separation opening \(\delta_f>\delta_0\), define
    !! \[r(\kappa)=1-d(\kappa)=
    !! \begin{cases}
    !! 1,&0\leq\kappa\leq\delta_0,\\[2pt]
    !! \displaystyle\frac{\delta_0(\delta_f-\kappa)}
    !! {\kappa(\delta_f-\delta_0)},&\delta_0<\kappa<\delta_f,\\[2pt]
    !! 0,&\kappa\geq\delta_f.
    !! \end{cases}\]
    !! Normal traction is \(t_{n+1}=K\,r(\kappa_{n+1})\delta_{n+1}\).
    !! On monotone loading: elastic \(K\delta\) up to the peak
    !! \(t_{\max}=K\delta_0\), then linear softening to zero at
    !! \(\delta_f\). Unloading/reloading at fixed \(\kappa_n\) is linear
    !! with degraded stiffness \(Kr(\kappa_n)\), without healing.
    !! The mode-I fracture energy per unit area for complete monotone
    !! opening is \(G_{\mathrm{Ic}}=K\delta_0\delta_f/2\) (energy per area).
    !!
    !! The consistent local derivative at fixed committed history is
    !! \[\frac{\partial t_{n+1}}{\partial\delta_{n+1}}=
    !! \begin{cases}
    !! K,&\kappa_{n+1}<\delta_0,\\[2pt]
    !! -K\delta_0/(\delta_f-\delta_0),
    !! &\delta_{n+1}>\kappa_n,\ \delta_0<\delta_{n+1}<\delta_f,\\[2pt]
    !! Kr(\kappa_n),&\delta_{n+1}<\kappa_n,\quad
    !! \delta_0<\kappa_n<\delta_f,\\[2pt]
    !! 0,&\kappa_n\geq\delta_f\ \text{or}\ \kappa_{n+1}>\delta_f.
    !! \end{cases}\]
    !! At damage onset, at first arrival at complete separation, and at
    !! the loading/unloading switch in the softening regime, the classical
    !! opening tangent is not unique. Requests including a tangent at these
    !! points return `cf_nondifferentiable` and must not be committed.
    !! Traction-only evaluations remain defined at the switch. Once complete
    !! separation is committed (\(\kappa_n\geq\delta_f\)), traction and its
    !! opening tangent remain zero for every admissible opening.
    !!
    !! For the selected parameters \(p=(K,\delta_0,\delta_f)\) and fixed
    !! prescribed opening history, \(\mathrm{d}\kappa_n/\mathrm{d}p=0\).
    !! For \(\delta_0<\kappa_{n+1}<\delta_f\), local analytical derivatives:
    !! \[\frac{\mathrm{d}t}{\mathrm{d}K}=\delta r,\qquad
    !! \frac{\mathrm{d}t}{\mathrm{d}\delta_0}=
    !! K\delta\frac{\delta_f(\delta_f-\kappa)}
    !! {\kappa(\delta_f-\delta_0)^2},\qquad
    !! \frac{\mathrm{d}t}{\mathrm{d}\delta_f}=
    !! K\delta\frac{\delta_0(\kappa-\delta_0)}
    !! {\kappa(\delta_f-\delta_0)^2}.\]
    !! In the intact regime only \(\mathrm{d}t/\mathrm{d}K=\delta\)
    !! is nonzero. After full failure all derivatives are zero.
    !! At positive opening, threshold parameter derivatives are undefined
    !! when the maximum history equals damage onset/completion, including
    !! during unloading. Requesting the corresponding onset or failure
    !! parameter returns `cf_nondifferentiable`. The other parameter
    !! derivatives remain defined. At zero opening all traction parameter
    !! derivatives vanish, including at a threshold history.
    !!
    !! Input `opening.normal` and committed history `opening.maximum`
    !! share the chosen length unit. Output `traction.normal` uses the
    !! chosen stress unit. Normal opening
    !! is nonnegative. **Compression, contact, friction and shear separation
    !! are not modeled:** negative opening is rejected, and this is NOT a
    !! mixed-mode cohesive law. No artificial viscosity or `dt` is required.
    !! Cohesive softening may make a global finite-element tangent indefinite.
    !! Constitutive derivatives are local and omit global solution effects.
    !!
    !! Adapted Mode I specialization of: P. P. Camanho and C. G. Dávila,
    !! *Mixed-Mode Decohesion Finite Elements for the Simulation of
    !! Delamination in Composite Materials*, NASA/TM-2002-211737, 2002.
    !! Section 2.2.1, equations (21), (24) and (25):
    !! https://ntrs.nasa.gov/api/citations/20020053651/downloads/20020053651.pdf .
    !! Its Mode I properties \(N\) and \(G_{\mathrm{Ic}}\) are
    !! reparameterized here as \(N=K\delta_0\) and
    !! \(G_{\mathrm{Ic}}=K\delta_0\delta_f/2\). The derivatives above
    !! hold the other independent parameters \(K\), \(\delta_0\) and \(\delta_f\) fixed.
    !! The published model also supports mixed modes and compression.
    !! Those mechanisms are deliberately excluded here.
    !! Code provenance: the scalar maximum-history update and its
    !! analytical derivatives are independently derived from this stated
    !! specialization. No third-party software code is reused in this revision.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument, cf_nondifferentiable
    implicit none (type, external)
    private

    public :: cohesive_mode_i_t

    type, extends(material_model_t) :: cohesive_mode_i_t
        !! Bilinear normal-separation law with irreversible maximum opening.
        !!
        !! Construct as `cohesive_mode_i_t(stiffness,onset_opening,failure_opening)`.
        !! Stiffness has dimensions of stress per length. Both openings
        !! share the chosen length unit. A `material_t` owns trial and
        !! committed maximum opening per point. Initial history is zero.
        real(wp) :: stiffness = 0.0_wp
            !! Elastic interface stiffness \(K>0\), stress per length.
        real(wp) :: onset_opening = 0.0_wp
            !! Opening \(\delta_0>0\) at the maximum tensile traction, in the chosen length unit.
        real(wp) :: failure_opening = 0.0_wp
            !! Opening \(\delta_f>\delta_0\) at complete separation, in the chosen length unit.
    contains
        procedure :: describe => describe_cohesive
            !! Declare the constitutive fields, parameter metadata and supported capabilities.
        procedure :: eval_batch => evaluate_cohesive
            !! Evaluate the batched constitutive response, optional tangent and trial history.
        procedure :: initialize_sensitivity => initialize_cohesive_sensitivity
            !! Initialize analytical sensitivities of the reference material history.
        procedure :: eval_sensitivity_batch => sensitivity_cohesive
            !! Evaluate exact parameter-path derivatives and corresponding trial history.
    end type cohesive_mode_i_t

contains

    subroutine describe_cohesive(self, spec, status)
        !! Validate bilinear parameters and register scalar interface fields.
        class(cohesive_mode_i_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema receiving field layouts, parameters and supported operations.
        integer, intent(out) :: status
            !! Operation status. `cf_ok` indicates a valid model or successful initialization.
        real(wp) :: slope
        real(wp) :: fracture_energy

        status = cf_invalid_argument
        if (.not. all(ieee_is_finite([self%stiffness,self%onset_opening, self%failure_opening]))) return
        if (self%stiffness <= 0.0_wp .or. self%onset_opening <= 0.0_wp) return
        if (self%failure_opening <= self%onset_opening) return
        slope = self%stiffness*self%onset_opening/ (self%failure_opening-self%onset_opening)
        fracture_energy = (self%stiffness*self%onset_opening)* (self%failure_opening/2.0_wp)
        if (.not. ieee_is_finite(slope)) return
        if (.not. ieee_is_finite(fracture_energy)) return
        if (slope <= 0.0_wp .or. fracture_energy <= 0.0_wp) return

        call spec%add_input("opening.normal", 1, "scalar", "normal_separation", "interface")
        call spec%add_output("traction.normal", 1, "scalar", "normal_traction", "interface")
        call spec%add_state("opening.maximum", 1, "scalar", "maximum_opening", "interface")
        call spec%add_parameter("stiffness")
        call spec%add_parameter("onset_opening")
        call spec%add_parameter("failure_opening")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_cohesive

    subroutine initialize_cohesive_sensitivity(self, initial_input, initial_state, selected_parameters, sensitivity, status)
        !! Virgin maximum opening is zero and parameter independent.
        class(cohesive_mode_i_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        real(wp), intent(in) :: initial_input(:,:)
            !! Initial packed material inputs used to establish the reference history.
        real(wp), intent(in) :: initial_state(:,:)
            !! Initial real-valued history `(nstate,npoints)`.
        integer, intent(in) :: selected_parameters(:)
            !! One-based IDs of selected model parameters, in sensitivity-axis order.
        real(wp), intent(out) :: sensitivity(:,:,:)
            !! Reference-state derivatives `(nstate,nselected,npoints)` initialized consistently with the model.
        integer, intent(out) :: status
            !! Operation status. `cf_ok` indicates a valid model or successful initialization.

        sensitivity = 0.0_wp
        status = cf_ok
    end subroutine initialize_cohesive_sensitivity

    subroutine evaluate_cohesive(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Produce a provisional traction and maximum-opening history.
        !!
        !! `x_new(1,npoints)`, `response(1,npoints)`, `state_trial(1,npoints)`.
        !! An invalid or nonsmooth derivative trial cannot be committed.
        class(cohesive_mode_i_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed input `(ninput,npoints)` if required, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `opening.normal` input `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued state `(nstate,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Trial real-valued state `(nstate,npoints)` returned by this evaluation.
        real(wp), intent(out) :: response(:,:)
            !! Packed `traction.normal` response `(noutput,npoints)` for the evaluated material points.
        real(wp), intent(in) :: dt
            !! Prescribed time increment. It is positive when the material requires time integration.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact discrete response Jacobian `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(inout) :: point_status(:)
            !! Per-point status values; successful evaluations must set each entry to `cf_ok`.

        call evaluate_points(self,x_new,state_old,state_trial,response, point_status,tangent)
    end subroutine evaluate_cohesive

    subroutine sensitivity_cohesive(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Exact sensitivity to stiffness/onset/failure at fixed opening history.
        !!
        !! `dsdp(1,nselected,npoints)` follows the order selected at `init`.
        !! The derivative of max-opening history is identically zero.
        class(cohesive_mode_i_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed input `(ninput,npoints)` if required, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `opening.normal` input `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued state `(nstate,npoints)`.
        integer, intent(in) :: selected_parameters(:)
            !! One-based IDs of selected model parameters, in sensitivity-axis order.
        real(wp), intent(in) :: sensitivity_old(:,:,:)
            !! Committed real-state derivatives `(nstate,nselected,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Trial real-valued state `(nstate,npoints)` returned by this evaluation.
        real(wp), intent(out) :: sensitivity_trial(:,:,:)
            !! Trial real-state derivatives `(nstate,nselected,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `traction.normal` response `(noutput,npoints)` for the evaluated material points.
        real(wp), intent(out) :: dsdp(:,:,:)
            !! Exact total response derivatives `(noutput,nselected,npoints)` at fixed prescribed input history.
        real(wp), intent(in) :: dt
            !! Prescribed time increment. It is positive when the material requires time integration.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact discrete response Jacobian `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(inout) :: point_status(:)
            !! Per-point status values; successful evaluations must set each entry to `cf_ok`.

        sensitivity_trial = sensitivity_old
        call evaluate_points(self,x_new,state_old,state_trial,response, point_status,tangent,selected_parameters,dsdp)
    end subroutine sensitivity_cohesive

    subroutine evaluate_points(self,x_new,state_old,state_trial,response, point_status,tangent,selected_parameters,dsdp)
        !! No dynamic allocation, no name lookup and no committed-state mutation.
        class(cohesive_mode_i_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `opening.normal` input `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued state `(nstate,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Trial real-valued state `(nstate,npoints)` returned by this evaluation.
        real(wp), intent(out) :: response(:,:)
            !! Packed `traction.normal` response `(noutput,npoints)` for the evaluated material points.
        integer, intent(inout) :: point_status(:)
            !! Per-point status values; successful evaluations must set each entry to `cf_ok`.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact discrete response Jacobian `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(in), optional :: selected_parameters(:)
            !! One-based IDs of selected model parameters, in sensitivity-axis order.
        real(wp), intent(out), optional :: dsdp(:,:,:)
            !! Exact total response derivatives `(noutput,nselected,npoints)` at fixed prescribed input history.
        real(wp) :: opening
        real(wp) :: previous_max
        real(wp) :: maximum
        real(wp) :: effective_opening
        real(wp) :: opening_fraction
        real(wp) :: remaining_fraction
        real(wp) :: difference
        real(wp) :: peak_traction
        real(wp) :: softening_slope
        real(wp) :: tolerance
        real(wp) :: derivative
        logical :: at_onset
        logical :: at_failure
        logical :: at_maximum
        integer :: ipoint
        integer :: iparam

        difference = self%failure_opening-self%onset_opening
        peak_traction = self%stiffness*self%onset_opening
        softening_slope = peak_traction/difference
        tolerance = 32.0_wp*epsilon(1.0_wp)
        point_loop: do ipoint = 1, size(x_new,2)
            point_status(ipoint) = cf_invalid_argument
            opening = x_new(1,ipoint)
            previous_max = state_old(1,ipoint)
            if (.not. ieee_is_finite(opening)) cycle point_loop
            if (.not. ieee_is_finite(previous_max)) cycle point_loop
            if (opening < 0.0_wp .or. previous_max < 0.0_wp) cycle point_loop
            maximum = max(previous_max,opening)
            at_onset = abs(maximum-self%onset_opening)/max(maximum,self%onset_opening) <= tolerance
            at_failure = abs(maximum-self%failure_opening)/max(maximum,self%failure_opening) <= tolerance
            at_maximum = .true.
            if (maximum > 0.0_wp) at_maximum = abs(opening-maximum)/maximum <= tolerance
            if (present(selected_parameters)) then
                if (((at_onset .and. any(selected_parameters == 2)) .or. (at_failure .and. any(selected_parameters == 3))) .and. opening > 0.0_wp) then
                    point_status(ipoint) = cf_nondifferentiable
                    cycle point_loop
                end if
            end if
            if (present(tangent)) then
                if (previous_max < self%failure_opening .and. (at_onset .or. at_failure) .and. at_maximum) then
                    point_status(ipoint) = cf_nondifferentiable
                    cycle point_loop
                end if
                if (previous_max > self%onset_opening .and. previous_max < self%failure_opening) then
                    if (abs(opening-previous_max)/maximum <= tolerance) then
                        point_status(ipoint) = cf_nondifferentiable
                        cycle point_loop
                    end if
                end if
            end if

            effective_opening = opening
            response(1,ipoint) = 0.0_wp
            if (maximum > self%onset_opening .and. maximum < self%failure_opening) then
                opening_fraction = opening/maximum
                remaining_fraction = (self%failure_opening-maximum)/difference
                effective_opening = (self%onset_opening*opening_fraction)*remaining_fraction
                response(1,ipoint) = (peak_traction*opening_fraction)*remaining_fraction
            else if (maximum >= self%failure_opening) then
                effective_opening = 0.0_wp
            else
                response(1,ipoint) = self%stiffness*opening
            end if
            if (.not. ieee_is_finite(response(1,ipoint))) cycle point_loop

            if (present(tangent)) then
                derivative = self%stiffness
                if (maximum > self%onset_opening .and. maximum < self%failure_opening) then
                    derivative = (peak_traction/maximum)*remaining_fraction
                    if (opening > previous_max) derivative = -softening_slope
                else if (maximum >= self%failure_opening) then
                    derivative = 0.0_wp
                end if
                if (.not. ieee_is_finite(derivative)) cycle point_loop
                tangent(1,1,ipoint) = derivative
            end if

            if (present(selected_parameters)) then
                parameter_loop: do iparam = 1, size(selected_parameters)
                    select case (selected_parameters(iparam))
                    case (1)
                        dsdp(1,iparam,ipoint) = effective_opening
                    case (2)
                        dsdp(1,iparam,ipoint) = 0.0_wp
                        if (maximum > self%onset_opening .and. maximum < self%failure_opening) then
                            dsdp(1,iparam,ipoint) = ((self%stiffness*opening_fraction)*remaining_fraction)*(self%failure_opening/difference)
                        end if
                    case (3)
                        dsdp(1,iparam,ipoint) = 0.0_wp
                        if (maximum > self%onset_opening .and. maximum < self%failure_opening) then
                            dsdp(1,iparam,ipoint) = (softening_slope*opening_fraction)*((maximum-self%onset_opening)/difference)
                        end if
                    case default
                        exit parameter_loop
                    end select
                end do parameter_loop
                if (iparam <= size(selected_parameters)) cycle point_loop
                if (.not. all(ieee_is_finite(dsdp(:,:,ipoint)))) cycle point_loop
            end if

            state_trial(1,ipoint) = maximum
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_points

end module constituf_cohesive_mode_i

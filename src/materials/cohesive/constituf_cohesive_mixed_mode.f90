! SPDX-License-Identifier: Apache-2.0

module constituf_cohesive_mixed_mode
    !! Three-dimensional, rate-independent mixed-mode cohesive interface.
    !!
    !! The input is the interface-local displacement jump
    !! \(\boldsymbol{u}=(u_n,u_s,u_t)\) (length units), with \(u_n>0\) in opening.
    !! Let \(\langle u_n\rangle_+=\max(u_n,0)\) and define the effective
    !! separation and committed-history update by
    !! \[\delta_{\mathrm{eff}}=
    !! \sqrt{\langle u_n\rangle_+^2+\beta^2(u_s^2+u_t^2)},
    !! \qquad
    !! \kappa_{n+1}=\max(\kappa_n,\delta_{\mathrm{eff}}).\]
    !! With penalty stiffness \(K>0\) (stress per length), onset separation
    !! \(0<\delta_0<\delta_f\) (length units) and dimensionless shear weight
    !! \(\beta>0\), the surviving stiffness fraction is
    !! \[r(\kappa)=1-D(\kappa)=
    !! \begin{cases}
    !! 1,&0\le\kappa\le\delta_0,\\
    !! \displaystyle\frac{\delta_0(\delta_f-\kappa)}
    !! {\kappa(\delta_f-\delta_0)},&
    !! \delta_0<\kappa<\delta_f,\\
    !! 0,&\kappa\ge\delta_f.
    !! \end{cases}\]
    !! The normal traction uses an undamaged, frictionless compressive
    !! contact penalty. Shear tractions retain the damage factor:
    !! \[\boldsymbol{t}
    !! =K\begin{bmatrix}
    !! r(\kappa)\langle u_n\rangle_++\min(u_n,0)\\
    !! r(\kappa)\beta^2u_s\\
    !! r(\kappa)\beta^2u_t
    !! \end{bmatrix}.\]
    !! Thus normal compression *alone* does not increase damage, but
    !! tangential separation can increase \(\kappa\) even under
    !! compression. The shear law is an elastic slip penalty and does
    !! not implement Coulomb friction. For proportional opening/shear
    !! separation, the damage-controlled work of separation is
    !! \(G_c=K\delta_0\delta_f/2\) (energy per area). No mode-mix-dependent
    !! fracture-energy criterion is implemented. Recoverable normal
    !! compression work is separate from this fracture energy.
    !!
    !! For loading in the smooth softening range with positive normal
    !! opening, let \(\boldsymbol{D}_\beta=
    !! \mathrm{diag}(1,\beta^2,\beta^2)\). Then the exact local
    !! algorithmic tangent is
    !! \[\boldsymbol{T}
    !! =K r\boldsymbol{D}_\beta+
    !! \frac{K r'(\delta_{\mathrm{eff}})}{\delta_{\mathrm{eff}}}
    !! (\boldsymbol{D}_\beta\boldsymbol{u})\otimes
    !! (\boldsymbol{D}_\beta\boldsymbol{u}),\qquad
    !! r'(\delta)=-\frac{\delta_0\delta_f}
    !! {(\delta_f-\delta_0)\delta^2}.\]
    !! For numerical evaluation in opening, write
    !! \(\boldsymbol{B}=\operatorname{diag}(1,\beta,\beta)\),
    !! \(\boldsymbol{v}=\boldsymbol{B}\boldsymbol{u}\) and
    !! \(\boldsymbol{n}=\boldsymbol{v}/\delta_{\mathrm{eff}}\).
    !! The equivalent loading tangent is
    !! \[\boldsymbol{T}=\boldsymbol{B}
    !! \left[Kr(\boldsymbol{I}-\boldsymbol{n}\boldsymbol{n}^T)
    !! -S\boldsymbol{n}\boldsymbol{n}^T\right]\boldsymbol{B},
    !! \qquad S=\frac{K\delta_0}{\delta_f-\delta_0}.\]
    !! This form retains the negative radial softening slope directly.
    !! During unloading, the tangent is \(Kr\boldsymbol{D}_\beta\).
    !! Under compression the normal diagonal is instead \(K\), with
    !! the normal component excluded from the effective separation.
    !! At damage onset/completion, the damaged contact switch and
    !! softening loading/unloading switches, the classical tangent can
    !! be nonunique. Tangent requests report `cf_nondifferentiable`.
    !! Traction-only evaluation remains defined at these points.
    !! Strict unloading below a committed onset/completion history has
    !! a unique opening tangent. A previously fully failed interface
    !! retains only the compressive contact tangent.
    !!
    !! Input `opening.jump` and output `traction.jump` are three-component
    !! vectors in the interface-local frame, with length and stress dimensions respectively.
    !! The real history `opening.maximum` is the maximum effective
    !! separation (length units). Only accepted trials update committed history.
    !! Initial history and its parameter derivatives are zero.
    !! Exact parameter sensitivities are tracked for stiffness, onset and
    !! failure openings and shear weight. In particular, the derivative
    !! of the maximum effective-opening history includes the shear weight.
    !! For fixed prescribed jumps, the current effective-opening derivative
    !! is zero except for \(\beta\), for which
    !! \[\frac{\partial\delta_{\mathrm{eff}}}{\partial\beta}
    !! =\frac{\beta(u_s^2+u_t^2)}{\delta_{\mathrm{eff}}},
    !! \qquad\delta_{\mathrm{eff}}>0.\]
    !! Loading adopts this derivative and unloading
    !! retains the committed derivative. At equal competing maxima the
    !! history derivative exists only when both candidates agree for
    !! every selected parameter. Otherwise tracking returns
    !! `cf_nondifferentiable`, even in an intact or fully failed regime.
    !! At \(\kappa=\delta_0\) or \(\kappa=\delta_f\), a selected parameter derivative of
    !! cohesive traction is unique only when \(\mathrm{d}\kappa/\mathrm{d}p\) equals the respective
    !! threshold derivative, unless the damaged-driving jump is zero.
    !! This condition also applies during unloading. The contact switch
    !! alone does not make fixed-jump parameter derivatives nonunique.
    !!
    !! Adapted framework: P. P. Camanho and C. G. Davila,
    !! *Mixed-Mode Decohesion Finite Elements for the Simulation of
    !! Delamination in Composite Materials*, NASA/TM-2002-211737 (2002).
    !! Sections 2.2.1--2.2.2, contact equations (26)/(43), bilinear
    !! damage (44) and maximum history (45):
    !! https://ntrs.nasa.gov/api/citations/20020053651/downloads/20020053651.pdf .
    !! This independently weighted adaptation uses fixed \(\delta_0\) and
    !! \(\delta_f\) with shear stiffness \(K\beta^2\). The report instead uses an
    !! unweighted Euclidean separation, equal mode penalties and
    !! mode-dependent onset/failure thresholds. Its \(\beta\) in equation (32)
    !! denotes the current mode-mix ratio, not this fixed shear parameter.
    !! Here the pure normal/shear strengths are \(K\delta_0\) and
    !! \(K\beta\delta_0\), the failure jumps are \(\delta_f\) and \(\delta_f/\beta\)
    !! and both modes have fracture energy \(K\delta_0\delta_f/2\).
    !! No Benzeggagh--Kenane or other mode-dependent energy law is used.
    !! Code provenance: this weighted specialization, discrete update and
    !! analytical derivatives are independently derived from the stated
    !! equations. No third-party software code is reused in this revision.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument, cf_nondifferentiable
    implicit none (type, external)
    private

    integer, parameter :: nvec=3
    public :: cohesive_mixed_mode_t

    type, extends(material_model_t) :: cohesive_mixed_mode_t
        !! Construct as `cohesive_mixed_mode_t(K,d0,df,beta)`.
        real(wp) :: stiffness=0.0_wp
            !! Undamaged normal stiffness \(K>0\), stress per length.
        real(wp) :: onset_opening=0.0_wp
            !! Effective opening at damage onset \(\delta_0>0\), in the chosen length unit.
        real(wp) :: failure_opening=0.0_wp
            !! Effective opening at complete damage \(\delta_f>\delta_0\), in the chosen length unit.
        real(wp) :: shear_weight=0.0_wp
            !! \(\beta>0\): initial shear/normal stiffness ratio is \(\beta^2\).
    contains
        procedure :: describe => describe_mixed
            !! Declare the constitutive fields, parameter metadata and supported capabilities.
        procedure :: eval_batch => evaluate_mixed
            !! Evaluate the batched constitutive response, optional tangent and trial history.
        procedure :: initialize_sensitivity => initialize_mixed_sensitivity
            !! Initialize analytical sensitivities of the reference material history.
        procedure :: eval_sensitivity_batch => sensitivity_mixed
            !! Evaluate exact parameter-path derivatives and corresponding trial history.
    end type cohesive_mixed_mode_t

contains

    subroutine describe_mixed(self,spec,status)
        !! Validate the parameters and declare `opening.jump` input and `traction.jump` output with their physical metadata.
        class(cohesive_mixed_mode_t),intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        type(material_spec_t),intent(out) :: spec
            !! Material schema receiving field layouts, parameters and supported operations.
        integer,intent(out) :: status
            !! Operation status. `cf_ok` indicates a valid model or successful initialization.

        status=cf_invalid_argument
        if (.not. all(ieee_is_finite([self%stiffness,self%onset_opening, self%failure_opening,self%shear_weight]))) return
        if (self%stiffness<=0.0_wp.or.self%onset_opening<=0.0_wp) return
        if (self%failure_opening<=self%onset_opening.or.self%shear_weight<=0.0_wp) return
        if (.not. ieee_is_finite((self%stiffness*self%shear_weight)*self%shear_weight)) return
        if (.not. ieee_is_finite(self%stiffness*self%failure_opening)) return

        call spec%add_input("opening.jump", nvec, "vector", "displacement_jump", "interface", [nvec])
        call spec%add_output("traction.jump", nvec, "vector", "interface_traction", "interface", [nvec])
        call spec%add_state("opening.maximum", 1, "scalar", "effective_opening", "interface")
        call spec%add_parameter("stiffness")
        call spec%add_parameter("onset_opening")
        call spec%add_parameter("failure_opening")
        call spec%add_parameter("shear_weight")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status=cf_ok
    end subroutine describe_mixed

    subroutine evaluate_mixed(self,x_old,x_new,state_old,state_trial,response,dt,tangent,point_status)
        !! Exact trial maximum opening and traction with optional \(3\times3\) tangent.
        class(cohesive_mixed_mode_t),intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        real(wp),intent(in) :: x_old(:,:)
            !! Previously accepted packed input `(ninput,npoints)` if required, otherwise an empty array.
        real(wp),intent(in) :: x_new(:,:)
            !! Prescribed `opening.jump` input `(ninput,npoints)`.
        real(wp),intent(in) :: state_old(:,:)
            !! Committed real-valued state `(nstate,npoints)`.
        real(wp),intent(out) :: state_trial(:,:)
            !! Trial real-valued state `(nstate,npoints)` returned by this evaluation.
        real(wp),intent(out) :: response(:,:)
            !! Packed `traction.jump` response `(noutput,npoints)` for the evaluated material points.
        real(wp),intent(in) :: dt
            !! Prescribed time increment. It is positive when the material requires time integration.
        real(wp),intent(out),optional :: tangent(:,:,:)
            !! Optional exact discrete response Jacobian `(noutput,ninput,npoints)` at fixed committed history.
        integer,intent(inout) :: point_status(:)
            !! Per-point status values; successful evaluations must set each entry to `cf_ok`.
        real(wp) :: u(nvec)
        real(wp) :: scaled_jump(nvec)
        real(wp) :: direction(nvec)
        real(wp) :: transverse_projection(nvec)
        real(wp) :: axis_scale(nvec)
        real(wp) :: effective
        real(wp) :: kappa
        real(wp) :: old_kappa
        real(wp) :: difference
        real(wp) :: peak_traction
        real(wp) :: remaining_fraction
        real(wp) :: degraded_stiffness
        real(wp) :: softening_slope
        real(wp) :: projection
        real(wp) :: tangent_entry
        real(wp) :: tolerance
        integer :: ipoint
        integer :: i
        integer :: j
        logical :: loading
        logical :: softening
        logical :: at_onset
        logical :: at_failure
        logical :: at_current

        state_trial=state_old
        difference=self%failure_opening-self%onset_opening
        peak_traction=self%stiffness*self%onset_opening
        axis_scale=[1.0_wp,self%shear_weight,self%shear_weight]
        tolerance=32.0_wp*epsilon(1.0_wp)
        point_loop: do ipoint=1,size(x_new,2)
            point_status(ipoint)=cf_invalid_argument
            u=x_new(:,ipoint)
            old_kappa=state_old(1,ipoint)
            if (.not. all(ieee_is_finite(u))) cycle point_loop
            if (.not. ieee_is_finite(old_kappa)) cycle point_loop
            if (old_kappa<0.0_wp) cycle point_loop

            scaled_jump=axis_scale*u
            scaled_jump(1)=max(u(1),0.0_wp)
            if (.not. all(ieee_is_finite(scaled_jump))) cycle point_loop
            effective=opening_norm(scaled_jump)
            if (.not. ieee_is_finite(effective)) cycle point_loop
            kappa=max(old_kappa,effective)
            softening=kappa>self%onset_opening.and.kappa<self%failure_opening
            loading=effective>old_kappa
            degraded_stiffness=self%stiffness
            if (kappa<=self%onset_opening) then
                response(:,ipoint)=(self%stiffness*axis_scale)*scaled_jump
            else if (kappa>=self%failure_opening) then
                degraded_stiffness=0.0_wp
                response(:,ipoint)=0.0_wp
            else
                remaining_fraction=(self%failure_opening-kappa)/difference
                degraded_stiffness=(peak_traction/kappa)*remaining_fraction
                response(:,ipoint)=((peak_traction*(scaled_jump/kappa))*remaining_fraction)*axis_scale
            end if
            if (u(1)<0.0_wp) response(1,ipoint)=self%stiffness*u(1)
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop

            if (present(tangent)) then
                at_onset=abs(kappa-self%onset_opening)/max(kappa,self%onset_opening)<=tolerance
                at_failure=abs(kappa-self%failure_opening)/max(kappa,self%failure_opening)<=tolerance
                at_current=.true.
                if (kappa>0.0_wp) at_current=abs(effective-kappa)/kappa<=tolerance
                if (old_kappa<self%failure_opening.and.(at_onset.or.at_failure).and.at_current) then
                    point_status(ipoint)=cf_nondifferentiable
                    cycle point_loop
                end if
                if (abs(u(1))<=0.0_wp.and.kappa>self%onset_opening) then
                    point_status(ipoint)=cf_nondifferentiable
                    cycle point_loop
                end if
                if (softening) then
                    if (abs(effective-old_kappa)/kappa<=tolerance) then
                        point_status(ipoint)=cf_nondifferentiable
                        cycle point_loop
                    end if
                end if
                tangent(:,:,ipoint)=0.0_wp
                do i=1,nvec
                    tangent(i,i,ipoint)=(degraded_stiffness*axis_scale(i))*axis_scale(i)
                end do
                if (loading.and.softening) then
                    direction=scaled_jump/effective
                    transverse_projection=[direction(2)**2+direction(3)**2,direction(1)**2+direction(3)**2,direction(1)**2+direction(2)**2]
                    softening_slope=peak_traction/difference
                    do j=1,nvec
                        do i=1,nvec
                            projection=direction(i)*direction(j)
                            tangent_entry=-softening_slope*projection
                            if (i==j) then
                                tangent_entry=tangent_entry+degraded_stiffness*transverse_projection(i)
                            else
                                tangent_entry=tangent_entry-degraded_stiffness*projection
                            end if
                            tangent(i,j,ipoint)=(tangent_entry*axis_scale(i))*axis_scale(j)
                        end do
                    end do
                end if
                if (u(1)<0.0_wp) tangent(1,1,ipoint)=self%stiffness
                if (.not. all(ieee_is_finite(tangent(:,:,ipoint)))) cycle point_loop
            end if
            state_trial(1,ipoint)=kappa
            point_status(ipoint)=cf_ok
        end do point_loop
    end subroutine evaluate_mixed


    subroutine initialize_mixed_sensitivity(self,initial_input,initial_state,selected_parameters,sensitivity,status)
        !! Virgin interface history is zero, independent of material parameters.
        class(cohesive_mixed_mode_t), intent(in) :: self
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

        sensitivity=0.0_wp
        status=cf_ok
    end subroutine initialize_mixed_sensitivity

    subroutine sensitivity_mixed(self,x_old,x_new,state_old,state_trial,selected_parameters, sensitivity_old,sensitivity_trial,response,dsdp,dt,tangent,point_status)
        !! Exact total derivative through \(\max(\kappa_n,\delta_{\mathrm{eff}})\) and the
        !! bilinear surviving fraction. The driving jumps remain prescribed.
        class(cohesive_mixed_mode_t), intent(in) :: self
            !! Material model instance containing the declared constitutive parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed input `(ninput,npoints)` if required, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Prescribed `opening.jump` input `(ninput,npoints)`.
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
            !! Packed `traction.jump` response `(noutput,npoints)` for the evaluated material points.
        real(wp), intent(out) :: dsdp(:,:,:)
            !! Exact total response derivatives `(noutput,nselected,npoints)` at fixed prescribed input history.
        real(wp), intent(in) :: dt
            !! Prescribed time increment. It is positive when the material requires time integration.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact discrete response Jacobian `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(inout) :: point_status(:)
            !! Per-point status values; successful evaluations must set each entry to `cf_ok`.
        real(wp) :: u(nvec)
        real(wp) :: scaled_jump(nvec)
        real(wp) :: direction(nvec)
        real(wp) :: weighted_fraction(nvec)
        real(wp) :: traction_derivative(nvec)
        real(wp) :: effective
        real(wp) :: kappa
        real(wp) :: kappa_derivative
        real(wp) :: candidate_derivative
        real(wp) :: derivative_scale
        real(wp) :: threshold_derivative
        real(wp) :: peak_traction
        real(wp) :: remaining_fraction
        real(wp) :: history_factor
        real(wp) :: softening_slope
        real(wp) :: d0
        real(wp) :: df
        real(wp) :: beta
        real(wp) :: difference
        real(wp) :: d_stiffness
        real(wp) :: d_onset
        real(wp) :: d_failure
        real(wp) :: d_beta
        real(wp) :: tolerance
        integer :: ipoint
        integer :: iparam
        logical :: loading
        logical :: history_tie
        logical :: at_onset
        logical :: at_failure

        call evaluate_mixed(self,x_old,x_new,state_old,state_trial,response,dt,tangent,point_status)
        sensitivity_trial=sensitivity_old
        dsdp=0.0_wp
        d0=self%onset_opening
        df=self%failure_opening
        beta=self%shear_weight
        difference=df-d0
        peak_traction=self%stiffness*d0
        tolerance=32.0_wp*epsilon(1.0_wp)
        point_loop: do ipoint=1,size(x_new,2)
            if (point_status(ipoint)/=cf_ok) cycle point_loop
            u=x_new(:,ipoint)
            scaled_jump=[max(u(1),0.0_wp),beta*u(2),beta*u(3)]
            effective=opening_norm(scaled_jump)
            kappa=state_trial(1,ipoint)
            at_onset=abs(kappa-d0)/max(kappa,d0)<=tolerance
            at_failure=abs(kappa-df)/max(kappa,df)<=tolerance
            if (at_onset.and.at_failure) then
                ! Preserve endpoint identity when the proximity bands overlap.
                at_onset=abs(kappa-d0)<=abs(kappa-df)
                at_failure=.not.at_onset
            end if
            history_tie=.true.
            if (kappa>0.0_wp) history_tie=abs(effective-state_old(1,ipoint))/kappa<=tolerance
            loading=effective>state_old(1,ipoint)
            direction=0.0_wp
            if (effective>0.0_wp) direction=scaled_jump/effective
            if (kappa>d0.and.kappa<df) then
                remaining_fraction=(df-kappa)/difference
                softening_slope=peak_traction/difference
                weighted_fraction=(scaled_jump/kappa)*[1.0_wp,beta,beta]
            end if
            parameter_loop: do iparam=1,size(selected_parameters)
                d_stiffness=0.0_wp
                d_onset=0.0_wp
                d_failure=0.0_wp
                d_beta=0.0_wp
                select case(selected_parameters(iparam))
                case(1)
                    d_stiffness=1.0_wp
                case(2)
                    d_onset=1.0_wp
                case(3)
                    d_failure=1.0_wp
                case(4)
                    d_beta=1.0_wp
                case default
                    point_status(ipoint)=cf_invalid_argument
                    exit parameter_loop
                end select

                candidate_derivative=0.0_wp
                if (d_beta>0.0_wp) candidate_derivative=dot_product(direction(2:nvec),u(2:nvec))
                if (history_tie) then
                    derivative_scale=max(abs(candidate_derivative),abs(sensitivity_old(1,iparam,ipoint)))
                    if (derivative_scale>0.0_wp) then
                        if (abs(candidate_derivative-sensitivity_old(1,iparam,ipoint))/derivative_scale>tolerance) then
                            point_status(ipoint)=cf_nondifferentiable
                            cycle point_loop
                        end if
                    end if
                end if
                kappa_derivative=sensitivity_old(1,iparam,ipoint)
                if (loading) kappa_derivative=candidate_derivative
                sensitivity_trial(1,iparam,ipoint)=kappa_derivative

                if ((at_onset.or.at_failure).and.effective>0.0_wp) then
                    threshold_derivative=d_onset
                    if (at_failure) threshold_derivative=d_failure
                    derivative_scale=max(abs(kappa_derivative),abs(threshold_derivative))
                    if (derivative_scale>0.0_wp) then
                        if (abs(kappa_derivative-threshold_derivative)/derivative_scale>tolerance) then
                            point_status(ipoint)=cf_nondifferentiable
                            cycle point_loop
                        end if
                    end if
                end if

                traction_derivative=0.0_wp
                if (effective>0.0_wp) then
                    if (kappa<=d0) then
                        if (d_stiffness>0.0_wp) traction_derivative=scaled_jump*[1.0_wp,beta,beta]
                        if (d_beta>0.0_wp) traction_derivative(2:nvec)=2.0_wp*(self%stiffness*scaled_jump(2:nvec))
                    else if (kappa<df) then
                        history_factor=-(peak_traction*(kappa_derivative/kappa))*(df/difference)
                        traction_derivative=history_factor*weighted_fraction
                        if (d_stiffness>0.0_wp) traction_derivative=traction_derivative+(d0*remaining_fraction)*weighted_fraction
                        if (d_onset>0.0_wp) traction_derivative=traction_derivative+((self%stiffness*remaining_fraction)*weighted_fraction)*(df/difference)
                        if (d_failure>0.0_wp) traction_derivative=traction_derivative+(softening_slope*weighted_fraction)*((kappa-d0)/difference)
                        if (d_beta>0.0_wp) traction_derivative(2:nvec)=traction_derivative(2:nvec)+(peak_traction*(scaled_jump(2:nvec)/kappa))*(2.0_wp*remaining_fraction)
                    end if
                end if
                if (u(1)<0.0_wp) traction_derivative(1)=d_stiffness*u(1)
                dsdp(:,iparam,ipoint)=traction_derivative
            end do parameter_loop
            if (point_status(ipoint)/=cf_ok) cycle point_loop
            if (.not. all(ieee_is_finite(dsdp(:,:,ipoint))) .or. .not. all(ieee_is_finite(sensitivity_trial(:,:,ipoint)))) then
                point_status(ipoint)=cf_invalid_argument
            end if
        end do point_loop
    end subroutine sensitivity_mixed

    pure function opening_norm(jump) result(effective)
        !! Euclidean length without squaring dimensional jump magnitudes.
        real(wp), intent(in) :: jump(nvec)
            !! Finite normal-positive and shear-weighted separation vector.
        real(wp) :: effective
        real(wp) :: scale

        scale=maxval(abs(jump))
        effective=0.0_wp
        if (scale>0.0_wp) effective=scale*sqrt(sum((jump/scale)**2))
    end function opening_norm

end module constituf_cohesive_mixed_mode

! SPDX-License-Identifier: Apache-2.0

module constituf_drucker_prager
    !! Three-dimensional, small-strain, associated Drucker--Prager plasticity.
    !!
    !! Cauchy stress \(\boldsymbol{\sigma}\) is tension positive (stress units).
    !! Infinitesimal strain \(\boldsymbol{\varepsilon}\) and stress are symmetric
    !! tensors. For either tensor \(\boldsymbol{A}\), the Mandel vector is
    !! \[\operatorname{mandel}(\boldsymbol{A})=
    !! [A_{11},A_{22},A_{33},\sqrt{2}A_{23},
    !! \sqrt{2}A_{13},\sqrt{2}A_{12}]^{\mathsf T}.\]
    !! Define the tensile-positive mean stress \(p=
    !! \mathrm{tr}(\boldsymbol{\sigma})/3\),
    !! \(\boldsymbol{s}=\mathrm{dev}(\boldsymbol{\sigma})\) and
    !! \(q=\sqrt{3/2}\|\boldsymbol{s}\|\). The yield condition is
    !! \[f=q+\beta p-(\sigma_y+H\alpha)\le0.\]
    !! Here \(\beta\ge0\) is a dimensionless pressure slope.
    !! \(\sigma_y>0\) is the initial yield intercept (stress units).
    !! \(H\ge0\) is the linear hardening modulus (stress units). And
    !! \(\alpha\) is the accumulated plastic multiplier (dimensionless).
    !! This sign convention makes yielding more accessible under
    !! hydrostatic tension. The cone opens toward hydrostatic compression.
    !!
    !! With elastic bulk modulus \(K>0\) and shear modulus \(G>0\)
    !! (both in the chosen stress unit), the associative flow direction on the smooth cone is
    !! \[\boldsymbol{n}=\frac{3\boldsymbol{s}}{2q},\qquad
    !! \Delta\boldsymbol{\varepsilon}^p
    !! =\Delta\gamma\left(\boldsymbol{n}
    !! +\frac{\beta}{3}\boldsymbol{I}\right),\qquad
    !! \Delta\alpha=\Delta\gamma.\]
    !! For the trial quantities \(p_{\mathrm{tr}}\),
    !! \(q_{\mathrm{tr}}\) and \(\boldsymbol{s}_{\mathrm{tr}}\),
    !! the smooth-cone backward-Euler return uses
    !! \[\Delta\gamma=
    !! \frac{q_{\mathrm{tr}}+\beta p_{\mathrm{tr}}
    !! -\sigma_y-H\alpha_n}{3G+K\beta^2+H},\]
    !! \[\boldsymbol{s}_{n+1}
    !! =\boldsymbol{s}_{\mathrm{tr}}
    !! -2G\Delta\gamma\boldsymbol{n},\qquad
    !! p_{n+1}=p_{\mathrm{tr}}-K\beta\Delta\gamma.\]
    !! The smooth-cone branch requires \(q_{\mathrm{tr}}>3G\Delta\gamma\).
    !! Otherwise the return is at the cone apex, where the algorithm
    !! sets \(\boldsymbol{s}_{n+1}=\boldsymbol{0}\) and uses
    !! \[\Delta\gamma=
    !! \frac{\beta p_{\mathrm{tr}}-\sigma_y-H\alpha_n}
    !! {K\beta^2+H},\]
    !! with \(\beta>0\) and pressure evaluated from consistency as
    !! \(p_{n+1}=[\sigma_y+H(\alpha_n+\Delta\gamma)]/\beta\).
    !! The apex plastic deviatoric strain is equal to the current
    !! total strain's deviator. The exact discrete algorithmic tangent
    !! is differentiated separately on elastic, smooth-cone and apex
    !! branches. Classical derivatives are not unique at elastic/plastic
    !! or smooth-cone/apex switches. Tangent and sensitivity requests at
    !! those switches report `cf_nondifferentiable`.
    !!
    !! History fields are `plastic.strain` (six dimensionless Mandel
    !! components) and `plastic.multiplier` (one dimensionless
    !! accumulated multiplier). A new handle starts with zero history
    !! and only accepted trials update it. Parameter sensitivities are
    !! computed analytically through accepted smooth-cone/apex branches,
    !! including accumulated plastic history. The model is not a
    !! Mohr--Coulomb law, a cap
    !! plasticity model, or a non-associated dilation rule.
    !!
    !! References: D. C. Drucker and W. Prager, *Soil Mechanics and
    !! Plastic Analysis or Limit Design*, Quarterly of Applied
    !! Mathematics **10** (1952), 157--165, DOI: 10.1090/qam/48291.
    !! Section 2, Eqs. (2)--(5), gives the yield condition and associated flow.
    !! Their \(aI_1+\sqrt{J_2}=k\) maps to
    !! \(\beta=3\sqrt3\,a\) and \(\sigma_y=\sqrt3\,k\).
    !! Linear isotropic hardening and the backward-Euler update are
    !! the present specialization. S. Sysala et al., *An improved
    !! return-mapping scheme for nonsmooth plastic potentials: PART I --
    !! the Haigh--Westergaard coordinates*, arXiv:1503.03605v2 (2015),
    !! Section 3, Eqs. (3.18)--(3.21) and (3.24)--(3.25), gives the
    !! smooth-cone/apex returns and their consistent tangents. The
    !! associated specialization uses \(\eta=\bar\eta=\beta/\sqrt3\),
    !! \(\xi=1/\sqrt3\), \(c_0=\sigma_y\) and
    !! \(\Delta\lambda_{\mathrm{ref}}=\sqrt3\Delta\gamma\), with
    !! \(\bar\varepsilon^p_{\mathrm{ref}}=\alpha\) and the linear hardening law.
    !! J. C. Simo and T. J. R. Hughes, *Computational Inelasticity*,
    !! Springer (1998), DOI: 10.1007/b98904. The discrete return algorithm
    !! is an independent implementation.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument, cf_nondifferentiable
    implicit none (type, external)
    private

    integer, parameter :: nmandel=6
    integer, parameter :: nstate=7
    public :: drucker_prager_t

    type, extends(material_model_t) :: drucker_prager_t
        !! Construct as `drucker_prager_t(K,G,sigma_y,beta,H)`, with consistent stress units.
        real(wp) :: bulk_modulus=0.0_wp
            !! Bulk modulus \(K>0\), in the chosen stress unit.
        real(wp) :: shear_modulus=0.0_wp
            !! Shear modulus \(G>0\), in the chosen stress unit.
        real(wp) :: yield_stress=0.0_wp
            !! Positive initial cohesion-like intercept \(\sigma_y\), in the chosen stress unit.
        real(wp) :: pressure_slope=0.0_wp
            !! Nonnegative Drucker--Prager slope \(\beta\), dimensionless.
        real(wp) :: hardening_modulus=0.0_wp
            !! Nonnegative isotropic hardening modulus \(H\), in the chosen stress unit.
    contains
        procedure :: describe => describe_dp
            !! Describe the material input, response, state and scalar parameter schema.
        procedure :: eval_batch => evaluate_dp
            !! Evaluate the batched material response and trial history.
        procedure :: initialize_sensitivity => initialize_dp_sensitivity
            !! Initialize analytical derivatives of the reference real history.
        procedure :: eval_sensitivity_batch => sensitivity_dp
            !! Evaluate exact selected-parameter sensitivities through the discrete material update.
    end type drucker_prager_t

contains

    subroutine describe_dp(self,spec,status)
        !! Validate model parameters and register `strain.small` input and `stress.cauchy` response fields.
        class(drucker_prager_t),intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        type(material_spec_t),intent(out) :: spec
            !! Material schema receiving the packed field descriptions and supported capabilities.
        integer,intent(out) :: status
            !! Status code; `cf_ok` indicates success.

        status=cf_invalid_argument
        if (.not. all(ieee_is_finite([self%bulk_modulus,self%shear_modulus, self%yield_stress,self%pressure_slope,self%hardening_modulus]))) return
        if (self%bulk_modulus<=0.0_wp.or.self%shear_modulus<=0.0_wp) return
        if (self%yield_stress<=0.0_wp.or.self%pressure_slope<0.0_wp.or.self%hardening_modulus<0.0_wp) return
        if (.not. ieee_is_finite(3.0_wp*self%shear_modulus+ self%bulk_modulus*self%pressure_slope**2+self%hardening_modulus)) return
        call spec%add_input("strain.small", nmandel, "mandel", "small_strain", "spatial", [3,3])
        call spec%add_output("stress.cauchy", nmandel, "mandel", "cauchy_stress", "spatial", [3,3])
        call spec%add_state("plastic.strain", nmandel, "mandel", "plastic_strain", "spatial", [3,3])
        call spec%add_state("plastic.multiplier", 1, "scalar", "accumulated_plastic_multiplier", "material")
        call spec%add_parameter("bulk_modulus")
        call spec%add_parameter("shear_modulus")
        call spec%add_parameter("yield_stress")
        call spec%add_parameter("pressure_slope")
        call spec%add_parameter("hardening_modulus")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status=cf_ok
    end subroutine describe_dp

    subroutine evaluate_dp(self,x_old,x_new,state_old,state_trial,response,dt,tangent,point_status)
        !! Trial update. Tangent is \(\partial\boldsymbol{\sigma}/\partial\boldsymbol{\varepsilon}\) at fixed committed state.
        class(drucker_prager_t),intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        real(wp),intent(in) :: x_old(:,:)
            !! Previously accepted packed inputs `(ninput,npoints)` when required, otherwise an empty array.
        real(wp),intent(in) :: x_new(:,:)
            !! Current prescribed `strain.small` inputs `(ninput,npoints)`.
        real(wp),intent(in) :: state_old(:,:)
            !! Committed real-valued history `(nstate,npoints)`.
        real(wp),intent(out) :: state_trial(:,:)
            !! Computed uncommitted real-valued history `(nstate,npoints)`.
        real(wp),intent(out) :: response(:,:)
            !! Packed `stress.cauchy` responses `(noutput,npoints)`.
        real(wp),intent(in) :: dt
            !! Prescribed time increment; a positive value is required for time-dependent laws.
        real(wp),intent(out),optional :: tangent(:,:,:)
            !! Optional exact algorithmic tangent `(noutput,ninput,npoints)` at fixed committed history.
        integer,intent(inout) :: point_status(:)
            !! Per-point outcome; write `cf_ok` only for a successful constitutive update.
        real(wp) :: projector(nmandel,nmandel)
        real(wp) :: identity(nmandel)
        real(wp) :: trial_dev(nmandel)
        real(wp) :: eps_dev(nmandel)
        real(wp) :: epsp_dev(nmandel)
        real(wp) :: direction(nmandel)
        real(wp) :: dev_final(nmandel)
        real(wp) :: trial_mean
        real(wp) :: qtrial
        real(wp) :: ftrial
        real(wp) :: yield_limit
        real(wp) :: denominator
        real(wp) :: apex_denominator
        real(wp) :: step
        real(wp) :: mean_final
        real(wp) :: dq
        real(wp) :: dlambda
        real(wp) :: ddev(nmandel)
        real(wp) :: dn(nmandel)
        real(wp) :: column(nmandel)
        real(wp) :: stress_scale
        real(wp) :: cone_correction
        real(wp) :: apex_scale
        integer :: ipoint
        integer :: i
        integer :: j
        logical :: plastic
        logical :: apex

        identity=[1.0_wp,1.0_wp,1.0_wp,0.0_wp,0.0_wp,0.0_wp]
        projector=0.0_wp
        do i=1,nmandel
            projector(i,i)=1.0_wp
        end do
        projector(1:3,1:3)=projector(1:3,1:3)-1.0_wp/3.0_wp
        denominator=3.0_wp*self%shear_modulus+ self%bulk_modulus*self%pressure_slope**2+self%hardening_modulus
        apex_denominator=self%bulk_modulus*self%pressure_slope**2+self%hardening_modulus
        state_trial=state_old

        point_loop: do ipoint=1,size(x_new,2)
            point_status(ipoint)=cf_invalid_argument
            if (.not. all(ieee_is_finite(x_new(:,ipoint)))) cycle point_loop
            if (.not. all(ieee_is_finite(state_old(:,ipoint)))) cycle point_loop
            if (state_old(nstate,ipoint)<0.0_wp) cycle point_loop

            eps_dev=x_new(:,ipoint)-sum(x_new(1:3,ipoint))*identity/3.0_wp
            epsp_dev=state_old(1:nmandel,ipoint)-sum(state_old(1:3,ipoint))*identity/3.0_wp
            trial_dev=2.0_wp*self%shear_modulus*(eps_dev-epsp_dev)
            trial_mean=self%bulk_modulus*(sum(x_new(1:3,ipoint))-sum(state_old(1:3,ipoint)))
            qtrial=sqrt(1.5_wp*dot_product(trial_dev,trial_dev))
            yield_limit=self%yield_stress+self%hardening_modulus*state_old(nstate,ipoint)
            ftrial=qtrial+self%pressure_slope*trial_mean-yield_limit
            if (.not. all(ieee_is_finite([trial_mean,qtrial,yield_limit,ftrial]))) cycle point_loop
            plastic=ftrial>0.0_wp
            apex=.false.
            step=0.0_wp
            direction=0.0_wp
            dev_final=trial_dev
            mean_final=trial_mean
            if (plastic) then
                step=ftrial/denominator
                apex=qtrial<=3.0_wp*self%shear_modulus*step
                if (apex) then
                    if (apex_denominator<=0.0_wp.or.self%pressure_slope<=0.0_wp) cycle point_loop
                    step=(self%pressure_slope*trial_mean-yield_limit)/apex_denominator
                    dev_final=0.0_wp
                    ! Consistency avoids cancellation between trial and correction pressures.
                    mean_final=(yield_limit+self%hardening_modulus*step)/self%pressure_slope
                else
                    direction=1.5_wp*trial_dev/qtrial
                    dev_final=trial_dev-2.0_wp*self%shear_modulus*step*direction
                    mean_final=trial_mean-self%bulk_modulus*self%pressure_slope*step
                end if
                if (.not. all(ieee_is_finite([step,mean_final]))) cycle point_loop
                if (step<0.0_wp) cycle point_loop
            end if

            response(:,ipoint)=dev_final+mean_final*identity
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop
            if (plastic) then
                state_trial(nstate,ipoint)=state_old(nstate,ipoint)+step
                if (apex) then
                    state_trial(1:nmandel,ipoint)=state_old(1:nmandel,ipoint)+ eps_dev-epsp_dev+self%pressure_slope*step*identity/3.0_wp
                else
                    state_trial(1:nmandel,ipoint)=state_old(1:nmandel,ipoint)+ step*(direction+self%pressure_slope*identity/3.0_wp)
                end if
                if (.not. all(ieee_is_finite(state_trial(:,ipoint)))) cycle point_loop
            end if

            if (present(tangent)) then
                stress_scale=max(qtrial,abs(self%pressure_slope*trial_mean),yield_limit)
                if (abs(ftrial)/stress_scale<=64.0_wp*epsilon(1.0_wp)) then
                    point_status(ipoint)=cf_nondifferentiable
                    cycle point_loop
                end if
                if (plastic) then
                    cone_correction=3.0_wp*self%shear_modulus*(ftrial/denominator)
                    apex_scale=max(qtrial,abs(cone_correction))
                    if (apex_scale<=0.0_wp) then
                        point_status(ipoint)=cf_nondifferentiable
                        cycle point_loop
                    end if
                    if (abs(qtrial-cone_correction)/apex_scale<=64.0_wp*epsilon(1.0_wp)) then
                        point_status(ipoint)=cf_nondifferentiable
                        cycle point_loop
                    end if
                end if

                do j=1,nmandel
                    ddev=2.0_wp*self%shear_modulus*projector(:,j)
                    column=ddev
                    if (j<=3) column=column+self%bulk_modulus*identity
                    if (plastic) then
                        if (apex) then
                            column=0.0_wp
                            if (j<=3) column=(self%bulk_modulus*(self%hardening_modulus/apex_denominator))*identity
                        else
                            dq=dot_product(direction,ddev)
                            dlambda=(dq+self%pressure_slope*self%bulk_modulus*identity(j))/denominator
                            dn=(1.5_wp/qtrial)*(ddev-trial_dev*dq/qtrial)
                            column=column-(2.0_wp*self%shear_modulus*direction+ self%bulk_modulus*self%pressure_slope*identity)*dlambda- 2.0_wp*self%shear_modulus*step*dn
                        end if
                    end if
                    tangent(:,j,ipoint)=column
                end do
                if (.not. all(ieee_is_finite(tangent(:,:,ipoint)))) cycle point_loop
            end if
            point_status(ipoint)=cf_ok
        end do point_loop
    end subroutine evaluate_dp


    subroutine initialize_dp_sensitivity(self,initial_input,initial_state,selected_parameters,sensitivity,status)
        !! Virgin zero plastic strain and multiplier are parameter independent.
        class(drucker_prager_t),intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        real(wp),intent(in) :: initial_input(:,:)
            !! Prescribed input state used for material initialization.
        real(wp),intent(in) :: initial_state(:,:)
            !! Initial real-valued history `(nstate,npoints)`.
        integer,intent(in) :: selected_parameters(:)
            !! One-based parameter IDs in the caller's selected derivative-axis order.
        real(wp),intent(out) :: sensitivity(:,:,:)
            !! Initial exact history derivatives `(nstate,nselected,npoints)`.
        integer,intent(out) :: status
            !! Status code; `cf_ok` indicates success.

        sensitivity=0.0_wp
        status=cf_ok
    end subroutine initialize_dp_sensitivity

    subroutine sensitivity_dp(self,x_old,x_new,state_old,state_trial,selected_parameters, sensitivity_old,sensitivity_trial,response,dsdp,dt,tangent,point_status)
        !! Exact derivative of the discrete elastic/cone/apex return through
        !! all committed plastic strain and multiplier sensitivities.
        class(drucker_prager_t),intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        real(wp),intent(in) :: x_old(:,:)
            !! Previously accepted packed inputs `(ninput,npoints)` when required, otherwise an empty array.
        real(wp),intent(in) :: x_new(:,:)
            !! Current prescribed `strain.small` inputs `(ninput,npoints)`.
        real(wp),intent(in) :: state_old(:,:)
            !! Committed real-valued history `(nstate,npoints)`.
        integer,intent(in) :: selected_parameters(:)
            !! One-based parameter IDs in the caller's selected derivative-axis order.
        real(wp),intent(in) :: sensitivity_old(:,:,:)
            !! Committed real-history parameter derivatives `(nstate,nselected,npoints)`.
        real(wp),intent(out) :: state_trial(:,:)
            !! Computed uncommitted real-valued history `(nstate,npoints)`.
        real(wp),intent(out) :: sensitivity_trial(:,:,:)
            !! Updated real-history parameter derivatives `(nstate,nselected,npoints)`.
        real(wp),intent(out) :: response(:,:)
            !! Packed `stress.cauchy` responses `(noutput,npoints)`.
        real(wp),intent(out) :: dsdp(:,:,:)
            !! Exact total response derivatives `(noutput,nselected,npoints)` at fixed prescribed driving histories.
        real(wp),intent(in) :: dt
            !! Prescribed time increment; a positive value is required for time-dependent laws.
        real(wp),intent(out),optional :: tangent(:,:,:)
            !! Optional exact algorithmic tangent `(noutput,ninput,npoints)` at fixed committed history.
        integer,intent(inout) :: point_status(:)
            !! Per-point outcome; write `cf_ok` only for a successful constitutive update.
        real(wp) :: identity(nmandel)
        real(wp) :: eps_dev(nmandel)
        real(wp) :: epsp_dev(nmandel)
        real(wp) :: trial_dev(nmandel)
        real(wp) :: direction(nmandel)
        real(wp) :: direction_derivative(nmandel)
        real(wp) :: previous_derivative(nmandel)
        real(wp) :: dev_derivative(nmandel)
        real(wp) :: stress_derivative(nmandel)
        real(wp) :: d_bulk
        real(wp) :: d_shear
        real(wp) :: d_yield
        real(wp) :: d_slope
        real(wp) :: d_hardening
        real(wp) :: mean_trial
        real(wp) :: qtrial
        real(wp) :: yield_limit
        real(wp) :: ftrial
        real(wp) :: denominator
        real(wp) :: apex_denominator
        real(wp) :: increment
        real(wp) :: di
        real(wp) :: dden
        real(wp) :: dmean
        real(wp) :: dq
        real(wp) :: dyield
        real(wp) :: dsdev(nmandel)
        real(wp) :: trace_strain
        real(wp) :: trace_old
        real(wp) :: slope
        real(wp) :: stress_scale
        real(wp) :: cone_correction
        real(wp) :: apex_scale
        integer :: ipoint
        integer :: iparam
        logical :: plastic
        logical :: apex

        call evaluate_dp(self,x_old,x_new,state_old,state_trial,response,dt,tangent,point_status)
        identity=[1.0_wp,1.0_wp,1.0_wp,0.0_wp,0.0_wp,0.0_wp]
        slope=self%pressure_slope
        denominator=3.0_wp*self%shear_modulus+self%bulk_modulus*slope**2+self%hardening_modulus
        apex_denominator=self%bulk_modulus*slope**2+self%hardening_modulus
        sensitivity_trial=sensitivity_old
        dsdp=0.0_wp

        point_loop: do ipoint=1,size(x_new,2)
            if (point_status(ipoint)/=cf_ok) cycle point_loop
            trace_strain=sum(x_new(1:3,ipoint))
            trace_old=sum(state_old(1:3,ipoint))
            eps_dev=x_new(:,ipoint)-trace_strain*identity/3.0_wp
            epsp_dev=state_old(1:nmandel,ipoint)-trace_old*identity/3.0_wp
            trial_dev=2.0_wp*self%shear_modulus*(eps_dev-epsp_dev)
            mean_trial=self%bulk_modulus*(trace_strain-trace_old)
            qtrial=sqrt(1.5_wp*dot_product(trial_dev,trial_dev))
            yield_limit=self%yield_stress+self%hardening_modulus*state_old(nstate,ipoint)
            ftrial=qtrial+slope*mean_trial-yield_limit
            plastic=ftrial>0.0_wp
            apex=.false.
            increment=0.0_wp
            direction=0.0_wp
            if (plastic) then
                increment=ftrial/denominator
                apex=qtrial<=3.0_wp*self%shear_modulus*increment
                if (apex) then
                    increment=(slope*mean_trial-yield_limit)/apex_denominator
                else
                    direction=1.5_wp*trial_dev/qtrial
                end if
            end if
            stress_scale=max(qtrial,abs(slope*mean_trial),yield_limit)
            if (abs(ftrial)/stress_scale<=64.0_wp*epsilon(1.0_wp)) then
                point_status(ipoint)=cf_nondifferentiable
                cycle point_loop
            end if
            if (plastic) then
                cone_correction=3.0_wp*self%shear_modulus*(ftrial/denominator)
                apex_scale=max(qtrial,abs(cone_correction))
                if (apex_scale<=0.0_wp) then
                    point_status(ipoint)=cf_nondifferentiable
                    cycle point_loop
                end if
                if (abs(qtrial-cone_correction)/apex_scale<=64.0_wp*epsilon(1.0_wp)) then
                    point_status(ipoint)=cf_nondifferentiable
                    cycle point_loop
                end if
            end if

            parameter_loop: do iparam=1,size(selected_parameters)
                d_bulk=0.0_wp
                d_shear=0.0_wp
                d_yield=0.0_wp
                d_slope=0.0_wp
                d_hardening=0.0_wp
                select case(selected_parameters(iparam))
                case(1)
                    d_bulk=1.0_wp
                case(2)
                    d_shear=1.0_wp
                case(3)
                    d_yield=1.0_wp
                case(4)
                    d_slope=1.0_wp
                case(5)
                    d_hardening=1.0_wp
                case default
                    point_status(ipoint)=cf_invalid_argument
                    exit parameter_loop
                end select

                previous_derivative=sensitivity_old(1:nmandel,iparam,ipoint)
                dev_derivative=previous_derivative-sum(previous_derivative(1:3))*identity/3.0_wp
                dsdev=2.0_wp*d_shear*(eps_dev-epsp_dev)- 2.0_wp*self%shear_modulus*dev_derivative
                dmean=d_bulk*(trace_strain-trace_old)- self%bulk_modulus*sum(previous_derivative(1:3))
                dyield=d_yield+d_hardening*state_old(nstate,ipoint)+ self%hardening_modulus*sensitivity_old(nstate,iparam,ipoint)
                if (plastic) then
                    if (apex) then
                        dden=d_bulk*slope**2+2.0_wp*self%bulk_modulus*slope*d_slope+d_hardening
                        di=(d_slope*mean_trial+slope*dmean-dyield-increment*dden)/apex_denominator
                        stress_derivative=0.0_wp
                        stress_derivative(1:3)=(dyield+d_hardening*increment+self%hardening_modulus*di-response(1,ipoint)*d_slope)/slope
                        sensitivity_trial(1:nmandel,iparam,ipoint)=previous_derivative-dev_derivative+ (d_slope*increment+slope*di)*identity/3.0_wp
                    else
                        dq=dot_product(direction,dsdev)
                        dden=3.0_wp*d_shear+d_bulk*slope**2+ 2.0_wp*self%bulk_modulus*slope*d_slope+d_hardening
                        di=(dq+d_slope*mean_trial+slope*dmean-dyield-increment*dden)/denominator
                        direction_derivative=(1.5_wp/qtrial)*(dsdev-trial_dev*dq/qtrial)
                        stress_derivative=dsdev-2.0_wp*d_shear*increment*direction- 2.0_wp*self%shear_modulus*(di*direction+increment*direction_derivative)
                        stress_derivative(1:3)=stress_derivative(1:3)+ dmean-d_bulk*slope*increment- self%bulk_modulus*(d_slope*increment+slope*di)
                        sensitivity_trial(1:nmandel,iparam,ipoint)=previous_derivative+ di*(direction+slope*identity/3.0_wp)+ increment*(direction_derivative+d_slope*identity/3.0_wp)
                    end if
                    sensitivity_trial(nstate,iparam,ipoint)= sensitivity_old(nstate,iparam,ipoint)+di
                else
                    stress_derivative=dsdev+dmean*identity
                end if
                dsdp(:,iparam,ipoint)=stress_derivative
            end do parameter_loop
            if (point_status(ipoint)/=cf_ok) cycle point_loop
            if (.not. all(ieee_is_finite(dsdp(:,:,ipoint))) .or. .not. all(ieee_is_finite(sensitivity_trial(:,:,ipoint)))) then
                point_status(ipoint)=cf_invalid_argument
            end if
        end do point_loop
    end subroutine sensitivity_dp

end module constituf_drucker_prager

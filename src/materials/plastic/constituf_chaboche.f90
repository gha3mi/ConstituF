! SPDX-License-Identifier: Apache-2.0

module constituf_chaboche
    !! Three-dimensional, small-strain \(J_2\) plasticity with nonlinear
    !! kinematic and linear isotropic hardening.
    !!
    !! This single-backstress Armstrong--Frederick specialization of
    !! Chaboche hardening uses tension-positive Cauchy stress
    !! \(\boldsymbol{\sigma}\) (stress units) and infinitesimal strain
    !! \(\boldsymbol{\varepsilon}\) (dimensionless). Each symmetric tensor
    !! \(\boldsymbol{A}\) is packed as
    !! \[\operatorname{mandel}(\boldsymbol{A})=
    !! [A_{11},A_{22},A_{33},\sqrt{2}A_{23},
    !! \sqrt{2}A_{13},\sqrt{2}A_{12}]^{\mathsf T}.\]
    !! The committed internal state comprises deviatoric plastic strain
    !! \(\boldsymbol{\varepsilon}^p_n\), deviatoric backstress
    !! \(\boldsymbol{X}_n\) (stress units) and accumulated plastic multiplier
    !! \(\alpha_n\ge0\).
    !!
    !! Let \(\boldsymbol{s}=\mathrm{dev}(\boldsymbol{\sigma})\) and
    !! \(q=\sqrt{3/2}\|\boldsymbol{s}-\boldsymbol{X}\|\).
    !! The yield condition is
    !! \[f=q-(\sigma_y+H\alpha)\le0,\]
    !! where \(\sigma_y>0\) is the initial yield stress (stress units) and
    !! \(H\ge0\) is the linear isotropic hardening modulus (stress units).
    !! The Armstrong--Frederick evolution law is
    !! \[\dot{\boldsymbol X}=(2/3)C\dot{\boldsymbol\varepsilon}^p
    !! -\gamma\boldsymbol X\dot\alpha.\]
    !! For a plastic increment \(r=\Delta\gamma\ge0\),
    !! the backward-Euler discrete updates are
    !! \[\boldsymbol{\varepsilon}^p_{n+1}
    !! =\boldsymbol{\varepsilon}^p_n+r\boldsymbol{n},\qquad
    !! \alpha_{n+1}=\alpha_n+r,\]
    !! \[\boldsymbol{X}_{n+1}
    !! =\frac{\boldsymbol{X}_n+(2/3)Cr\boldsymbol{n}}
    !! {1+\gamma r},\qquad
    !! \boldsymbol{s}_{n+1}=\boldsymbol{s}_{\mathrm{tr}}
    !! -2Gr\boldsymbol{n},\]
    !! \[\boldsymbol{n}
    !! =\frac{3}{2}\frac{\boldsymbol{s}_{n+1}
    !! -\boldsymbol{X}_{n+1}}{q_{n+1}}.\]
    !! Here \(K>0\), \(G>0\) (stress units) are the elastic bulk and
    !! shear moduli, \(C\ge0\) (stress units) is the kinematic modulus and
    !! \(\gamma\ge0\) is the dimensionless dynamic-recovery coefficient.
    !! The elastic trial deviator is
    !! \(\boldsymbol{s}_{\mathrm{tr}}=
    !! 2G\,\mathrm{dev}(\boldsymbol{\varepsilon}
    !! -\boldsymbol{\varepsilon}^p_n)\).
    !!
    !! A safeguarded scalar Newton solve uses
    !! \[R(r)=\sqrt{3/2}
    !! \left\|\boldsymbol{s}_{\mathrm{tr}}
    !! -\frac{\boldsymbol{X}_n}{1+\gamma r}\right\|
    !! -\left(3G+\frac{C}{1+\gamma r}+H\right)r
    !! -\sigma_y-H\alpha_n=0.\]
    !! The algorithmic tangent differentiates the discrete return map,
    !! including the changes of \(r\), \(\boldsymbol{X}\) and
    !! \(\boldsymbol{n}\). Analytical parameter sensitivities use
    !! implicit differentiation of the same converged scalar residual,
    !! including all accepted backstress and plastic-strain history.
    !! Return convergence and yield-switch bands use relative stress scales.
    !! The classical tangent is undefined at yield activation. Tangent or
    !! sensitivity requests at that switch report `cf_nondifferentiable`.
    !! Only successful accepted steps commit history.
    !!
    !! The model is rate-independent, isothermal and limited to
    !! infinitesimal kinematics and a single backstress. It does not
    !! claim a multi-backstress fit, ratcheting calibration or
    !! finite-strain behavior.
    !!
    !! References: J.-L. Chaboche, *Constitutive Equations for Cyclic
    !! Plasticity and Cyclic Viscoplasticity*, International Journal of
    !! Plasticity **5** (1989), 247--302,
    !! DOI: 10.1016/0749-6419(89)90015-6. P. J. Armstrong and
    !! C. O. Frederick, CEGB Report RD/B/N731 (1966).
    !! C. Suchocki, *On determination of material parameters in cyclic
    !! plasticity*, Acta Mechanica 236 (2025), 4569--4602,
    !! DOI: 10.1007/s00707-025-04395-6. Section 2.1, Eqs. (3)--(5),
    !! and Appendix B, Eqs. (66)--(72), give the yield, flow and
    !! backward-Euler backstress equations. This implementation uses one
    !! backstress and the linear isotropic law \(R(\alpha)=H\alpha\).
    !! The numerical integration and tangent are independently implemented.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument, cf_nondifferentiable
    implicit none (type, external)
    private

    integer, parameter :: nmandel=6
    integer, parameter :: nstate=13
    public :: chaboche_t

    type, extends(material_model_t) :: chaboche_t
        !! Construct as `chaboche_t(K,G,sy,H,C,gamma)` with moduli in the chosen stress unit.
        real(wp) :: bulk_modulus=0.0_wp
            !! \(K>0\), bulk modulus (stress units).
        real(wp) :: shear_modulus=0.0_wp
            !! \(G>0\), shear modulus (stress units).
        real(wp) :: yield_stress=0.0_wp
            !! \(\sigma_y>0\), initial yield stress (stress units).
        real(wp) :: hardening_modulus=0.0_wp
            !! \(H\geq0\), isotropic linear hardening modulus (stress units).
        real(wp) :: kinematic_modulus=0.0_wp
            !! \(C\geq0\), Armstrong--Frederick kinematic modulus (stress units).
        real(wp) :: recovery=0.0_wp
            !! \(\gamma\geq0\), dynamic-recovery coefficient (dimensionless).
    contains
        procedure :: describe => describe_chaboche
            !! Describe the material input, response, state and scalar parameter schema.
        procedure :: eval_batch => evaluate_chaboche
            !! Evaluate the batched material response and trial history.
        procedure :: initialize_sensitivity => initialize_chaboche_sensitivity
            !! Initialize analytical derivatives of the reference real history.
        procedure :: eval_sensitivity_batch => sensitivity_chaboche
            !! Evaluate exact selected-parameter sensitivities through the discrete material update.
    end type chaboche_t

contains

    subroutine describe_chaboche(self,spec,status)
        !! Validate model parameters and register `strain.small` input and `stress.cauchy` response fields.
        class(chaboche_t),intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        type(material_spec_t),intent(out) :: spec
            !! Material schema receiving the packed field descriptions and supported capabilities.
        integer,intent(out) :: status
            !! Status code; `cf_ok` indicates success.

        status=cf_invalid_argument
        if (.not. all(ieee_is_finite([self%bulk_modulus,self%shear_modulus, self%yield_stress,self%hardening_modulus,self%kinematic_modulus,self%recovery]))) return
        if (self%bulk_modulus<=0.0_wp.or.self%shear_modulus<=0.0_wp.or.self%yield_stress<=0.0_wp) return
        if (self%hardening_modulus<0.0_wp.or.self%kinematic_modulus<0.0_wp.or.self%recovery<0.0_wp) return
        if (.not. ieee_is_finite(3.0_wp*self%shear_modulus+ self%kinematic_modulus+self%hardening_modulus)) return

        call spec%add_input("strain.small", nmandel, "mandel", "small_strain", "spatial", [3,3])
        call spec%add_output("stress.cauchy", nmandel, "mandel", "cauchy_stress", "spatial", [3,3])
        call spec%add_state("plastic.strain", nmandel, "mandel", "plastic_strain", "spatial", [3,3])
        call spec%add_state("plastic.backstress", nmandel, "mandel", "backstress", "spatial", [3,3])
        call spec%add_state("plastic.accumulated", 1, "scalar", "equivalent_plastic_strain", "material")
        call spec%add_parameter("bulk_modulus")
        call spec%add_parameter("shear_modulus")
        call spec%add_parameter("yield_stress")
        call spec%add_parameter("hardening_modulus")
        call spec%add_parameter("kinematic_modulus")
        call spec%add_parameter("recovery")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status=cf_ok
    end subroutine describe_chaboche

    subroutine evaluate_chaboche(self,x_old,x_new,state_old,state_trial,response,dt,tangent,point_status)
        !! All trial evaluations start from committed plastic history.
        class(chaboche_t),intent(in) :: self
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
        real(wp) :: eps_dev(nmandel)
        real(wp) :: s_tr(nmandel)
        real(wp) :: back(nmandel)
        real(wp) :: u(nmandel)
        real(wp) :: direction(nmandel)
        real(wp) :: dur(nmandel)
        real(wp) :: dstrain(nmandel)
        real(wp) :: du(nmandel)
        real(wp) :: dn(nmandel)
        real(wp) :: column(nmandel)
        real(wp) :: qtrial
        real(wp) :: ftrial
        real(wp) :: yield_limit
        real(wp) :: r
        real(wp) :: lower
        real(wp) :: upper
        real(wp) :: rmid
        real(wp) :: gr
        real(wp) :: gprime
        real(wp) :: ghi
        real(wp) :: q
        real(wp) :: fac
        real(wp) :: denom
        real(wp) :: dq
        real(wp) :: dr
        real(wp) :: stress_scale
        real(wp) :: trace_strain
        integer :: ipoint
        integer :: i
        integer :: j
        integer :: iter
        logical :: converged
        logical :: plastic

        identity=[1.0_wp,1.0_wp,1.0_wp,0.0_wp,0.0_wp,0.0_wp]
        projector=0.0_wp
        do i=1,nmandel
            projector(i,i)=1.0_wp
        end do
        projector(1:3,1:3)=projector(1:3,1:3)-1.0_wp/3.0_wp
        state_trial=state_old

        point_loop: do ipoint=1,size(x_new,2)
            point_status(ipoint)=cf_invalid_argument
            if (.not. all(ieee_is_finite(x_new(:,ipoint)))) cycle point_loop
            if (.not. all(ieee_is_finite(state_old(:,ipoint)))) cycle point_loop
            if (state_old(nstate,ipoint)<0.0_wp) cycle point_loop
            trace_strain=sum(x_new(1:3,ipoint))
            eps_dev=x_new(:,ipoint)-trace_strain*identity/3.0_wp
            back=state_old(7:12,ipoint)
            s_tr=2.0_wp*self%shear_modulus*(eps_dev-state_old(1:6,ipoint))
            u=s_tr-back
            qtrial=sqrt(1.5_wp*dot_product(u,u))
            yield_limit=self%yield_stress+self%hardening_modulus*state_old(nstate,ipoint)
            ftrial=qtrial-yield_limit
            if (.not. all(ieee_is_finite([qtrial,yield_limit,ftrial]))) cycle point_loop
            stress_scale=max(qtrial,yield_limit)
            plastic=ftrial>0.0_wp
            r=0.0_wp
            direction=0.0_wp

            if (plastic) then
                lower=0.0_wp
                upper=ftrial/(3.0_wp*self%shear_modulus+self%hardening_modulus)
                converged=.false.
                bracket_loop: do iter=1,64
                    call residual(self,s_tr,back,yield_limit,upper,ghi,gprime)
                    if (.not. ieee_is_finite(ghi)) exit bracket_loop
                    if (ghi<=0.0_wp) then
                        converged=.true.
                        exit bracket_loop
                    end if
                    upper=2.0_wp*upper
                    if (.not. ieee_is_finite(upper)) exit bracket_loop
                end do bracket_loop
                if (.not. converged) cycle point_loop

                r=0.5_wp*(lower+upper)
                converged=.false.
                solve_loop: do iter=1,64
                    call residual(self,s_tr,back,yield_limit,r,gr,gprime)
                    if (.not. all(ieee_is_finite([gr,gprime]))) exit solve_loop
                    if (abs(gr)/stress_scale<=32.0_wp*epsilon(1.0_wp)) then
                        converged=.true.
                        exit solve_loop
                    end if
                    if (gr>0.0_wp) then
                        lower=r
                    else
                        upper=r
                    end if
                    rmid=0.5_wp*(lower+upper)
                    if (gprime<0.0_wp) then
                        denom=r-gr/gprime
                        if (ieee_is_finite(denom).and.denom>lower.and.denom<=upper) rmid=denom
                    end if
                    r=rmid
                end do solve_loop
                if (.not. converged) cycle point_loop

                fac=1.0_wp+self%recovery*r
                u=s_tr-back/fac
                q=sqrt(1.5_wp*dot_product(u,u))
                if (.not. ieee_is_finite(q)) cycle point_loop
                if (q<=0.0_wp) cycle point_loop
                direction=1.5_wp*u/q
                state_trial(1:6,ipoint)=state_old(1:6,ipoint)+r*direction
                state_trial(7:12,ipoint)=(back+(2.0_wp/3.0_wp)* self%kinematic_modulus*r*direction)/fac
                state_trial(nstate,ipoint)=state_old(nstate,ipoint)+r
                if (.not. all(ieee_is_finite(state_trial(:,ipoint)))) cycle point_loop
            end if

            response(:,ipoint)=s_tr-2.0_wp*self%shear_modulus*r*direction+ self%bulk_modulus*trace_strain*identity
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop
            if (present(tangent)) then
                if (abs(ftrial)/stress_scale<=64.0_wp*epsilon(1.0_wp)) then
                    point_status(ipoint)=cf_nondifferentiable
                    cycle point_loop
                end if
                if (plastic) then
                    dur=self%recovery*back/fac**2
                    gprime=dot_product(direction,dur)-3.0_wp*self%shear_modulus- self%kinematic_modulus/fac+self%kinematic_modulus*self%recovery*r/fac**2- self%hardening_modulus
                    if (abs(gprime)<=tiny(gprime)) cycle point_loop
                end if
                do j=1,nmandel
                    dstrain=2.0_wp*self%shear_modulus*projector(:,j)
                    column=dstrain
                    if (j<=3) column=column+self%bulk_modulus*identity
                    if (plastic) then
                        dr=-dot_product(direction,dstrain)/gprime
                        du=dstrain+dur*dr
                        dq=dot_product(direction,du)
                        dn=(1.5_wp/q)*(du-u*dq/q)
                        column=column-2.0_wp*self%shear_modulus*(dr*direction+r*dn)
                    end if
                    tangent(:,j,ipoint)=column
                end do
                if (.not. all(ieee_is_finite(tangent(:,:,ipoint)))) cycle point_loop
            end if
            point_status(ipoint)=cf_ok
        end do point_loop
    end subroutine evaluate_chaboche

    pure subroutine residual(self,s_tr,back,yield_limit,r,g,gprime)
        !! Nonlinear scalar consistency and analytical derivative.
        class(chaboche_t),intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        real(wp),intent(in) :: s_tr(nmandel)
            !! Elastic trial deviatoric stress in Mandel components.
        real(wp),intent(in) :: back(nmandel)
            !! Committed kinematic backstress in the six-component Mandel convention.
        real(wp),intent(in) :: yield_limit
            !! Current yield threshold including the accepted isotropic hardening contribution.
        real(wp),intent(in) :: r
            !! Trial scalar plastic multiplier increment in the local return equation.
        real(wp),intent(out) :: g
            !! Scalar Chaboche consistency residual at the trial multiplier.
        real(wp),intent(out) :: gprime
            !! Derivative of the Chaboche consistency residual with respect to `r`.
        real(wp) :: u(nmandel)
        real(wp) :: fac
        real(wp) :: q
        real(wp) :: dq

        fac=1.0_wp+self%recovery*r
        u=s_tr-back/fac
        q=sqrt(1.5_wp*dot_product(u,u))
        dq=0.0_wp
        if (q>0.0_wp) dq=dot_product(1.5_wp*u/q,self%recovery*back/fac**2)
        g=q-(3.0_wp*self%shear_modulus+self%kinematic_modulus/fac+ self%hardening_modulus)*r-yield_limit
        gprime=dq-3.0_wp*self%shear_modulus-self%kinematic_modulus/fac+ self%kinematic_modulus*self%recovery*r/fac**2-self%hardening_modulus
    end subroutine residual


    subroutine initialize_chaboche_sensitivity(self,initial_input,initial_state,selected_parameters,sensitivity,status)
        !! Initialize analytical history sensitivities for the virgin material state.
        class(chaboche_t),intent(in) :: self
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
    end subroutine initialize_chaboche_sensitivity

    subroutine sensitivity_chaboche(self,x_old,x_new,state_old,state_trial,selected_parameters, sensitivity_old,sensitivity_trial,response,dsdp,dt,tangent,point_status)
        !! Differentiate the converged \(R(r)=0\) by
        !! \(\mathrm{d}r/\mathrm{d}p=-R_p/R_r\), with
        !! \(R_p\) including the complete committed backstress/history path.
        class(chaboche_t),intent(in) :: self
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
        real(wp) :: strain_trial(nmandel)
        real(wp) :: back(nmandel)
        real(wp) :: u(nmandel)
        real(wp) :: direction(nmandel)
        real(wp) :: dstrain(nmandel)
        real(wp) :: dback(nmandel)
        real(wp) :: du_fixed(nmandel)
        real(wp) :: du(nmandel)
        real(wp) :: ddirection(nmandel)
        real(wp) :: d_bulk
        real(wp) :: d_shear
        real(wp) :: d_yield
        real(wp) :: d_hardening
        real(wp) :: d_kinematic
        real(wp) :: d_recovery
        real(wp) :: r
        real(wp) :: fac
        real(wp) :: qtrial
        real(wp) :: ftrial
        real(wp) :: yield_limit
        real(wp) :: q
        real(wp) :: gprime
        real(wp) :: dr
        real(wp) :: dfac
        real(wp) :: dfixed
        real(wp) :: dbound
        real(wp) :: trace_strain
        real(wp) :: derivative_q
        real(wp) :: derivative_fac
        real(wp) :: derivative_increment
        real(wp) :: new_back(nmandel)
        real(wp) :: stress_derivative(nmandel)
        real(wp) :: stress_scale
        integer :: ipoint
        integer :: iparam
        logical :: plastic

        call evaluate_chaboche(self,x_old,x_new,state_old,state_trial,response,dt,tangent,point_status)
        identity=[1.0_wp,1.0_wp,1.0_wp,0.0_wp,0.0_wp,0.0_wp]
        sensitivity_trial=sensitivity_old
        dsdp=0.0_wp
        point_loop: do ipoint=1,size(x_new,2)
            if (point_status(ipoint)/=cf_ok) cycle point_loop
            trace_strain=sum(x_new(1:3,ipoint))
            eps_dev=x_new(:,ipoint)-trace_strain*identity/3.0_wp
            strain_trial=2.0_wp*self%shear_modulus*(eps_dev-state_old(1:6,ipoint))
            back=state_old(7:12,ipoint)
            u=strain_trial-back
            qtrial=sqrt(1.5_wp*dot_product(u,u))
            yield_limit=self%yield_stress+self%hardening_modulus*state_old(nstate,ipoint)
            ftrial=qtrial-yield_limit
            stress_scale=max(qtrial,yield_limit)
            if (abs(ftrial)/stress_scale<=64.0_wp*epsilon(1.0_wp)) then
                point_status(ipoint)=cf_nondifferentiable
                cycle point_loop
            end if
            plastic=ftrial>0.0_wp
            r=state_trial(nstate,ipoint)-state_old(nstate,ipoint)
            fac=1.0_wp+self%recovery*r
            q=0.0_wp
            direction=0.0_wp
            gprime=0.0_wp
            new_back=state_trial(7:12,ipoint)
            if (plastic) then
                u=strain_trial-back/fac
                q=sqrt(1.5_wp*dot_product(u,u))
                if (q<=0.0_wp) then
                    point_status(ipoint)=cf_invalid_argument
                    cycle point_loop
                end if
                direction=1.5_wp*u/q
                gprime=dot_product(direction,self%recovery*back/fac**2) - 3.0_wp*self%shear_modulus-self%kinematic_modulus/fac+ self%kinematic_modulus*self%recovery*r/fac**2-self%hardening_modulus
                if (abs(gprime)<=tiny(1.0_wp)) then
                    point_status(ipoint)=cf_invalid_argument
                    cycle point_loop
                end if
            end if

            parameter_loop: do iparam=1,size(selected_parameters)
                d_bulk=0.0_wp
                d_shear=0.0_wp
                d_yield=0.0_wp
                d_hardening=0.0_wp
                d_kinematic=0.0_wp
                d_recovery=0.0_wp
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
                    d_kinematic=1.0_wp
                case(6)
                    d_recovery=1.0_wp
                case default
                    point_status(ipoint)=cf_invalid_argument
                    exit parameter_loop
                end select
                dstrain=2.0_wp*d_shear*(eps_dev-state_old(1:nmandel,ipoint))- 2.0_wp*self%shear_modulus*sensitivity_old(1:nmandel,iparam,ipoint)
                stress_derivative=dstrain+d_bulk*trace_strain*identity
                if (plastic) then
                    dback=sensitivity_old(7:12,iparam,ipoint)
                    du_fixed=dstrain-dback/fac+back*(d_recovery*r)/fac**2
                    derivative_q=dot_product(direction,du_fixed)
                    derivative_fac=3.0_wp*d_shear+d_kinematic/fac- self%kinematic_modulus*d_recovery*r/fac**2+d_hardening
                    dbound=d_yield+d_hardening*state_old(nstate,ipoint)+ self%hardening_modulus*sensitivity_old(nstate,iparam,ipoint)
                    dfixed=derivative_q-r*derivative_fac-dbound
                    dr=-dfixed/gprime
                    dfac=d_recovery*r+self%recovery*dr
                    du=dstrain-dback/fac+back*dfac/fac**2
                    derivative_increment=dot_product(direction,du)
                    ddirection=(1.5_wp/q)*(du-u*derivative_increment/q)
                    sensitivity_trial(1:nmandel,iparam,ipoint)= sensitivity_old(1:nmandel,iparam,ipoint)+dr*direction+r*ddirection
                    sensitivity_trial(7:12,iparam,ipoint)= (dback+(2.0_wp/3.0_wp)*(d_kinematic*r*direction+ self%kinematic_modulus*(dr*direction+r*ddirection)))/fac- new_back*dfac/fac
                    sensitivity_trial(nstate,iparam,ipoint)=sensitivity_old(nstate,iparam,ipoint)+dr
                    stress_derivative=stress_derivative-2.0_wp*d_shear*r*direction- 2.0_wp*self%shear_modulus*(dr*direction+r*ddirection)
                end if
                dsdp(:,iparam,ipoint)=stress_derivative
            end do parameter_loop
            if (point_status(ipoint)/=cf_ok) cycle point_loop
            if (.not. all(ieee_is_finite(dsdp(:,:,ipoint))) .or. .not. all(ieee_is_finite(sensitivity_trial(:,:,ipoint)))) then
                point_status(ipoint)=cf_invalid_argument
            end if
        end do point_loop
    end subroutine sensitivity_chaboche

end module constituf_chaboche

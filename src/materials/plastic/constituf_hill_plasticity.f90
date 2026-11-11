! SPDX-License-Identifier: Apache-2.0

module constituf_hill_plasticity
    !! Three-dimensional associative Hill-1948 anisotropic plasticity.
    !!
    !! An isothermal, small-strain solid with 3D isotropic elastic moduli
    !! \(K>0\) and \(G>0\) uses the Hill quadratic stress metric
    !! \[q^2=F(\sigma_{22}-\sigma_{33})^2+
    !! G_H(\sigma_{33}-\sigma_{11})^2+
    !! H_H(\sigma_{11}-\sigma_{22})^2+
    !! 2L\sigma_{23}^2+2M\sigma_{13}^2+2N\sigma_{12}^2
    !! =\boldsymbol\sigma^T\boldsymbol Q\boldsymbol\sigma.\]
    !! The implemented positive-coefficient subclass uses dimensionless
    !! Hill coefficients \((F,G_H,H_H,L,M,N)\). At virgin yield these equal
    !! \(Y_0^2\) times Hill's original inverse-stress-squared coefficients.
    !! The Mandel shear entries of \(\boldsymbol{Q}\) are \((L,M,N)\).
    !! Hill \(G_H,H_H\) are *not* the elastic shear or hardening moduli.
    !! The isotropic von Mises limit is \((1/2,1/2,1/2,3/2,3/2,3/2)\).
    !! Input \(\boldsymbol{\varepsilon}\) and output Cauchy stress are symmetric
    !! spatial tensors. For either tensor \(\boldsymbol{A}\), the Mandel vector is
    !! \[\operatorname{mandel}(\boldsymbol{A})=
    !! [A_{11},A_{22},A_{33},\sqrt{2}A_{23},
    !! \sqrt{2}A_{13},\sqrt{2}A_{12}]^{\mathsf T}.\]
    !! Tensile components are positive. Strain is dimensionless.
    !! Stress and moduli are in the chosen stress unit.
    !!
    !! With accumulated plastic multiplier \(\alpha\) and initial positive
    !! yield stress \(Y_0\), linear isotropic hardening \(h\geq0\) defines
    !! \[f=q-(Y_0+h\alpha)\leq0,\qquad
    !! \dot{\boldsymbol\varepsilon}^{p}
    !! =\dot\gamma\,\boldsymbol n,\quad
    !! \boldsymbol n=\partial_{\boldsymbol\sigma}q
    !! =\boldsymbol Q\boldsymbol\sigma/q,\quad
    !! \dot\alpha=\dot\gamma\geq0.\]
    !! The plastic multiplier is defined via the Hill work-conjugate
    !! equivalent measure. It equals the \(J_2\) equivalent plastic strain
    !! increment in the isotropic limit.
    !!
    !! The backward-Euler return from the committed plastic strain and
    !! hardening is solved by a bracketed scalar root in \(t=\Delta\gamma/q\):
    !! \[\boldsymbol\sigma(t)=
    !! (\boldsymbol I+2Gt\boldsymbol Q)^{-1}
    !! \boldsymbol\sigma_\mathrm{trial},\qquad
    !! q(t)(1-ht)=Y_0+h\alpha_n.\]
    !! Plastic strain advances by \(\Delta\gamma\boldsymbol{n}\), \(\alpha\) by \(\Delta\gamma\).
    !! Hydrostatic stress remains elastic because \(\boldsymbol{Q}\) annihilates the Mandel identity vector \(\boldsymbol{I}\).
    !! The elastic pressure is retained separately. Only deviatoric stress
    !! enters the numerical return solve and Hill metric products.
    !! Positive Hill coefficients make \(\boldsymbol{Q}\) positive on deviatoric stress.
    !! The monotone scalar residual admits a unique return for \(Y_0>0\).
    !! The algorithm uses pivoted, fixed-size linear solves and bisection
    !! to a kind-scaled residual tolerance. Failure returns a point status.
    !!
    !! The *consistent discrete* tangent, at fixed committed history,
    !! follows the exact implicit linearization. Define
    !! \[\boldsymbol H_n=(\boldsymbol Q-
    !! \boldsymbol n\otimes\boldsymbol n)/q,\qquad
    !! \boldsymbol A=\boldsymbol I+
    !! 2G\Delta\gamma\boldsymbol H_n,\quad
    !! \boldsymbol B=\boldsymbol A^{-1}\boldsymbol C_e.\]
    !! Then
    !! \[\boldsymbol C_\mathrm{alg}
    !! =\boldsymbol B-
    !! \frac{(\boldsymbol B\boldsymbol n)\otimes
    !! (\boldsymbol n^T\boldsymbol B)}
    !! {h+\boldsymbol n^T\boldsymbol B\boldsymbol n}.\]
    !! This is *not* a finite-difference tangent. At a yield activation
    !! switch the classical tangent is not unique. Use off-switch trials
    !! to verify it.
    !!
    !! Optional analytical parameter sensitivities differentiate the exact
    !! *discrete* return through all accepted plastic histories while the
    !! prescribed strain history is fixed. The selected scalar parameters
    !! are `bulk_modulus`, `shear_modulus`, `yield_stress`,
    !! `hardening_modulus`, followed by `hill_f`, `hill_g`, `hill_h`,
    !! `hill_l`, `hill_m`, `hill_n`. The latter six are dimensionless.
    !! For a parameter \(p\), write the converged plastic constraint as
    !! \[\boldsymbol R=\boldsymbol\sigma-\boldsymbol\sigma_\mathrm{trial}
    !! +\Delta\gamma\,\boldsymbol C_e\boldsymbol n=\boldsymbol0,\quad
    !! g=q-Y_0-h(\alpha_n+\Delta\gamma)=0.\]
    !! Let \(\mathbb A=\boldsymbol I+
    !! \Delta\gamma\,\boldsymbol C_e\mathbb H_n\), where
    !! \(\mathbb H_n=(\boldsymbol Q-\boldsymbol n\otimes\boldsymbol n)/q\).
    !! Define the explicit Hill metric derivative
    !! \[q_p^{(0)}=\boldsymbol\sigma^T\boldsymbol Q_p
    !! \boldsymbol\sigma/(2q),\quad
    !! \boldsymbol n_p^{(0)}=
    !! \boldsymbol Q_p\boldsymbol\sigma/q
    !! -\boldsymbol n\,q_p^{(0)}/q.\]
    !! With \(D_p\boldsymbol\varepsilon_n^p\) and \(D_p\alpha_n\)
    !! carried from the last committed step, define
    !! \[\boldsymbol v=\boldsymbol C_{e,p}
    !! (\boldsymbol\varepsilon-\boldsymbol\varepsilon_n^p)
    !! -\boldsymbol C_e D_p\boldsymbol\varepsilon_n^p
    !! -\Delta\gamma(\boldsymbol C_{e,p}\boldsymbol n+
    !! \boldsymbol C_e\boldsymbol n_p^{(0)}),\quad
    !! \boldsymbol b=\mathbb A^{-1}\boldsymbol v,\quad
    !! \boldsymbol w=\mathbb A^{-1}\boldsymbol C_e\boldsymbol n.\]
    !! The **total discrete path derivative** is
    !! \[D_p\Delta\gamma=
    !! \frac{\boldsymbol n\cdot\boldsymbol b+q_p^{(0)}
    !! -Y_{0,p}-(\alpha_n+\Delta\gamma)h_p-hD_p\alpha_n}
    !! {h+\boldsymbol n\cdot\boldsymbol w},\quad
    !! D_p\boldsymbol\sigma=\boldsymbol b-\boldsymbol w
    !! D_p\Delta\gamma.\]
    !! The plastic-strain derivative advances by
    !! \[D_p\boldsymbol\varepsilon^p_{n+1}=
    !! D_p\boldsymbol\varepsilon^p_n+
    !! (D_p\Delta\gamma)\boldsymbol n+
    !! \Delta\gamma(\mathbb H_nD_p\boldsymbol\sigma+
    !! \boldsymbol n_p^{(0)}),\]
    !! and \(D_p\alpha_{n+1}=D_p\alpha_n+D_p\Delta\gamma\).
    !! For elastic steps, the existing state derivatives are retained and
    !! \(D_p\boldsymbol\sigma=\boldsymbol C_{e,p}
    !! (\boldsymbol\varepsilon-\boldsymbol\varepsilon^p_n)-
    !! \boldsymbol C_eD_p\boldsymbol\varepsilon^p_n\).
    !! At active-set switches these derivatives need not be unique.
    !! Sensitivity-enabled trials report `cf_nondifferentiable` in a
    !! kind-scaled yield-surface band. The 3D strain trajectory and
    !! increments are held fixed, so no coupled structural sensitivity
    !! is claimed.
    !!
    !! Limitations: fixed material axes aligned with Cartesian basis.
    !! No rotations of anisotropy axes, finite strain, kinematic hardening,
    !! rate dependence, or plane-stress constraint solve.
    !!
    !! Reference: R. Hill, "A theory of the yielding and plastic flow of
    !! anisotropic metals", Proc. R. Soc. Lond. A 193 (1948), 281--297,
    !! DOI: 10.1098/rspa.1948.0045. Section 3, Eq. (5), gives the quadratic
    !! yield surface. Section 4, Eq. (9), gives its associated flow.
    !! The fixed-ratio linear hardening and backward-Euler integration
    !! are independently derived extensions. J. C. Simo and T. J. R. Hughes,
    !! Computational Inelasticity, Springer (1998), DOI: 10.1007/b98904.
    !! Code provenance: independent implementation from those equations.
    !! No third-party source code has been translated or copied.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument, cf_evaluation_failed, cf_nondifferentiable
    implicit none (type, external)
    private

    integer, parameter :: nmandel = 6
    integer, parameter :: nstate = 7
    integer, parameter :: max_root_iterations = 200
    public :: hill_plasticity_t

    type, extends(material_model_t) :: hill_plasticity_t
        !! Hill-1948 associated plastic solid with linear isotropic hardening.
        !!
        !! Construct with `hill_plasticity_t(bulk_modulus, shear_modulus, yield_stress, hardening_modulus, hill_coefficients)`.
        !! The six entries are \((F,G_H,H_H,L,M,N)\).
        !! All parameters must be finite and positive, except \(h\geq0\).
        !! Default coefficients give the isotropic \(J_2\) surface.
        real(wp) :: bulk_modulus = 0.0_wp
            !! Bulk modulus \(K>0\) in the chosen stress unit.
        real(wp) :: shear_modulus = 0.0_wp
            !! Elastic shear modulus \(G>0\) in the chosen stress unit.
        real(wp) :: yield_stress = 0.0_wp
            !! Reference positive yield stress \(Y_0\) in the chosen stress unit.
        real(wp) :: hardening_modulus = 0.0_wp
            !! Linear isotropic hardening \(h\geq0\) in the chosen stress unit.
        real(wp) :: hill_coefficients(nmandel) = [0.5_wp,0.5_wp,0.5_wp,1.5_wp,1.5_wp,1.5_wp]
            !! Dimensionless Hill coefficients \((F,G_H,H_H,L,M,N)\).
    contains
        procedure :: describe => describe_hill
            !! Describe the material input, response, state and scalar parameter schema.
        procedure :: eval_batch => evaluate_hill
            !! Evaluate the batched material response and trial history.
        procedure :: initialize_sensitivity => initialize_hill_sensitivity
            !! Initialize analytical derivatives of the reference real history.
        procedure :: eval_sensitivity_batch => sensitivity_hill
            !! Evaluate exact selected-parameter sensitivities through the discrete material update.
    end type hill_plasticity_t

contains

    subroutine describe_hill(self, spec, status)
        !! Validate model parameters and register `strain.small` input and `stress.cauchy` response fields.
        class(hill_plasticity_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema receiving the packed field descriptions and supported capabilities.
        integer, intent(out) :: status
            !! Status code; `cf_ok` indicates success.

        status = cf_invalid_argument
        if (.not. ieee_is_finite(self%bulk_modulus) .or. .not. ieee_is_finite(self%shear_modulus)) return
        if (.not. ieee_is_finite(self%yield_stress) .or. .not. ieee_is_finite(self%hardening_modulus)) return
        if (.not. all(ieee_is_finite(self%hill_coefficients))) return
        if (self%bulk_modulus <= 0.0_wp .or. self%shear_modulus <= 0.0_wp) return
        if (self%yield_stress <= 0.0_wp .or. self%hardening_modulus < 0.0_wp) return
        if (any(self%hill_coefficients <= 0.0_wp)) return
        if (.not. ieee_is_finite(self%bulk_modulus + 4.0_wp*self%shear_modulus/3.0_wp)) return

        call spec%add_input("strain.small", nmandel, "mandel", "small_strain", "spatial", [3,3])
        call spec%add_output("stress.cauchy", nmandel, "mandel", "cauchy_stress", "spatial", [3,3])
        call spec%add_state("plastic.strain", nmandel, "mandel", "plastic_strain", "spatial", [3,3])
        call spec%add_state("plastic.accumulated", 1, "scalar", "equivalent_plastic_strain", "material")
        call spec%add_parameter("bulk_modulus")
        call spec%add_parameter("shear_modulus")
        call spec%add_parameter("yield_stress")
        call spec%add_parameter("hardening_modulus")
        call spec%add_parameter("hill_f")
        call spec%add_parameter("hill_g")
        call spec%add_parameter("hill_h")
        call spec%add_parameter("hill_l")
        call spec%add_parameter("hill_m")
        call spec%add_parameter("hill_n")
        call spec%set_has_sensitivity(.true.)
        call spec%set_has_tangent(.true.)
        status = cf_ok
    end subroutine describe_hill

    subroutine initialize_hill_sensitivity(self, initial_input, initial_state, selected_parameters, sensitivity, status)
        !! A virgin zero-plastic-strain state is parameter-independent.
        class(hill_plasticity_t), intent(in) :: self
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
    end subroutine initialize_hill_sensitivity

    subroutine sensitivity_hill(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Fused forward differentiation of the accepted Hill return map.
        !!
        !! Response/state updates reuse the *same* base batch kernel and
        !! local return solution. Linearization uses analytic Jacobians
        !! and never perturbs parameters in production material evaluation.
        !! dsdp(6,nselected,npoints) and sensitivity_trial(7,nselected,npoints)
        !! preserve selected parameter order and atomic trial semantics.
        class(hill_plasticity_t), intent(in) :: self
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
        real(wp), intent(out) :: response(:,:)
            !! Packed `stress.cauchy` responses `(noutput,npoints)`.
        real(wp), intent(out) :: sensitivity_trial(:,:,:)
            !! Updated real-history parameter derivatives `(nstate,nselected,npoints)`.
        real(wp), intent(out) :: dsdp(:,:,:)
            !! Exact total response derivatives `(noutput,nselected,npoints)` at fixed prescribed driving histories.
        real(wp), intent(in) :: dt
            !! Prescribed time increment; a positive value is required for time-dependent laws.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic tangent `(noutput,ninput,npoints)` at fixed committed history.
        integer, intent(inout) :: point_status(:)
            !! Per-point outcome; write `cf_ok` only for a successful constitutive update.
        real(wp) :: qmat(nmandel,nmandel)
        real(wp) :: cel(nmandel,nmandel)
        real(wp) :: a(nmandel,nmandel)
        real(wp) :: inverse(nmandel,nmandel)
        real(wp) :: identity(nmandel,nmandel)
        real(wp) :: hessian(nmandel,nmandel)
        real(wp) :: trial(nmandel)
        real(wp) :: sigma(nmandel)
        real(wp) :: flow(nmandel)
        real(wp) :: bflow(nmandel)
        real(wp) :: dcel(nmandel,nmandel)
        real(wp) :: dmetric(nmandel,nmandel)
        real(wp) :: coefficient(nmandel)
        real(wp) :: explicit_normal(nmandel)
        real(wp) :: driving(nmandel)
        real(wp) :: response_derivative(nmandel)
        real(wp) :: delta_state(nmandel)
        real(wp) :: explicit_metric(nmandel)
        real(wp) :: equivalent_trial
        real(wp) :: yield_level
        real(wp) :: surface
        real(wp) :: equivalent
        real(wp) :: gamma
        real(wp) :: alpha
        real(wp) :: param_q
        real(wp) :: param_yield
        real(wp) :: param_hardening
        real(wp) :: dgamma
        real(wp) :: denominator
        real(wp) :: lambda
        integer :: ipoint
        integer :: iparam
        integer :: i
        integer :: j
        integer :: parameter_id
        logical :: plastic_branch
        logical :: ok

        call evaluate_hill(self,x_old,x_new,state_old,state_trial,response,dt,tangent,point_status)

        call hill_metric(self%hill_coefficients,qmat)
        lambda = self%bulk_modulus - 2.0_wp*self%shear_modulus/3.0_wp
        cel = 0.0_wp
        cel(1:3,1:3) = lambda
        do i = 1, nmandel
            cel(i,i) = cel(i,i)+2.0_wp*self%shear_modulus
        end do
        identity = 0.0_wp
        do i = 1, nmandel
            identity(i,i) = 1.0_wp
        end do

        point_loop: do ipoint = 1, size(x_new,2)
            if (point_status(ipoint) /= cf_ok) cycle point_loop
            point_status(ipoint) = cf_invalid_argument
            if (.not. all(ieee_is_finite(sensitivity_old(:,:,ipoint)))) cycle point_loop

            alpha = state_old(nstate,ipoint)
            yield_level = self%yield_stress+self%hardening_modulus*alpha
            trial = x_new(:,ipoint)-state_old(1:nmandel,ipoint)
            trial(1:3) = trial(1:3)-sum(trial(1:3))/3.0_wp
            trial = 2.0_wp*self%shear_modulus*trial
            equivalent_trial = sqrt(hill_quadratic(qmat,trial))
            if (.not. ieee_is_finite(equivalent_trial)) cycle point_loop
            surface = equivalent_trial-yield_level
            if (abs(surface)/max(self%yield_stress,equivalent_trial,abs(yield_level)) <= 32.0_wp*epsilon(1.0_wp)) then
                point_status(ipoint) = cf_nondifferentiable
                cycle point_loop
            end if
            plastic_branch = surface > 0.0_wp
            sensitivity_trial(:,:,ipoint) = sensitivity_old(:,:,ipoint)
            if (plastic_branch) then
                ! Recover the elastic deviator without subtracting a large pressure from the response.
                sigma = x_new(:,ipoint)-state_trial(1:nmandel,ipoint)
                sigma(1:3) = sigma(1:3)-sum(sigma(1:3))/3.0_wp
                sigma = 2.0_wp*self%shear_modulus*sigma
                equivalent = sqrt(hill_quadratic(qmat,sigma))
                if (equivalent <= 0.0_wp .or. .not. ieee_is_finite(equivalent)) cycle point_loop
                gamma = state_trial(nstate,ipoint)-alpha
                call apply_hill_metric(qmat,sigma,flow)
                flow = flow/equivalent
                if (.not. all(ieee_is_finite(flow)) .or. gamma < 0.0_wp) cycle point_loop
                do j = 1, nmandel
                    do i = 1, nmandel
                        hessian(i,j) = (qmat(i,j)-flow(i)*flow(j))/equivalent
                    end do
                end do
                a = identity+2.0_wp*self%shear_modulus*gamma*hessian
                call solve_six(a,identity,inverse,ok)
                if (.not. ok) then
                    point_status(ipoint) = cf_evaluation_failed
                    cycle point_loop
                end if
                bflow = 2.0_wp*self%shear_modulus*matmul(inverse,flow)
                denominator = self%hardening_modulus+dot_product(flow,bflow)
                if (.not. ieee_is_finite(denominator) .or. denominator <= 0.0_wp) then
                    point_status(ipoint) = cf_evaluation_failed
                    cycle point_loop
                end if
            end if

            parameter_loop: do iparam = 1, size(selected_parameters)
                parameter_id = selected_parameters(iparam)
                dcel = 0.0_wp
                param_yield = 0.0_wp
                param_hardening = 0.0_wp
                if (parameter_id == 1) then
                    dcel(1:3,1:3) = 1.0_wp
                else if (parameter_id == 2) then
                    dcel(1:3,1:3) = -2.0_wp/3.0_wp
                    do i = 1, nmandel
                        dcel(i,i) = dcel(i,i)+2.0_wp
                    end do
                else if (parameter_id == 3) then
                    param_yield = 1.0_wp
                else if (parameter_id == 4) then
                    param_hardening = 1.0_wp
                else if (parameter_id < 5 .or. parameter_id > 10) then
                    exit parameter_loop
                end if

                delta_state = sensitivity_old(1:nmandel,iparam,ipoint)
                driving = matmul(dcel,x_new(:,ipoint)-state_old(1:nmandel,ipoint)) - matmul(cel,delta_state)
                response_derivative = driving
                if (plastic_branch) then
                    coefficient = 0.0_wp
                    if (parameter_id >= 5) coefficient(parameter_id-4) = 1.0_wp
                    call hill_metric(coefficient,dmetric)
                    call apply_hill_metric(dmetric,sigma,explicit_metric)
                    param_q = hill_quadratic(dmetric,sigma)/(2.0_wp*equivalent)
                    explicit_normal = (explicit_metric-flow*param_q)/equivalent
                    driving = driving-2.0_wp*gamma*self%shear_modulus*explicit_normal
                    if (parameter_id == 2) driving = driving-2.0_wp*gamma*flow
                    response_derivative = matmul(inverse,driving)
                    dgamma = (dot_product(flow,response_derivative)+param_q-param_yield- (alpha+gamma)*param_hardening-self%hardening_modulus*sensitivity_old(nstate,iparam,ipoint))/denominator
                    response_derivative = response_derivative-bflow*dgamma
                    sensitivity_trial(1:nmandel,iparam,ipoint) = delta_state + dgamma*flow+gamma*(matmul(hessian,response_derivative)+explicit_normal)
                    sensitivity_trial(nstate,iparam,ipoint) = sensitivity_old(nstate,iparam,ipoint)+dgamma
                end if
                dsdp(:,iparam,ipoint) = response_derivative
            end do parameter_loop
            if (iparam <= size(selected_parameters)) cycle point_loop
            if (.not. all(ieee_is_finite(dsdp(:,:,ipoint)))) cycle point_loop
            if (.not. all(ieee_is_finite(sensitivity_trial(:,:,ipoint)))) cycle point_loop
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine sensitivity_hill

    subroutine evaluate_hill(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Evaluate `stress.cauchy` and update the trial material state for the prescribed inputs.
        class(hill_plasticity_t), intent(in) :: self
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
        real(wp) :: qmat(nmandel,nmandel)
        real(wp) :: cel(nmandel,nmandel)
        real(wp) :: trial(nmandel)
        real(wp) :: sigma(nmandel)
        real(wp) :: flow(nmandel)
        real(wp) :: a(nmandel,nmandel)
        real(wp) :: b(nmandel,nmandel)
        real(wp) :: hn(nmandel,nmandel)
        real(wp) :: bb(nmandel,nmandel)
        real(wp) :: left(nmandel)
        real(wp) :: right(nmandel)
        real(wp) :: alpha
        real(wp) :: yield_level
        real(wp) :: qtrial
        real(wp) :: equivalent
        real(wp) :: low
        real(wp) :: high
        real(wp) :: t
        real(wp) :: residual
        real(wp) :: gamma
        real(wp) :: denom
        real(wp) :: stress_scale
        real(wp) :: lambda
        real(wp) :: upper_residual
        real(wp) :: trace_elastic
        real(wp) :: pressure
        real(wp) :: temp(nmandel)
        integer :: ipoint
        integer :: i
        integer :: j
        integer :: iter
        logical :: ok
        logical :: converged

        call hill_metric(self%hill_coefficients, qmat)
        lambda = self%bulk_modulus - 2.0_wp*self%shear_modulus/3.0_wp
        cel = 0.0_wp
        cel(1:3,1:3) = lambda
        do i = 1, nmandel
            cel(i,i) = cel(i,i) + 2.0_wp*self%shear_modulus
        end do
        point_status = cf_invalid_argument
        point_loop: do ipoint = 1, size(x_new,2)
            if (.not. all(ieee_is_finite(x_new(:,ipoint)))) cycle point_loop
            if (.not. all(ieee_is_finite(state_old(:,ipoint)))) cycle point_loop
            alpha = state_old(nstate,ipoint)
            if (alpha < 0.0_wp) cycle point_loop
            yield_level = self%yield_stress + self%hardening_modulus*alpha
            if (.not. ieee_is_finite(yield_level)) cycle point_loop
            temp = x_new(:,ipoint)-state_old(1:nmandel,ipoint)
            trace_elastic = sum(temp(1:3))
            pressure = self%bulk_modulus*trace_elastic
            if (.not. ieee_is_finite(pressure)) cycle point_loop
            temp(1:3) = temp(1:3)-trace_elastic/3.0_wp
            trial = 2.0_wp*self%shear_modulus*temp
            if (.not. all(ieee_is_finite(trial))) cycle point_loop
            qtrial = sqrt(hill_quadratic(qmat,trial))
            if (.not. ieee_is_finite(qtrial)) cycle point_loop

            sigma = trial
            state_trial(:,ipoint) = state_old(:,ipoint)
            if (qtrial <= yield_level) then
                response(:,ipoint) = sigma
                response(1:3,ipoint) = response(1:3,ipoint)+pressure
                if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop
                if (present(tangent)) tangent(:,:,ipoint) = cel
                point_status(ipoint) = cf_ok
                cycle point_loop
            end if

            low = 0.0_wp
            high = 1.0_wp/(2.0_wp*self%shear_modulus)
            ok = .false.
            upper_residual = huge(1.0_wp)
            bracket_loop: do iter = 1, max_root_iterations
                call stress_at_multiplier(high, self%shear_modulus, qmat, trial, sigma, equivalent, ok)
                if (.not. ok) exit bracket_loop
                upper_residual = equivalent*(1.0_wp-self%hardening_modulus*high)-yield_level
                if (.not. ieee_is_finite(upper_residual)) then
                    ok = .false.
                    exit bracket_loop
                end if
                if (upper_residual <= 0.0_wp) exit bracket_loop
                if (high >= huge(high)/2.0_wp) then
                    ok = .false.
                    exit bracket_loop
                end if
                high = 2.0_wp*high
            end do bracket_loop
            if (.not. ok .or. upper_residual > 0.0_wp) then
                point_status(ipoint) = cf_evaluation_failed
                cycle point_loop
            end if

            stress_scale = max(qtrial,yield_level)
            converged = .false.
            return_loop: do iter = 1, max_root_iterations
                t = low + 0.5_wp*(high-low)
                call stress_at_multiplier(t, self%shear_modulus, qmat, trial, sigma, equivalent, ok)
                if (.not. ok) exit return_loop
                residual = equivalent*(1.0_wp-self%hardening_modulus*t)-yield_level
                if (.not. ieee_is_finite(residual)) exit return_loop
                if (abs(residual)/stress_scale <= 64.0_wp*epsilon(1.0_wp)) then
                    converged = .true.
                    exit return_loop
                end if
                if (residual > 0.0_wp) then
                    low = t
                else
                    high = t
                end if
            end do return_loop
            if (.not. converged .or. .not. ok) then
                point_status(ipoint) = cf_evaluation_failed
                cycle point_loop
            end if
            if (equivalent <= 0.0_wp) then
                point_status(ipoint) = cf_evaluation_failed
                cycle point_loop
            end if
            gamma = t*equivalent
            call apply_hill_metric(qmat,sigma,flow)
            flow = flow/equivalent
            state_trial(1:nmandel,ipoint) = state_old(1:nmandel,ipoint) + gamma*flow
            state_trial(nstate,ipoint) = alpha + gamma
            response(:,ipoint) = sigma
            response(1:3,ipoint) = response(1:3,ipoint)+pressure
            if (.not. all(ieee_is_finite(state_trial(:,ipoint))) .or. .not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop

            if (present(tangent)) then
                do j = 1, nmandel
                    do i = 1, nmandel
                        hn(i,j) = (qmat(i,j)-flow(i)*flow(j))/equivalent
                    end do
                end do
                a = 2.0_wp*self%shear_modulus*gamma*hn
                do i = 1, nmandel
                    a(i,i) = a(i,i)+1.0_wp
                end do
                call solve_six(a, cel, bb, ok)
                if (.not. ok) then
                    point_status(ipoint) = cf_evaluation_failed
                    cycle point_loop
                end if
                left = matmul(bb,flow)
                right = matmul(flow,bb)
                denom = self%hardening_modulus + dot_product(flow,left)
                if (.not. ieee_is_finite(denom) .or. denom <= 0.0_wp) then
                    point_status(ipoint) = cf_evaluation_failed
                    cycle point_loop
                end if
                b = bb
                do j = 1, nmandel
                    do i = 1, nmandel
                        b(i,j) = bb(i,j)-left(i)*right(j)/denom
                    end do
                end do
                if (.not. all(ieee_is_finite(b))) then
                    point_status(ipoint) = cf_evaluation_failed
                    cycle point_loop
                end if
                tangent(:,:,ipoint) = b
            end if
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_hill

    pure subroutine hill_metric(coef, q)
        !! Assemble the symmetric Mandel Hill yield metric from its six anisotropy coefficients.
        real(wp), intent(in) :: coef(nmandel)
            !! Six dimensionless Hill anisotropy coefficients used to build the quadratic yield metric.
        real(wp), intent(out) :: q(nmandel,nmandel)
            !! Computed Hill quadratic yield metric.

        q = 0.0_wp
        q(1,1) = coef(2)+coef(3)
        q(2,2) = coef(1)+coef(3)
        q(3,3) = coef(1)+coef(2)
        q(1,2) = -coef(3)
        q(2,1) = q(1,2)
        q(1,3) = -coef(2)
        q(3,1) = q(1,3)
        q(2,3) = -coef(1)
        q(3,2) = q(2,3)
        q(4,4) = coef(4)
        q(5,5) = coef(5)
        q(6,6) = coef(6)
    end subroutine hill_metric

    pure function hill_quadratic(metric,stress) result(value)
        !! Evaluate the pressure-independent quadratic form without hydrostatic cancellation.
        real(wp), intent(in) :: metric(nmandel,nmandel)
        real(wp), intent(in) :: stress(nmandel)
        real(wp) :: value

        value = -metric(2,3)*(stress(2)-stress(3))**2-metric(1,3)*(stress(3)-stress(1))**2-metric(1,2)*(stress(1)-stress(2))**2+metric(4,4)*stress(4)**2+metric(5,5)*stress(5)**2+metric(6,6)*stress(6)**2
    end function hill_quadratic

    pure subroutine apply_hill_metric(metric,stress,product)
        !! Apply the Hill metric using stress differences so hydrostatic stress contributes zero.
        real(wp), intent(in) :: metric(nmandel,nmandel)
        real(wp), intent(in) :: stress(nmandel)
        real(wp), intent(out) :: product(nmandel)

        product(1) = -metric(1,2)*(stress(1)-stress(2))-metric(1,3)*(stress(1)-stress(3))
        product(2) = -metric(1,2)*(stress(2)-stress(1))-metric(2,3)*(stress(2)-stress(3))
        product(3) = -metric(1,3)*(stress(3)-stress(1))-metric(2,3)*(stress(3)-stress(2))
        product(4) = metric(4,4)*stress(4)
        product(5) = metric(5,5)*stress(5)
        product(6) = metric(6,6)*stress(6)
    end subroutine apply_hill_metric

    subroutine stress_at_multiplier(t, shear, metric, trial, stress, equivalent, ok)
        !! Solve the Hill plastic corrector for a prescribed multiplier and compute the equivalent stress.
        real(wp), intent(in) :: t
            !! Scalar trial plastic multiplier in the Hill return-map solve.
        real(wp), intent(in) :: shear
            !! Positive shear modulus \(G\) used in the Hill plastic corrector.
        real(wp), intent(in) :: metric(nmandel,nmandel)
            !! Symmetric six-by-six Hill quadratic yield metric.
        real(wp), intent(in) :: trial(nmandel)
            !! Elastic predictor stress before the plastic corrector.
        real(wp), intent(out) :: stress(nmandel)
            !! Corrected Cauchy stress in the model's Mandel convention.
        real(wp), intent(out) :: equivalent
            !! Equivalent stress obtained from the Hill quadratic yield metric.
        logical, intent(out) :: ok
            !! True when the local linear solve and resulting stress state are admissible.
        real(wp) :: a(nmandel,nmandel)
        real(wp) :: rhs(nmandel,1)
        real(wp) :: solved(nmandel,1)
        real(wp) :: pressure
        integer :: i

        a = 2.0_wp*shear*t*metric
        do i = 1, nmandel
            a(i,i) = a(i,i)+1.0_wp
        end do
        rhs(:,1) = trial
        pressure = sum(trial(1:3))/3.0_wp
        rhs(1:3,1) = rhs(1:3,1)-pressure
        call solve_six(a,rhs,solved,ok)
        if (.not. ok) return
        stress = solved(:,1)
        equivalent = sqrt(hill_quadratic(metric,stress))
        stress(1:3) = stress(1:3)+pressure
        ok = ieee_is_finite(equivalent) .and. all(ieee_is_finite(stress))
    end subroutine stress_at_multiplier

    subroutine solve_six(matrix, rhs, solution, ok)
        !! Pivoted dense 6 by 6 solve. All work arrays have fixed stack size.
        real(wp), intent(in) :: matrix(nmandel,nmandel)
            !! Dense coefficient matrix of the local linear system.
        real(wp), intent(in) :: rhs(:,:)
            !! Matrix of right-hand sides for the fixed six-by-six linear solve.
        real(wp), intent(out) :: solution(:,:)
            !! Matrix of solution vectors for each supplied right-hand side.
        logical, intent(out) :: ok
            !! True when the local linear solve and resulting stress state are admissible.
        real(wp) :: a(nmandel,nmandel)
        real(wp) :: work(nmandel,size(rhs,2))
        real(wp) :: swap_row(nmandel)
        real(wp) :: swap_rhs(size(rhs,2))
        real(wp) :: factor
        integer :: k
        integer :: i
        integer :: pivot

        a = matrix
        work = rhs
        ok = .false.
        if (.not. all(ieee_is_finite(a))) return
        do k = 1, nmandel
            pivot = k-1+maxloc(abs(a(k:nmandel,k)),dim=1)
            if (.not. ieee_is_finite(a(pivot,k))) return
            if (abs(a(pivot,k)) <= tiny(1.0_wp)) return
            if (pivot /= k) then
                swap_row = a(k,:)
                a(k,:) = a(pivot,:)
                a(pivot,:) = swap_row
                swap_rhs = work(k,:)
                work(k,:) = work(pivot,:)
                work(pivot,:) = swap_rhs
            end if
            do i = k+1, nmandel
                factor = a(i,k)/a(k,k)
                a(i,k:nmandel) = a(i,k:nmandel)-factor*a(k,k:nmandel)
                work(i,:) = work(i,:)-factor*work(k,:)
            end do
        end do
        solution = work
        do k = nmandel, 1, -1
            do i = k+1, nmandel
                solution(k,:) = solution(k,:)-a(k,i)*solution(i,:)
            end do
            solution(k,:) = solution(k,:)/a(k,k)
        end do
        ok = all(ieee_is_finite(solution))
    end subroutine solve_six

end module constituf_hill_plasticity

! SPDX-License-Identifier: Apache-2.0

module constituf_j2_plane
    !! Two-dimensional small-strain \(J_2\) plasticity under plane constraints.
    !!
    !! This material adapts the *existing* three-dimensional \(J_2\) radial-return
    !! algorithm, preserving its seven internal history components and exact
    !! consistent tangent. The in-plane Mandel input and output have physical
    !! shape \(2\times2\), with the in-plane strain packed as
    !! \(\boldsymbol{\varepsilon}_a=
    !! [\varepsilon_{11},\varepsilon_{22},\sqrt{2}\varepsilon_{12}]^T\)
    !! and Cauchy stress packed in the same order. Strain is dimensionless
    !! and stress shares the moduli's stress unit. Normal tensile components
    !! are positive.
    !!
    !! In plane strain, \(\varepsilon_{33}=\varepsilon_{13}=\varepsilon_{23}=0\)
    !! and the returned tangent is the principal submatrix with Mandel
    !! indices \((1,2,6)\) of the 3D algorithmic tangent.
    !! Out-of-plane stress \(\sigma_{33}\) may be nonzero and is not part
    !! of the 2D response.
    !!
    !! In plane stress, \(\sigma_{33}=\sigma_{13}=\sigma_{23}=0\).
    !! For prescribed in-plane strain \(\boldsymbol{\varepsilon}_a\) and
    !! committed 3D plastic history \(\boldsymbol{s}_n\), the nonlinear
    !! local scalar constraint is
    !! \[R(\varepsilon_{33})=
    !! \sigma_{33}(\boldsymbol{\varepsilon}_a,\varepsilon_{33},
    !! \boldsymbol{s}_n)=0.\]
    !! Newton iterations use the *3D discrete* algorithmic derivative
    !! \(C_{33}=\partial\sigma_{33}/\partial\varepsilon_{33}\), where
    !! \(C_{ij}\) denotes the 3D Mandel tangent matrix and index \(3\)
    !! selects the transverse normal component. Every iteration evaluates
    !! from the same committed plastic state. No trial is committed
    !! prematurely. The accepted in-plane algorithmic tangent is the exact
    !! Schur complement
    !! \[C^{\mathrm{ps}}_{ab}=
    !! C_{ab}-C_{a3}C_{33}^{-1}C_{3b},\qquad
    !! a,b\in\{1,2,6\}.\]
    !! This is not a truncation of the 3D tangent. Out-of-plane strain is
    !! solved consistently with the actual nonlinear return mapping.
    !!
    !! The real state includes the 3D plastic-strain Mandel vector,
    !! accumulated equivalent plastic strain and accepted \(\varepsilon_{33}\).
    !! Successful `material_t%commit` accepts all eight components together.
    !! Revert and failed trials leave committed history unchanged. For a
    !! newly initialized handle the initial state is virgin, with zero
    !! plastic strain and zero out-of-plane strain.
    !!
    !! Inputs: `strain.small` (3 Mandel components, dimensionless).
    !! Response: `stress.cauchy` (3 Mandel components, stress units).
    !! The 3D parameters are bulk modulus \(K>0\), shear modulus \(G>0\),
    !! initial yield stress \(\sigma_y>0\) and hardening modulus \(H\geq0\),
    !! all in the chosen stress unit. Analytical sensitivities with respect
    !! to a selected parameter \(p\) differentiate the underlying 3D return
    !! and, for plane stress, the implicit transverse constraint:
    !! \[\frac{\mathrm{d}\varepsilon_{33}}{\mathrm{d}p}
    !! =-\frac{1}{C_{33}}
    !! \left.\frac{\mathrm{d}\sigma_{33}}{\mathrm{d}p}
    !! \right|_{\varepsilon_{33}}.\]
    !! The stress derivative on the right includes the accepted plastic
    !! history path while holding the current transverse strain fixed.
    !! Prescribed in-plane strains are held fixed throughout.
    !! The transverse correction is also propagated into the complete
    !! plastic state path.
    !! The material is isothermal, rate-independent and associative,
    !! with isotropic hardening and infinitesimal rotations only.
    !! Finite-strain and anisotropic plane-stress laws are not implemented.
    !!
    !! If the local constraint cannot converge within 40 iterations or
    !! the 3D tangent loses its positive scalar pivot \(C_{33}>0\), the point fails
    !! with `cf_evaluation_failed`. No extrapolated or secant tangent is used.
    !!
    !! Underlying 3D return and tangent: J. C. Simo and R. L. Taylor,
    !! "Consistent Tangent Operators for Rate-Independent Elastoplasticity",
    !! Computer Methods in Applied Mechanics and Engineering 48 (1985),
    !! 101--118, DOI: 10.1016/0045-7825(85)90070-2. Sections 3 and 4,
    !! Eqs. (3.8) and (4.12)--(4.13), specialized to zero back stress and
    !! linear isotropic hardening. J. C. Simo and T. J. R. Hughes,
    !! Computational Inelasticity, Springer (1998), DOI: 10.1007/b98904.
    !! The scalar transverse-strain Newton solve and Schur elimination
    !! described above are this adapter's construction around that 3D
    !! discrete update.
    !! Code provenance: independently authored adapter of ConstituF's own
    !! constituf_j2_plasticity. No third-party code was reused.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument, cf_evaluation_failed
    use constituf_j2_plasticity, only: j2_plasticity_t
    implicit none (type, external)
    private

    integer, parameter :: nplane = 3
    integer, parameter :: nmandel = 6
    integer, parameter :: nbase_state = 7
    integer, parameter :: nstate = 8
    integer, parameter :: max_plane_iterations = 40
    integer, parameter :: in_plane(nplane) = [1,2,6]

    public :: j2_plane_t

    type, extends(material_model_t) :: j2_plane_t
        !! Build using `j2_plane_t(bulk_modulus, shear_modulus, yield_stress, hardening_modulus, plane_stress)`.
        !! The final logical defaults to `.false.` (plane strain).
        !! All moduli have stress units and are cloned at `init`.
        real(wp) :: bulk_modulus = 0.0_wp
            !! Positive bulk modulus \(K\).
        real(wp) :: shear_modulus = 0.0_wp
            !! Positive 3D shear modulus \(G\).
        real(wp) :: yield_stress = 0.0_wp
            !! Positive initial von Mises yield stress \(\sigma_y\).
        real(wp) :: hardening_modulus = 0.0_wp
            !! Nonnegative isotropic hardening modulus \(H\).
        logical :: plane_stress = .false.
            !! Enforce zero out-of-plane stress if true.
    contains
        procedure :: describe => describe_j2_plane
            !! Describe the material input, response, state and scalar parameter schema.
        procedure :: eval_batch => evaluate_j2_plane
            !! Evaluate the batched material response and trial history.
        procedure :: initialize_sensitivity => initialize_j2_plane_sensitivity
            !! Initialize analytical derivatives of the reference real history.
        procedure :: eval_sensitivity_batch => sensitivity_j2_plane
            !! Evaluate exact selected-parameter sensitivities through the discrete material update.
    end type j2_plane_t

contains

    subroutine describe_j2_plane(self, spec, status)
        !! Reject invalid physical parameters at model initialization.
        class(j2_plane_t), intent(in) :: self
            !! Model instance containing the fixed constitutive parameters.
        type(material_spec_t), intent(out) :: spec
            !! Material schema receiving the packed field descriptions and supported capabilities.
        integer, intent(out) :: status
            !! Status code; `cf_ok` indicates success.
        type(j2_plasticity_t) :: base
        type(material_spec_t) :: three_d

        base = j2_plasticity_t(self%bulk_modulus,self%shear_modulus,self%yield_stress,self%hardening_modulus)
        call base%describe(three_d,status)
        if (status /= cf_ok) return

        call spec%add_input("strain.small", nplane, "mandel", "small_strain", "spatial", [2,2])
        call spec%add_output("stress.cauchy", nplane, "mandel", "cauchy_stress", "spatial", [2,2])
        call spec%add_state("plastic.strain", nmandel, "mandel", "plastic_strain", "spatial", [3,3])
        call spec%add_state("plastic.accumulated", 1, "scalar", "equivalent_plastic_strain", "material")
        call spec%add_state("plane.out_of_plane_strain", 1, "scalar", "small_strain", "spatial")
        call spec%add_parameter("bulk_modulus")
        call spec%add_parameter("shear_modulus")
        call spec%add_parameter("yield_stress")
        call spec%add_parameter("hardening_modulus")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_j2_plane

    subroutine evaluate_j2_plane(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Allocation-free local Newton solve and consistent Schur tangent.
        !!
        !! All local evaluations start from committed `state_old`. Only
        !! `state_trial` receives the converged result. The 3D \(J_2\) evaluator
        !! is called with one point at a time using fixed-size stack arrays.
        !! Even without a requested 2D tangent, plane stress needs the 3D
        !! algorithmic tangent to converge its out-of-plane constraint.
        class(j2_plane_t), intent(in) :: self
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
        type(j2_plasticity_t) :: base
        real(wp) :: old3(0,0)
        real(wp) :: strain3(nmandel,1)
        real(wp) :: sigma3(nmandel,1)
        real(wp) :: plastic3(nbase_state,1)
        real(wp) :: jac3(nmandel,nmandel,1)
        real(wp) :: e33
        real(wp) :: lambda
        real(wp) :: epsp_trace
        real(wp) :: scale
        real(wp) :: residual
        real(wp) :: pivot
        real(wp) :: condensed(nplane,nplane)
        integer :: local_status(1)
        integer :: ipoint
        integer :: iteration
        integer :: i
        integer :: j
        logical :: converged

        base = j2_plasticity_t(self%bulk_modulus,self%shear_modulus,self%yield_stress,self%hardening_modulus)
        lambda = self%bulk_modulus-2.0_wp*self%shear_modulus/3.0_wp
        point_status = cf_invalid_argument
        point_loop: do ipoint = 1, size(x_new,2)
            if (.not. all(ieee_is_finite(x_new(:,ipoint)))) cycle point_loop
            if (.not. all(ieee_is_finite(state_old(:,ipoint)))) cycle point_loop
            strain3 = 0.0_wp
            strain3(in_plane,1) = x_new(:,ipoint)
            e33 = 0.0_wp
            if (self%plane_stress) then
                epsp_trace = sum(state_old(1:3,ipoint))
                e33 = (2.0_wp*self%shear_modulus*state_old(3,ipoint) - lambda*(x_new(1,ipoint)+x_new(2,ipoint)-epsp_trace))/ (lambda+2.0_wp*self%shear_modulus)
            end if
            if (.not. ieee_is_finite(e33)) cycle point_loop

            converged = .false.
            newton_loop: do iteration = 1, max_plane_iterations
                strain3(3,1) = e33
                local_status = cf_evaluation_failed
                call base%eval_batch(old3,strain3,state_old(1:nbase_state,ipoint:ipoint), plastic3,sigma3,0.0_wp,jac3,local_status)
                if (local_status(1) /= cf_ok) exit newton_loop
                if (self%plane_stress) then
                    residual = sigma3(3,1)
                    ! Normalize by a positive stress scale without an absolute pressure floor.
                    scale = max(self%yield_stress,maxval(abs(sigma3(:,1))))
                    if (abs(residual)/scale <= 128.0_wp*epsilon(1.0_wp)) then
                        converged = .true.
                        exit newton_loop
                    end if
                    pivot = jac3(3,3,1)
                    if (.not. ieee_is_finite(pivot)) exit newton_loop
                    if (pivot <= tiny(1.0_wp)) exit newton_loop
                    e33 = e33-residual/pivot
                    if (.not. ieee_is_finite(e33)) exit newton_loop
                else
                    converged = .true.
                    exit newton_loop
                end if
            end do newton_loop
            if (.not. converged) then
                point_status(ipoint) = cf_evaluation_failed
                cycle point_loop
            end if

            state_trial(1:nbase_state,ipoint) = plastic3(:,1)
            state_trial(nstate,ipoint) = e33
            response(:,ipoint) = sigma3(in_plane,1)
            if (.not. all(ieee_is_finite(state_trial(:,ipoint)))) cycle point_loop
            if (.not. all(ieee_is_finite(response(:,ipoint)))) cycle point_loop
            if (present(tangent)) then
                if (self%plane_stress) then
                    pivot = jac3(3,3,1)
                    if (.not. ieee_is_finite(pivot)) cycle point_loop
                    if (pivot <= tiny(1.0_wp)) cycle point_loop
                end if
                do j = 1, nplane
                    do i = 1, nplane
                        condensed(i,j) = jac3(in_plane(i),in_plane(j),1)
                        if (self%plane_stress) then
                            condensed(i,j) = condensed(i,j) - jac3(in_plane(i),3,1)*jac3(3,in_plane(j),1)/pivot
                        end if
                    end do
                end do
                if (.not. all(ieee_is_finite(condensed))) cycle point_loop
                tangent(:,:,ipoint) = condensed
            end if
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_j2_plane


    subroutine initialize_j2_plane_sensitivity(self,initial_input,initial_state,selected_parameters,sensitivity,status)
        !! Initialize analytical history sensitivities for the virgin material state.
        class(j2_plane_t),intent(in) :: self
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
    end subroutine initialize_j2_plane_sensitivity

    subroutine sensitivity_j2_plane(self,x_old,x_new,state_old,state_trial,selected_parameters, sensitivity_old,sensitivity_trial,response,dsdp,dt,tangent,point_status)
        !! Exact discrete \(J_2\) path derivative plus implicit plane-stress
        !! elimination. All Newton iterations use fixed committed history.
        class(j2_plane_t),intent(in) :: self
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
        type(j2_plasticity_t) :: base
        real(wp) :: old3(0,0)
        real(wp) :: strain3(nmandel,1)
        real(wp) :: sigma3(nmandel,1)
        real(wp) :: plastic3(nbase_state,1)
        real(wp) :: jac3(nmandel,nmandel,1)
        real(wp) :: base_sens(nbase_state,size(selected_parameters),1)
        real(wp) :: base_dsdp(nmandel,size(selected_parameters),1)
        real(wp) :: strain_dev(nmandel)
        real(wp) :: trial_dev(nmandel)
        real(wp) :: normal(nmandel)
        real(wp) :: derivative_trial(nmandel)
        real(wp) :: derivative_normal(nmandel)
        real(wp) :: plastic_input_derivative(nbase_state)
        real(wp) :: identity(nmandel)
        real(wp) :: qtrial
        real(wp) :: yield_function
        real(wp) :: gamma
        real(wp) :: dgamma
        real(wp) :: dq
        real(wp) :: denominator
        real(wp) :: transverse_derivative
        real(wp) :: pivot
        real(wp) :: trace_strain
        integer :: local_status(1)
        integer :: ipoint
        integer :: iparam
        integer :: i
        integer :: j

        call evaluate_j2_plane(self,x_old,x_new,state_old,state_trial,response,dt,tangent,point_status)
        base=j2_plasticity_t(self%bulk_modulus,self%shear_modulus,self%yield_stress,self%hardening_modulus)
        identity=[1.0_wp,1.0_wp,1.0_wp,0.0_wp,0.0_wp,0.0_wp]
        sensitivity_trial=sensitivity_old
        dsdp=0.0_wp
        point_loop: do ipoint=1,size(x_new,2)
            if (point_status(ipoint)/=cf_ok) cycle point_loop
            strain3=0.0_wp
            strain3(in_plane,1)=x_new(:,ipoint)
            strain3(3,1)=state_trial(nstate,ipoint)
            local_status=cf_evaluation_failed
            call base%eval_sensitivity_batch(old3,strain3,state_old(1:nbase_state,ipoint:ipoint), plastic3,selected_parameters,sensitivity_old(1:nbase_state,:,ipoint:ipoint), base_sens,sigma3,base_dsdp,0.0_wp,jac3,local_status)
            if (local_status(1)/=cf_ok) then
                point_status(ipoint)=local_status(1)
                cycle point_loop
            end if

            plastic_input_derivative=0.0_wp
            if (self%plane_stress) then
                trace_strain=sum(strain3(1:3,1))
                strain_dev=strain3(:,1)-trace_strain*identity/3.0_wp
                trial_dev=2.0_wp*self%shear_modulus*(strain_dev-state_old(1:nmandel,ipoint))
                qtrial=sqrt(1.5_wp*dot_product(trial_dev,trial_dev))
                yield_function=qtrial-self%yield_stress- self%hardening_modulus*state_old(nbase_state,ipoint)
                if (yield_function>0.0_wp) then
                    denominator=3.0_wp*self%shear_modulus+self%hardening_modulus
                    gamma=yield_function/denominator
                    normal=1.5_wp*trial_dev/qtrial
                    derivative_trial=-2.0_wp*self%shear_modulus*identity/3.0_wp
                    derivative_trial(3)=derivative_trial(3)+2.0_wp*self%shear_modulus
                    dq=dot_product(normal,derivative_trial)
                    dgamma=dq/denominator
                    derivative_normal=(1.5_wp/qtrial)*(derivative_trial-trial_dev*dq/qtrial)
                    plastic_input_derivative(1:nmandel)=dgamma*normal+gamma*derivative_normal
                    plastic_input_derivative(nbase_state)=dgamma
                end if
                pivot=jac3(3,3,1)
                if (.not. ieee_is_finite(pivot)) then
                    point_status(ipoint)=cf_evaluation_failed
                    cycle point_loop
                end if
                if (pivot<=tiny(1.0_wp)) then
                    point_status(ipoint)=cf_evaluation_failed
                    cycle point_loop
                end if
            end if
            parameter_loop: do iparam=1,size(selected_parameters)
                transverse_derivative=0.0_wp
                if (self%plane_stress) transverse_derivative=-base_dsdp(3,iparam,1)/pivot
                sensitivity_trial(1:nbase_state,iparam,ipoint)=base_sens(:,iparam,1)+ plastic_input_derivative*transverse_derivative
                sensitivity_trial(nstate,iparam,ipoint)=transverse_derivative
                do i=1,nplane
                    j=in_plane(i)
                    dsdp(i,iparam,ipoint)=base_dsdp(j,iparam,1)+ jac3(j,3,1)*transverse_derivative
                end do
            end do parameter_loop
            if (.not. all(ieee_is_finite(dsdp(:,:,ipoint))) .or. .not. all(ieee_is_finite(sensitivity_trial(:,:,ipoint)))) then
                point_status(ipoint)=cf_evaluation_failed
            end if
        end do point_loop
    end subroutine sensitivity_j2_plane

end module constituf_j2_plane

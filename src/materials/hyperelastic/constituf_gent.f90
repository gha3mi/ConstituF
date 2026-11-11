! SPDX-License-Identifier: Apache-2.0

module constituf_gent
    !! Compressible isotropic Gent hyperelasticity in three dimensions.
    !!
    !! For deformation gradient \(\boldsymbol{F}\), \(J=\det\boldsymbol{F}>0\),
    !! \(I_1=\boldsymbol{F}:\boldsymbol{F}\) and
    !! \(\bar I_1=J^{-2/3}I_1\), this implementation uses energy per unit
    !! reference volume
    !! \[\Psi=-\frac{\mu J_m}{2}
    !! \ln\left(1-\frac{\bar I_1-3}{J_m}\right)
    !! +\frac{K}{2}(J-1)^2.\]
    !! Here \(\mu>0\) and \(K>0\) (stress units) are initial shear and bulk moduli.
    !! \(J_m>0\) (dimensionless) is the limiting isochoric invariant excess.
    !! The quadratic compressibility penalty is a ConstituF adaptation,
    !! not part of Gent's original incompressible model.
    !! As \(J_m\to\infty\), the energy tends to
    !! \(\Psi=\mu(\bar I_1-3)/2+K(J-1)^2/2\).
    !!
    !! Let \(s=1-(\bar I_1-3)/J_m>0\),
    !! \(\boldsymbol{A}=\boldsymbol{F}^{-T}\) and
    !! \(\boldsymbol{Q}=\boldsymbol{F}-(I_1/3)\boldsymbol{A}\). Then
    !! \[\boldsymbol{P}=\frac{\mu J^{-2/3}}{s}\boldsymbol{Q}
    !! +K(J-1)J\boldsymbol{A}.\]
    !! \(\boldsymbol{P}=\partial\Psi/\partial\boldsymbol{F}\) is the first
    !! Piola--Kirchhoff stress (stress units), work-conjugate to dimensionless \(\boldsymbol{F}\).
    !! The optional \(9\times9\) material tangent is the exact derivative
    !! \(T_{iJkL}=\partial P_{iJ}/\partial F_{kL}\), at fixed parameters:
    !! \[T_{iJkL}=\frac{\mu J^{-2/3}}{s}
    !! \left[\delta_{ik}\delta_{JL}
    !! -\frac{2}{3}F_{kL}A_{iJ}+\frac{I_1}{3}A_{iL}A_{kJ}
    !! -\frac{2}{3}A_{kL}Q_{iJ}
    !! +\frac{2J^{-2/3}}{J_ms}Q_{iJ}Q_{kL}\right]
    !! +KJ\left[(2J-1)A_{kL}A_{iJ}
    !! -(J-1)A_{iL}A_{kJ}\right].\]
    !! Optional exact *local* sensitivities at fixed \(\boldsymbol{F}\) are
    !! \[\partial_\mu\boldsymbol{P}=J^{-2/3}\boldsymbol{Q}/s,\quad
    !! \partial_{J_m}\boldsymbol{P}=
    !! -\mu J^{-2/3}(\bar I_1-3)\boldsymbol{Q}/(J_m^2s^2),\quad
    !! \partial_K\boldsymbol{P}=(J-1)J\boldsymbol{A}.\]
    !!
    !! Input `deformation.f` and output `stress.piola` are full \(3\times3\)
    !! spatial-reference tensors, flattened in Fortran column-major component-index order
    !! \((11,21,31,12,22,32,13,23,33)\), not Mandel/Voigt.
    !! This stateless, isothermal, rate-independent law has no time-step
    !! dependence and provides no Cauchy stress, energy output or global
    !! equilibrium sensitivities. Inverted deformations, locking limit
    !! \(s\leq0\) and nonrepresentable values fail without accepted
    !! history changes. No general large-strain ellipticity is guaranteed.
    !!
    !! Scientific source (not source code): A. N. Gent, "A New Constitutive
    !! Relation for Rubber", Rubber Chemistry and Technology 69(1) (1996),
    !! 59--61, DOI: 10.5254/1.3538357.
    !! The uncoupled isochoric Gent energy is also documented in FEBio,
    !! "Gent", Description:
    !! https://febiosoftware.github.io/febio-docs/features/features/solid_material_gent/.
    !! This implementation independently combines that isochoric law with
    !! the quadratic volume penalty stated above.
    !! No third-party implementation is copied.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: ndim = 3
    integer, parameter :: nfull = ndim*ndim

    public :: gent_t

    type, extends(material_model_t) :: gent_t
        !! Gent model, constructed as `gent_t(mu,j_m,bulk_modulus)`.
        real(wp) :: shear_modulus = 0.0_wp
            !! Initial shear modulus \(\mu>0\), in the chosen stress unit.
        real(wp) :: limiting_parameter = 0.0_wp
            !! Locking parameter \(J_m>0\), dimensionless.
        real(wp) :: bulk_modulus = 0.0_wp
            !! Initial bulk modulus \(K>0\), in the chosen stress unit.
    contains
        procedure :: describe => describe_gent
            !! Declare input and output fields, parameter names and constitutive capabilities.
        procedure :: eval_batch => evaluate_gent
            !! Evaluate the batched constitutive response and optional consistent tangent.
        procedure :: eval_sensitivity_batch => sensitivity_gent
            !! Evaluate exact scalar-parameter derivatives of the constitutive response.
    end type gent_t

contains

    subroutine describe_gent(self, spec, status)
        !! Describe full \(\boldsymbol{F}\) and \(\boldsymbol{P}\) tensors and exactly three scalar parameters.
        class(gent_t), intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        type(material_spec_t), intent(out) :: spec
            !! Packed field schema, material parameter metadata and declared capabilities.
        integer, intent(out) :: status
            !! Operation status; `cf_ok` indicates success.

        status = cf_invalid_argument
        if (.not. all(ieee_is_finite([self%shear_modulus, self%limiting_parameter,self%bulk_modulus]))) return
        if (self%shear_modulus <= 0.0_wp) return
        if (self%limiting_parameter <= 0.0_wp) return
        if (self%bulk_modulus <= 0.0_wp) return

        call spec%add_input("deformation.f", nfull, "tensor_full", "deformation_gradient", "spatial_reference", [ndim,ndim])
        call spec%add_output("stress.piola", nfull, "tensor_full", "first_piola_kirchhoff", "spatial_reference", [ndim,ndim])
        call spec%add_parameter("shear_modulus")
        call spec%add_parameter("limiting_parameter")
        call spec%add_parameter("bulk_modulus")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_gent

    subroutine evaluate_gent(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Stateless batch evaluation. Caller controls commit and revert.
        class(gent_t), intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed input `(ninput,npoints)` when needed, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `deformation.f` input `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued material history `(nstate,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Trial real-valued material history `(nstate,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `stress.piola` output `(noutput,npoints)`.
        real(wp), intent(in) :: dt
            !! Prescribed time increment; required to be positive only for time-dependent models.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic response Jacobian `(noutput,ninput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Status of each point; a successful update must set its code to `cf_ok`.

        state_trial = state_old
        call evaluate_points(self, x_new, response, point_status, tangent)
    end subroutine evaluate_gent

    subroutine sensitivity_gent(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Exact local derivatives in selected order. Parameter IDs: \(\mu\) is `1`, \(J_m\) is `2` and \(K\) is `3`.
        !! `dsdp` has shape `(9,nselected,npoints)`. No sensitivity state.
        class(gent_t), intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        real(wp), intent(in) :: x_old(:,:)
            !! Previously accepted packed input `(ninput,npoints)` when needed, otherwise an empty array.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `deformation.f` input `(ninput,npoints)`.
        real(wp), intent(in) :: state_old(:,:)
            !! Committed real-valued material history `(nstate,npoints)`.
        integer, intent(in) :: selected_parameters(:)
            !! One-based IDs of scalar parameters in the requested sensitivity-axis order.
        real(wp), intent(in) :: sensitivity_old(:,:,:)
            !! Committed derivatives of real material history `(nstate,nselected,npoints)`.
        real(wp), intent(out) :: state_trial(:,:)
            !! Trial real-valued material history `(nstate,npoints)`.
        real(wp), intent(out) :: sensitivity_trial(:,:,:)
            !! Trial derivatives of real material history `(nstate,nselected,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `stress.piola` output `(noutput,npoints)`.
        real(wp), intent(out) :: dsdp(:,:,:)
            !! Exact discrete response derivatives `(noutput,nselected,npoints)` at fixed input histories.
        real(wp), intent(in) :: dt
            !! Prescribed time increment; required to be positive only for time-dependent models.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic response Jacobian `(noutput,ninput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Status of each point; a successful update must set its code to `cf_ok`.

        state_trial = state_old
        sensitivity_trial = sensitivity_old
        call evaluate_points(self, x_new, response, point_status, tangent, selected_parameters, dsdp)
    end subroutine sensitivity_gent

    subroutine evaluate_points(self, x_new, response, point_status, tangent, selected_parameters, dsdp)
        !! Exact energy gradient/Hessian. No allocation in the point loop.
        class(gent_t), intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        real(wp), intent(in) :: x_new(:,:)
            !! Current prescribed `deformation.f` input `(ninput,npoints)`.
        real(wp), intent(out) :: response(:,:)
            !! Packed `stress.piola` output `(noutput,npoints)`.
        integer, intent(inout) :: point_status(:)
            !! Status of each point; a successful update must set its code to `cf_ok`.
        real(wp), intent(out), optional :: tangent(:,:,:)
            !! Optional exact algorithmic response Jacobian `(noutput,ninput,npoints)`.
        integer, intent(in), optional :: selected_parameters(:)
            !! One-based IDs of scalar parameters in the requested sensitivity-axis order.
        real(wp), intent(out), optional :: dsdp(:,:,:)
            !! Exact discrete response derivatives `(noutput,nselected,npoints)` at fixed input histories.
        real(wp) :: f(ndim,ndim)
        real(wp) :: cof(ndim,ndim)
        real(wp) :: a(ndim,ndim)
        real(wp) :: q(ndim,ndim)
        real(wp) :: p(ndim,ndim)
        real(wp) :: jac
        real(wp) :: jm23
        real(wp) :: i1
        real(wp) :: excess
        real(wp) :: slack
        real(wp) :: factor
        real(wp) :: strengthening
        real(wp) :: delta
        real(wp) :: dq
        real(wp) :: value
        integer :: ipoint
        integer :: i
        integer :: j
        integer :: k
        integer :: l
        integer :: row
        integer :: column
        integer :: iparam

        point_loop: do ipoint = 1, size(x_new,2)
            point_status(ipoint) = cf_invalid_argument
            if (.not. all(ieee_is_finite(x_new(:,ipoint)))) cycle point_loop
            f = reshape(x_new(:,ipoint),[ndim,ndim])
            call cofactor_matrix(f,cof)
            jac = sum(f(1,:)*cof(1,:))
            if (.not. ieee_is_finite(jac)) cycle point_loop
            if (jac <= 0.0_wp) cycle point_loop
            a = cof/jac
            if (.not. all(ieee_is_finite(a))) cycle point_loop
            jm23 = jac**(-2.0_wp/3.0_wp)
            i1 = sum(f*f)
            excess = jm23*i1-3.0_wp
            if (.not. all(ieee_is_finite([jm23,i1,excess]))) cycle point_loop
            slack = 1.0_wp-excess/self%limiting_parameter
            if (.not. ieee_is_finite(slack)) cycle point_loop
            if (slack <= 0.0_wp) cycle point_loop
            factor = jm23/slack
            if (.not. ieee_is_finite(factor)) cycle point_loop
            q = f-(i1/3.0_wp)*a
            p = self%shear_modulus*factor*q + self%bulk_modulus*(jac-1.0_wp)*jac*a
            if (.not. all(ieee_is_finite(p))) cycle point_loop

            if (present(tangent)) then
                strengthening = 2.0_wp*jm23/(self%limiting_parameter*slack)
                if (.not. ieee_is_finite(strengthening)) cycle point_loop
                do j = 1, ndim
                    do i = 1, ndim
                        row = i+ndim*(j-1)
                        do l = 1, ndim
                            do k = 1, ndim
                                column = k+ndim*(l-1)
                                delta = 0.0_wp
                                if (i == k .and. j == l) delta = 1.0_wp
                                dq = delta-(2.0_wp/3.0_wp)*f(k,l)*a(i,j) + (i1/3.0_wp)*a(i,l)*a(k,j)
                                value = self%shear_modulus*factor* (dq-(2.0_wp/3.0_wp)*a(k,l)*q(i,j) + strengthening*q(i,j)*q(k,l)) + self%bulk_modulus*jac* ((2.0_wp*jac-1.0_wp)*a(k,l)*a(i,j) - (jac-1.0_wp)*a(i,l)*a(k,j))
                                tangent(row,column,ipoint) = value
                            end do
                        end do
                    end do
                end do
                if (.not. all(ieee_is_finite(tangent(:,:,ipoint)))) cycle point_loop
            end if

            if (present(selected_parameters)) then
                parameter_loop: do iparam = 1, size(selected_parameters)
                    do j = 1, ndim
                        row = ndim*(j-1)
                        select case (selected_parameters(iparam))
                        case (1)
                            dsdp(row+1:row+ndim,iparam,ipoint) = factor*q(:,j)
                        case (2)
                            dsdp(row+1:row+ndim,iparam,ipoint) = -self%shear_modulus*factor*q(:,j)* (excess/(self%limiting_parameter*slack))/ self%limiting_parameter
                        case (3)
                            dsdp(row+1:row+ndim,iparam,ipoint) = (jac-1.0_wp)*jac*a(:,j)
                        case default
                            exit parameter_loop
                        end select
                    end do
                end do parameter_loop
                if (iparam <= size(selected_parameters)) cycle point_loop
                if (.not. all(ieee_is_finite(dsdp(:,:,ipoint)))) cycle point_loop
            end if
            response(:,ipoint) = reshape(p,[nfull])
            point_status(ipoint) = cf_ok
        end do point_loop
    end subroutine evaluate_points

    pure subroutine cofactor_matrix(f, cof)
        !! Explicit 3D cofactor matrix: \(J\boldsymbol{F}^{-T}\).
        real(wp), intent(in) :: f(ndim,ndim)
            !! Three-dimensional deformation gradient \(\boldsymbol{F}\).
        real(wp), intent(out) :: cof(ndim,ndim)
            !! Cofactor matrix \(\operatorname{cof}\boldsymbol{F}=J\boldsymbol{F}^{-T}\).

        cof(1,1) = f(2,2)*f(3,3)-f(2,3)*f(3,2)
        cof(1,2) = f(2,3)*f(3,1)-f(2,1)*f(3,3)
        cof(1,3) = f(2,1)*f(3,2)-f(2,2)*f(3,1)
        cof(2,1) = f(1,3)*f(3,2)-f(1,2)*f(3,3)
        cof(2,2) = f(1,1)*f(3,3)-f(1,3)*f(3,1)
        cof(2,3) = f(1,2)*f(3,1)-f(1,1)*f(3,2)
        cof(3,1) = f(1,2)*f(2,3)-f(1,3)*f(2,2)
        cof(3,2) = f(1,3)*f(2,1)-f(1,1)*f(2,3)
        cof(3,3) = f(1,1)*f(2,2)-f(1,2)*f(2,1)
    end subroutine cofactor_matrix

end module constituf_gent

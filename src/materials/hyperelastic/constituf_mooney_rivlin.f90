! SPDX-License-Identifier: Apache-2.0

module constituf_mooney_rivlin
    !! Compressible three-dimensional isotropic Mooney--Rivlin hyperelasticity.
    !!
    !! For deformation gradient \(\boldsymbol{F}\), \(J=\det\boldsymbol{F}>0\)
    !! and \(\boldsymbol{C}=\boldsymbol{F}^{T}\boldsymbol{F}\), define
    !! \(I_1=\mathrm{tr}\boldsymbol{C}\) and
    !! \(I_2=[I_1^2-\mathrm{tr}(\boldsymbol{C}^2)]/2\). The adopted energy
    !! density per reference volume is
    !! \[\Psi=c_{10}(J^{-2/3}I_1-3)+c_{01}(J^{-4/3}I_2-3)
    !! +\frac{K}{2}(J-1)^2.\]
    !! This is an explicitly chosen *compressible isochoric/volumetric
    !! extension* of the classical incompressible Mooney--Rivlin theory.
    !! The quadratic volume penalty is not attributed to the original papers.
    !! Parameters \(c_{10}>0\), \(c_{01}\geq0\) and \(K>0\) are in the chosen stress unit.
    !! The infinitesimal shear modulus is \(\mu_0=2(c_{10}+c_{01})\).
    !!
    !! With \(\boldsymbol{A}=\boldsymbol{F}^{-T}\) and
    !! \(\boldsymbol{D}=\boldsymbol{F}\boldsymbol{C}\), set
    !! \[\boldsymbol{Q}_1=\boldsymbol{F}-I_1\boldsymbol{A}/3,\quad
    !! \boldsymbol{Q}_2=I_1\boldsymbol{F}-\boldsymbol{D}
    !! -(2I_2/3)\boldsymbol{A}.\]
    !! The first Piola--Kirchhoff stress, conjugate to \(\boldsymbol{F}\), is
    !! \[\boldsymbol{P}=2c_{10}J^{-2/3}\boldsymbol{Q}_1
    !! +2c_{01}J^{-4/3}\boldsymbol{Q}_2
    !! +K(J-1)J\boldsymbol{A}.\]
    !! The optional exact \(9\times9\) tangent differentiates this *discrete algebraic
    !! response* at fixed parameters: \(T_{iJkL}=\partial P_{iJ}/\partial
    !! F_{kL}\). It is computed analytically, not by finite differences.
    !! The optional local parameter derivatives at fixed deformation are
    !! \[\partial_{c_{10}}\boldsymbol{P}=2J^{-2/3}\boldsymbol{Q}_1,\quad
    !! \partial_{c_{01}}\boldsymbol{P}=2J^{-4/3}\boldsymbol{Q}_2,\quad
    !! \partial_K\boldsymbol{P}=(J-1)J\boldsymbol{A}.\]
    !!
    !! Input `deformation.f` and output `stress.piola` are \(3\times3\) full tensors
    !! flattened in Fortran column-major component-index order
    !! \((11,21,31,12,22,32,13,23,33)\).
    !! Stress is positive in tension, in the chosen stress unit. \(\boldsymbol{F}\) is dimensionless. \(\boldsymbol{P}\) is
    !! spatial-reference. This is objective, stateless, isothermal and
    !! rate-independent. No guaranteed global ellipticity at large strains,
    !! incompressible mixed formulation, energy output, or finite-element
    !! equilibrium sensitivities are provided. Inverted/degenerate or
    !! numerically nonrepresentable trial states are rejected.
    !!
    !! Scientific origins (not code sources):
    !! M. Mooney, "A Theory of Large Elastic Deformation", Journal of
    !! Applied Physics 11 (1940), 582--592, DOI 10.1063/1.1712836.
    !! R. S. Rivlin, "Large Elastic Deformations of Isotropic Materials IV",
    !! Philosophical Transactions of the Royal Society A 241 (1948),
    !! 379--397, DOI 10.1098/rsta.1948.0024.
    !! The adopted compressible energy is also given in Dassault Systèmes,
    !! *Abaqus 2026 Material Reference*, "Hyperelastic Behavior of Rubberlike
    !! Materials", "Mooney-Rivlin Form":
    !! https://docs.software.vt.edu/abaqusv2026/English/SIMACAEMATRefMap/simamat-c-hyperelastic.htm.
    !! Its parameters correspond to \(K=2/D_1\) and \(J^{\mathrm{el}}=J\)
    !! here because there is no thermal deformation.
    !! The implementation is independently derived from the stated energy.
    !! No published implementation code is imported.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: ndim = 3
    integer, parameter :: nfull = ndim*ndim

    public :: mooney_rivlin_t

    type, extends(material_model_t) :: mooney_rivlin_t
        !! Stateless finite-strain material with a quadratic volumetric penalty.
        real(wp) :: c10 = 0.0_wp
            !! First modified-invariant coefficient \(c_{10}>0\), in the chosen stress unit.
        real(wp) :: c01 = 0.0_wp
            !! Second modified-invariant coefficient \(c_{01}\geq0\), in the chosen stress unit.
        real(wp) :: bulk_modulus = 0.0_wp
            !! Initial volumetric bulk modulus \(K>0\), in the chosen stress unit.
    contains
        procedure :: describe => describe_mooney_rivlin
            !! Declare input and output fields, parameter names and constitutive capabilities.
        procedure :: eval_batch => evaluate_mooney_rivlin
            !! Evaluate the batched constitutive response and optional consistent tangent.
        procedure :: eval_sensitivity_batch => sensitivity_mooney_rivlin
            !! Evaluate exact scalar-parameter derivatives of the constitutive response.
    end type mooney_rivlin_t

contains

    subroutine describe_mooney_rivlin(self, spec, status)
        !! Describe the full deformation/first-Piola pair and scalar parameters.
        class(mooney_rivlin_t), intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        type(material_spec_t), intent(out) :: spec
            !! Packed field schema, material parameter metadata and declared capabilities.
        integer, intent(out) :: status
            !! Operation status; `cf_ok` indicates success.
        real(wp) :: shear_sum

        status = cf_invalid_argument
        if (.not. all(ieee_is_finite([self%c10,self%c01,self%bulk_modulus]))) return
        if (self%c10 <= 0.0_wp .or. self%c01 < 0.0_wp) return
        if (self%bulk_modulus <= 0.0_wp) return
        shear_sum = self%c10 + self%c01
        if (.not. ieee_is_finite(shear_sum)) return
        if (.not. ieee_is_finite(4.0_wp*shear_sum)) return

        call spec%add_input("deformation.f", nfull, "tensor_full", "deformation_gradient", "spatial_reference", [ndim,ndim])
        call spec%add_output("stress.piola", nfull, "tensor_full", "first_piola_kirchhoff", "spatial_reference", [ndim,ndim])
        call spec%add_parameter("c10")
        call spec%add_parameter("c01")
        call spec%add_parameter("bulk_modulus")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status = cf_ok
    end subroutine describe_mooney_rivlin

    subroutine evaluate_mooney_rivlin(self, x_old, x_new, state_old, state_trial, response, dt, tangent, point_status)
        !! Evaluate a stateless 3D deformation batch. Commit controls the handle.
        class(mooney_rivlin_t), intent(in) :: self
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
    end subroutine evaluate_mooney_rivlin

    subroutine sensitivity_mooney_rivlin(self, x_old, x_new, state_old, state_trial, selected_parameters, sensitivity_old, sensitivity_trial, response, dsdp, dt, tangent, point_status)
        !! Exact local derivatives in the order selected during material init.
        !!
        !! Parameter indices are `c10=1`, `c01=2` and `bulk_modulus=3`. A zero-sized
        !! sensitivity history is preserved because the law is stateless.
        class(mooney_rivlin_t), intent(in) :: self
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
    end subroutine sensitivity_mooney_rivlin

    subroutine evaluate_points(self, x_new, response, point_status, tangent, selected_parameters, dsdp)
        !! Analytic response and optional Hessian without allocations or I/O.
        class(mooney_rivlin_t), intent(in) :: self
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
        real(wp) :: c(ndim,ndim)
        real(wp) :: b(ndim,ndim)
        real(wp) :: d(ndim,ndim)
        real(wp) :: q1(ndim,ndim)
        real(wp) :: q2(ndim,ndim)
        real(wp) :: p(ndim,ndim)
        real(wp) :: jac
        real(wp) :: jm23
        real(wp) :: jm43
        real(wp) :: i1
        real(wp) :: i2
        real(wp) :: dq1
        real(wp) :: dq2
        real(wp) :: dd
        real(wp) :: delta
        real(wp) :: di2
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
            call cofactor_matrix(f, cof)
            jac = sum(f(1,:)*cof(1,:))
            if (.not. ieee_is_finite(jac)) cycle point_loop
            if (jac <= 0.0_wp) cycle point_loop
            a = cof/jac
            if (.not. all(ieee_is_finite(a))) cycle point_loop
            jm23 = jac**(-2.0_wp/3.0_wp)
            jm43 = jac**(-4.0_wp/3.0_wp)
            if (.not. all(ieee_is_finite([jm23,jm43]))) cycle point_loop

            c = matmul(transpose(f),f)
            b = matmul(f,transpose(f))
            d = matmul(f,c)
            i1 = sum(f*f)
            i2 = 0.5_wp*(i1*i1-sum(c*c))
            if (.not. all(ieee_is_finite([i1,i2]))) cycle point_loop
            q1 = f - (i1/3.0_wp)*a
            q2 = i1*f - d - (2.0_wp*i2/3.0_wp)*a
            p = 2.0_wp*self%c10*jm23*q1 + 2.0_wp*self%c01*jm43*q2 + self%bulk_modulus*(jac-1.0_wp)*jac*a
            if (.not. all(ieee_is_finite(p))) cycle point_loop

            if (present(tangent)) then
                do j = 1, ndim
                    do i = 1, ndim
                        row = i + ndim*(j-1)
                        do l = 1, ndim
                            do k = 1, ndim
                                column = k + ndim*(l-1)
                                delta = 0.0_wp
                                if (i == k .and. j == l) delta = 1.0_wp
                                dd = 0.0_wp
                                if (i == k) dd = c(l,j)
                                dd = dd + f(i,l)*f(k,j)
                                if (j == l) dd = dd + b(i,k)
                                di2 = i1*f(k,l)-d(k,l)
                                dq1 = delta - (2.0_wp/3.0_wp)*f(k,l)*a(i,j) + (i1/3.0_wp)*a(i,l)*a(k,j)
                                dq2 = 2.0_wp*f(k,l)*f(i,j) + i1*delta - dd - (4.0_wp/3.0_wp)*di2*a(i,j) + (2.0_wp/3.0_wp)*i2*a(i,l)*a(k,j)
                                value = 2.0_wp*self%c10*jm23* (dq1-(2.0_wp/3.0_wp)*a(k,l)*q1(i,j)) + 2.0_wp*self%c01*jm43* (dq2-(4.0_wp/3.0_wp)*a(k,l)*q2(i,j)) + self%bulk_modulus*jac* ((2.0_wp*jac-1.0_wp)*a(k,l)*a(i,j) - (jac-1.0_wp)*a(i,l)*a(k,j))
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
                            dsdp(row+1:row+ndim,iparam,ipoint) = 2.0_wp*jm23*q1(:,j)
                        case (2)
                            dsdp(row+1:row+ndim,iparam,ipoint) = 2.0_wp*jm43*q2(:,j)
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

end module constituf_mooney_rivlin

! SPDX-License-Identifier: Apache-2.0

module constituf_ogden
    !! Compressible, isotropic, single-term Ogden hyperelasticity in 3D.
    !!
    !! Let \(\boldsymbol{F}\) denote the deformation gradient with
    !! \(J=\det\boldsymbol{F}>0\). Its principal stretches are
    !! \(\lambda_i>0\), with isochoric stretches
    !! \(\bar{\lambda}_i=J^{-1/3}\lambda_i\).
    !! The adopted energy per unit reference volume is
    !! \[\Psi(\boldsymbol{F})
    !! =\frac{2\mu}{a^2}
    !! \left(\sum_{i=1}^3\bar{\lambda}_i^a-3\right)
    !! +\frac{K}{2}(J-1)^2.\]
    !! Here \(\mu>0\) and \(K>0\) have units of stress, and the
    !! dimensionless exponent \(a\) satisfies \(|a|\ge0.5\).
    !! Negative exponents are supported. The \(a\to0\) limit and
    !! multiple Ogden spectral terms are not implemented.
    !!
    !! Define \(\boldsymbol{C}=\boldsymbol{F}^{T}\boldsymbol{F}\),
    !! \(p=a/2\), \(\boldsymbol{M}=\boldsymbol{C}^{p-1}\),
    !! \(s=\mathrm{tr}(\boldsymbol{C}^{p})\),
    !! \(\boldsymbol{A}=\boldsymbol{F}^{-T}\) and
    !! \(\boldsymbol{H}=\boldsymbol{F}\boldsymbol{M}
    !! -s\boldsymbol{A}/3\). The first Piola--Kirchhoff stress is
    !! \[\boldsymbol{P}
    !! =\frac{2\mu}{a}J^{-a/3}\boldsymbol{H}
    !! +K(J-1)J\boldsymbol{A}.\]
    !! Its exact material Jacobian
    !! \(\mathbb{A}_{\mathrm{alg}}=
    !! \partial\boldsymbol{P}/\partial\boldsymbol{F}\)
    !! uses the spectral Frechet derivative of
    !! \(\boldsymbol{C}^{p-1}\). In the eigenbasis of
    !! \(\boldsymbol{C}\), divided differences of
    !! \(f(x)=x^{p-1}\) multiply the components of
    !! \(\mathrm{d}\boldsymbol{C}\). Repeated eigenvalues use
    !! the analytic limit \(f'(x)\). Nearby distinct eigenvalues use
    !! cancellation-free divided differences rather than a repeated-root
    !! approximation.
    !! A local \(3\times3\) Jacobi eigensolver performs the spectral
    !! decomposition without a LAPACK dependency or finite differences.
    !! Ill-conditioned deformations may limit eigenvalue accuracy.
    !!
    !! Partial analytical parameter derivatives at fixed
    !! \(\boldsymbol{F}\) are available for shear modulus \(\mu\),
    !! exponent \(a\) and bulk modulus \(K\). The input field
    !! `deformation.f` (dimensionless) and output field
    !! `stress.piola` (stress units) are full \(3\times3\) tensors,
    !! packed in Fortran column-major component-index order
    !! \((11,21,31,12,22,32,13,23,33)\) in the
    !! spatial-reference frame. Nonfinite or inverted deformations and
    !! non-positive-definite trial \(\boldsymbol{C}\) fail at the
    !! affected point. The material is stateless, isothermal and
    !! rate-independent.
    !!
    !! Reference: R. W. Ogden, *Large Deformation Isotropic Elasticity:
    !! On the Correlation of Theory and Experiment for Incompressible
    !! Rubberlike Solids*, Proceedings of the Royal Society A **326**
    !! (1972), 565--584, DOI: 10.1098/rspa.1972.0026. The volumetric
    !! penalty and spectral implementation are independent adaptations.
    !! The normalization matches Dassault Systèmes, *Abaqus 2026 Material
    !! Reference*, "Hyperelastic Behavior of Rubberlike Materials", "Ogden Form",
    !! with one term, initial shear modulus \(\mu\) and \(K=2/D_1\):
    !! https://docs.software.vt.edu/abaqusv2026/English/SIMACAEMATRefMap/simamat-c-hyperelastic.htm.
    !! For the distinct/repeated eigenvalue derivative rule see N. J. Higham,
    !! "What Is a Fréchet Derivative?" (2020):
    !! https://nhigham.com/2020/06/23/what-is-a-frechet-derivative/.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer,parameter :: ndim=3
    integer,parameter :: nfull=9
    public :: ogden_t

    type, extends(material_model_t) :: ogden_t
        !! One Ogden term. Construct as `ogden_t(mu,exponent,bulk_modulus)`.
        real(wp) :: shear_modulus=0.0_wp
            !! \(\mu>0\), infinitesimal shear modulus (stress units).
        real(wp) :: exponent=0.0_wp
            !! General spectral exponent \(a\), \(|a|\geq0.5\) (dimensionless).
        real(wp) :: bulk_modulus=0.0_wp
            !! \(K>0\), quadratic volumetric penalty (stress units).
    contains
        procedure :: describe => describe_ogden
            !! Declare input and output fields, parameter names and constitutive capabilities.
        procedure :: eval_batch => evaluate_ogden
            !! Evaluate the batched constitutive response and optional consistent tangent.
        procedure :: eval_sensitivity_batch => sensitivity_ogden
            !! Evaluate exact scalar-parameter derivatives of the constitutive response.
    end type ogden_t

contains

    subroutine describe_ogden(self,spec,status)
        !! Validate constitutive parameters and register `deformation.f` input and `stress.piola` response fields.
        class(ogden_t),intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        type(material_spec_t),intent(out) :: spec
            !! Packed field schema, material parameter metadata and declared capabilities.
        integer,intent(out) :: status
            !! Operation status; `cf_ok` indicates success.

        status=cf_invalid_argument
        if (.not. all(ieee_is_finite([self%shear_modulus,self%exponent,self%bulk_modulus]))) return
        if (self%shear_modulus<=0.0_wp.or.self%bulk_modulus<=0.0_wp) return
        if (abs(self%exponent)<0.5_wp) return
        call spec%add_input("deformation.f", nfull, "tensor_full", "deformation_gradient", "spatial_reference", [ndim,ndim])
        call spec%add_output("stress.piola", nfull, "tensor_full", "first_piola_kirchhoff", "spatial_reference", [ndim,ndim])
        call spec%add_parameter("shear_modulus")
        call spec%add_parameter("exponent")
        call spec%add_parameter("bulk_modulus")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status=cf_ok
    end subroutine describe_ogden

    subroutine evaluate_ogden(self,x_old,x_new,state_old,state_trial,response,dt,tangent,point_status)
        !! Evaluate `stress.piola` from `deformation.f` over the current point batch.
        class(ogden_t),intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        real(wp),intent(in) :: x_old(:,:)
            !! Previously accepted packed input `(ninput,npoints)` when needed, otherwise an empty array.
        real(wp),intent(in) :: x_new(:,:)
            !! Current prescribed `deformation.f` input `(ninput,npoints)`.
        real(wp),intent(in) :: state_old(:,:)
            !! Committed real-valued material history `(nstate,npoints)`.
        real(wp),intent(out) :: state_trial(:,:)
            !! Trial real-valued material history `(nstate,npoints)`.
        real(wp),intent(out) :: response(:,:)
            !! Packed `stress.piola` output `(noutput,npoints)`.
        real(wp),intent(in) :: dt
            !! Prescribed time increment; required to be positive only for time-dependent models.
        real(wp),intent(out),optional :: tangent(:,:,:)
            !! Optional exact algorithmic response Jacobian `(noutput,ninput,npoints)`.
        integer,intent(inout) :: point_status(:)
            !! Status of each point; a successful update must set its code to `cf_ok`.

        state_trial=state_old
        call evaluate_points(self,x_new,response,point_status,tangent)
    end subroutine evaluate_ogden

    subroutine sensitivity_ogden(self,x_old,x_new,state_old,state_trial,selected_parameters,sensitivity_old,sensitivity_trial,response,dsdp, dt,tangent,point_status)
        !! Evaluate exact derivatives of `stress.piola` with respect to the selected scalar parameters.
        class(ogden_t),intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        real(wp),intent(in) :: x_old(:,:)
            !! Previously accepted packed input `(ninput,npoints)` when needed, otherwise an empty array.
        real(wp),intent(in) :: x_new(:,:)
            !! Current prescribed `deformation.f` input `(ninput,npoints)`.
        real(wp),intent(in) :: state_old(:,:)
            !! Committed real-valued material history `(nstate,npoints)`.
        integer,intent(in) :: selected_parameters(:)
            !! One-based IDs of scalar parameters in the requested sensitivity-axis order.
        real(wp),intent(in) :: sensitivity_old(:,:,:)
            !! Committed derivatives of real material history `(nstate,nselected,npoints)`.
        real(wp),intent(out) :: state_trial(:,:)
            !! Trial real-valued material history `(nstate,npoints)`.
        real(wp),intent(out) :: sensitivity_trial(:,:,:)
            !! Trial derivatives of real material history `(nstate,nselected,npoints)`.
        real(wp),intent(out) :: response(:,:)
            !! Packed `stress.piola` output `(noutput,npoints)`.
        real(wp),intent(out) :: dsdp(:,:,:)
            !! Exact discrete response derivatives `(noutput,nselected,npoints)` at fixed input histories.
        real(wp),intent(in) :: dt
            !! Prescribed time increment; required to be positive only for time-dependent models.
        real(wp),intent(out),optional :: tangent(:,:,:)
            !! Optional exact algorithmic response Jacobian `(noutput,ninput,npoints)`.
        integer,intent(inout) :: point_status(:)
            !! Status of each point; a successful update must set its code to `cf_ok`.

        state_trial=state_old
        sensitivity_trial=sensitivity_old
        call evaluate_points(self,x_new,response,point_status,tangent,selected_parameters,dsdp)
    end subroutine sensitivity_ogden

    subroutine evaluate_points(self,x,response,point_status,tangent,selected_parameters,dsdp)
        !! Compute `stress.piola` and optional exact input and parameter derivatives for each point.
        class(ogden_t),intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        real(wp),intent(in) :: x(:,:)
            !! Current prescribed `deformation.f` inputs for the point batch.
        real(wp),intent(out) :: response(:,:)
            !! Packed `stress.piola` output `(noutput,npoints)`.
        integer,intent(inout) :: point_status(:)
            !! Status of each point; a successful update must set its code to `cf_ok`.
        real(wp),intent(out),optional :: tangent(:,:,:)
            !! Optional exact algorithmic response Jacobian `(noutput,ninput,npoints)`.
        integer,intent(in),optional :: selected_parameters(:)
            !! One-based IDs of scalar parameters in the requested sensitivity-axis order.
        real(wp),intent(out),optional :: dsdp(:,:,:)
            !! Exact discrete response derivatives `(noutput,nselected,npoints)` at fixed input histories.
        real(wp) :: f(ndim,ndim)
        real(wp) :: a(ndim,ndim)
        real(wp) :: cof(ndim,ndim)
        real(wp) :: c(ndim,ndim)
        real(wp) :: v(ndim,ndim)
        real(wp) :: lam(ndim)
        real(wp) :: power(ndim)
        real(wp) :: pmat(ndim,ndim)
        real(wp) :: dpower(ndim,ndim)
        real(wp) :: fm(ndim,ndim)
        real(wp) :: h(ndim,ndim)
        real(wp) :: stress(ndim,ndim)
        real(wp) :: e(ndim,ndim)
        real(wp) :: dd(ndim,ndim)
        real(wp) :: diff(ndim,ndim)
        real(wp) :: fv(ndim,ndim)
        real(wp) :: jac
        real(wp) :: p
        real(wp) :: alpha
        real(wp) :: w
        real(wp) :: s
        real(wp) :: sp
        real(wp) :: val
        real(wp) :: volume
        real(wp) :: der(ndim,ndim)
        real(wp) :: dtemp(ndim,ndim)
        real(wp) :: power_exponent
        real(wp) :: eigen_gap
        real(wp) :: eigen_delta
        integer :: ipoint
        integer :: i
        integer :: j
        integer :: k
        integer :: l
        integer :: m
        integer :: n
        integer :: row
        integer :: col
        integer :: iparam
        logical :: converged

        alpha=self%exponent
        p=alpha/2.0_wp
        power_exponent=p-1.0_wp
        point_loop: do ipoint=1,size(x,2)
            point_status(ipoint)=cf_invalid_argument
            if (.not. all(ieee_is_finite(x(:,ipoint)))) cycle point_loop
            f=reshape(x(:,ipoint),[ndim,ndim])
            call cofactor_matrix(f,cof)
            jac=sum(f(1,:)*cof(1,:))
            if (.not. ieee_is_finite(jac)) cycle point_loop
            if (jac<=0.0_wp) cycle point_loop
            a=cof/jac
            if (.not. all(ieee_is_finite(a))) cycle point_loop
            c=matmul(transpose(f),f)
            call symmetric_eigen(c,lam,v,converged)
            if (.not. converged) cycle point_loop
            if (any(lam<=0.0_wp).or..not.all(ieee_is_finite(lam))) cycle point_loop

            power=lam**(p-1.0_wp)
            s=sum(lam**p)
            sp=0.5_wp*sum(lam**p*log(lam))
            w=2.0_wp*self%shear_modulus/alpha*jac**(-alpha/3.0_wp)
            if (.not. all(ieee_is_finite([s,sp,w]))) cycle point_loop
            pmat=0.0_wp
            dpower=0.0_wp
            do k=1,ndim
                do j=1,ndim
                    do i=1,ndim
                        pmat(i,j)=pmat(i,j)+power(k)*v(i,k)*v(j,k)
                        dpower(i,j)=dpower(i,j)+0.5_wp*power(k)*log(lam(k))*v(i,k)*v(j,k)
                    end do
                end do
            end do
            fm=matmul(f,pmat)
            h=fm-s*a/3.0_wp
            stress=w*h+self%bulk_modulus*(jac-1.0_wp)*jac*a
            if (.not. all(ieee_is_finite(stress))) cycle point_loop

            if (present(tangent)) then
                diff=0.0_wp
                do j=1,ndim
                    do i=1,j
                        eigen_gap=lam(i)-lam(j)
                        if (abs(power_exponent)<=0.0_wp) then
                            diff(i,j)=0.0_wp
                        else if (abs(eigen_gap)<=0.0_wp) then
                            diff(i,j)=power_exponent*lam(i)**(power_exponent-1.0_wp)
                        else if (abs(eigen_gap)<=0.5_wp*max(lam(i),lam(j))) then
                            eigen_delta=atanh((0.5_wp*eigen_gap)/(0.5_wp*lam(i)+0.5_wp*lam(j)))
                            if (abs(power_exponent*eigen_delta)<=0.5_wp) then
                                ! Exact divided difference without subtracting nearby powers.
                                diff(i,j)=lam(j)**(power_exponent-1.0_wp)*exp((power_exponent-1.0_wp)*eigen_delta)*(sinh(power_exponent*eigen_delta)/sinh(eigen_delta))
                            else
                                diff(i,j)=(power(i)-power(j))/eigen_gap
                            end if
                        else
                            diff(i,j)=(power(i)-power(j))/eigen_gap
                        end if
                        diff(j,i)=diff(i,j)
                    end do
                end do
                fv=matmul(f,v)
                do l=1,ndim
                    do k=1,ndim
                        col=k+ndim*(l-1)
                        do n=1,ndim
                            do m=1,ndim
                                e(m,n)=v(l,m)*fv(k,n)+fv(k,m)*v(l,n)
                            end do
                        end do
                        dd=matmul(v,matmul(diff*e,transpose(v)))
                        der=matmul(f,dd)
                        do j=1,ndim
                            do i=1,ndim
                                row=i+ndim*(j-1)
                                val=der(i,j)
                                if (i==k) val=val+pmat(l,j)
                                val=val-(alpha/3.0_wp)*fm(k,l)*a(i,j)+ (s/3.0_wp)*a(i,l)*a(k,j)- (alpha/3.0_wp)*a(k,l)*h(i,j)
                                volume=self%bulk_modulus*jac* ((2.0_wp*jac-1.0_wp)*a(k,l)*a(i,j)- (jac-1.0_wp)*a(i,l)*a(k,j))
                                tangent(row,col,ipoint)=w*val+volume
                            end do
                        end do
                    end do
                end do
                if (.not. all(ieee_is_finite(tangent(:,:,ipoint)))) cycle point_loop
            end if

            if (present(selected_parameters)) then
                dtemp=matmul(f,dpower)-sp*a/3.0_wp
                do iparam=1,size(selected_parameters)
                    select case(selected_parameters(iparam))
                    case(1)
                        dsdp(:,iparam,ipoint)=reshape((w/self%shear_modulus)*h,[nfull])
                    case(2)
                        dsdp(:,iparam,ipoint)=reshape(w*(dtemp- (1.0_wp/alpha+log(jac)/3.0_wp)*h),[nfull])
                    case(3)
                        dsdp(:,iparam,ipoint)=reshape((jac-1.0_wp)*jac*a,[nfull])
                    case default
                        cycle point_loop
                    end select
                end do
                if (.not. all(ieee_is_finite(dsdp(:,:,ipoint)))) cycle point_loop
            end if
            response(:,ipoint)=reshape(stress,[nfull])
            point_status(ipoint)=cf_ok
        end do point_loop
    end subroutine evaluate_points

    pure subroutine symmetric_eigen(c,lam,v,converged)
        !! Cyclic Jacobi diagonalization with eigenvectors in columns of `v`.
        real(wp),intent(in) :: c(ndim,ndim)
            !! Symmetric right Cauchy--Green deformation tensor \(\boldsymbol{C}=\boldsymbol{F}^{T}\boldsymbol{F}\).
        real(wp),intent(out) :: lam(ndim)
            !! Eigenvalues of the input symmetric tensor.
        real(wp),intent(out) :: v(ndim,ndim)
            !! Orthogonal eigenvectors of the input symmetric tensor.
        logical,intent(out) :: converged
            !! True when the symmetric eigendecomposition converged.
        real(wp) :: mat(ndim,ndim)
        real(wp) :: theta
        real(wp) :: t
        real(wp) :: cs
        real(wp) :: sn
        real(wp) :: app
        real(wp) :: aqq
        real(wp) :: apq
        real(wp) :: akp
        real(wp) :: akq
        real(wp) :: vkp
        real(wp) :: vkq
        real(wp) :: off
        real(wp) :: tol
        integer :: sweep
        integer :: i
        integer :: j
        integer :: k
        mat=c
        v=0.0_wp
        do i=1,ndim
            v(i,i)=1.0_wp
        end do
        converged=.false.
        sweep_loop: do sweep=1,64
            off=max(abs(mat(1,2)),abs(mat(1,3)),abs(mat(2,3)))
            tol=epsilon(1.0_wp)*max(maxval(abs(mat)),tiny(1.0_wp))
            if (off<=tol) then
                converged=.true.
                exit sweep_loop
            end if
            do i=1,ndim-1
                plane_loop: do j=i+1,ndim
                    apq=mat(i,j)
                    if (abs(apq)<=tol) cycle plane_loop
                    theta=(mat(j,j)-mat(i,i))/(2.0_wp*apq)
                    t=sign(1.0_wp,theta)/(abs(theta)+sqrt(1.0_wp+theta**2))
                    cs=1.0_wp/sqrt(1.0_wp+t**2)
                    sn=t*cs
                    app=mat(i,i)
                    aqq=mat(j,j)
                    mat(i,i)=app-t*apq
                    mat(j,j)=aqq+t*apq
                    mat(i,j)=0.0_wp
                    mat(j,i)=0.0_wp
                    do k=1,ndim
                        if (k/=i.and.k/=j) then
                            akp=mat(k,i)
                            akq=mat(k,j)
                            mat(k,i)=cs*akp-sn*akq
                            mat(i,k)=mat(k,i)
                            mat(k,j)=sn*akp+cs*akq
                            mat(j,k)=mat(k,j)
                        end if
                        vkp=v(k,i)
                        vkq=v(k,j)
                        v(k,i)=cs*vkp-sn*vkq
                        v(k,j)=sn*vkp+cs*vkq
                    end do
                end do plane_loop
            end do
        end do sweep_loop
        do i=1,ndim
            lam(i)=mat(i,i)
        end do
    end subroutine symmetric_eigen

    pure subroutine cofactor_matrix(f,cof)
        !! Evaluate the cofactor matrix of a three-dimensional deformation gradient.
        real(wp),intent(in) :: f(ndim,ndim)
            !! Three-dimensional deformation gradient \(\boldsymbol{F}\).
        real(wp),intent(out) :: cof(ndim,ndim)
            !! Cofactor matrix \(\operatorname{cof}\boldsymbol{F}=J\boldsymbol{F}^{-T}\).
        cof(1,1)=f(2,2)*f(3,3)-f(2,3)*f(3,2)
        cof(1,2)=f(2,3)*f(3,1)-f(2,1)*f(3,3)
        cof(1,3)=f(2,1)*f(3,2)-f(2,2)*f(3,1)
        cof(2,1)=f(1,3)*f(3,2)-f(1,2)*f(3,3)
        cof(2,2)=f(1,1)*f(3,3)-f(1,3)*f(3,1)
        cof(2,3)=f(1,2)*f(3,1)-f(1,1)*f(3,2)
        cof(3,1)=f(1,2)*f(2,3)-f(1,3)*f(2,2)
        cof(3,2)=f(1,3)*f(2,1)-f(1,1)*f(2,3)
        cof(3,3)=f(1,1)*f(2,2)-f(1,2)*f(2,1)
    end subroutine cofactor_matrix

end module constituf_ogden

! SPDX-License-Identifier: Apache-2.0

module constituf_yeoh
    !! Compressible, isotropic, three-term Yeoh hyperelasticity in 3D.
    !!
    !! For the deformation gradient \(\boldsymbol{F}\) and
    !! \(J=\det\boldsymbol{F}>0\), define
    !! \(I_1=\boldsymbol{F}:\boldsymbol{F}\),
    !! \(\bar I_1=J^{-2/3}I_1\) and \(s=\bar I_1-3\ge0\).
    !! The adopted strain energy per reference volume is
    !! \[\Psi=c_{10}s+c_{20}s^2+c_{30}s^3
    !! +\frac{K}{2}(J-1)^2.\]
    !! Parameters \(c_{10}>0\), \(c_{20}\ge0\), \(c_{30}\ge0\) and
    !! \(K>0\) all have units of stress. The nonnegative higher-order coefficients
    !! restrict calibration to the hardening branch of this model.
    !! The infinitesimal shear modulus is \(\mu_0=2c_{10}\).
    !! Setting \(c_{20}=c_{30}=0\) gives an isochoric neo-Hookean energy
    !! with the same quadratic volume penalty.
    !! With \(w=c_{10}+2c_{20}s+3c_{30}s^2\) and
    !! \(\boldsymbol{Q}=\boldsymbol{F}-
    !! (I_1/3)\boldsymbol{F}^{-T}\), the first
    !! Piola--Kirchhoff stress is
    !! \[\boldsymbol{P}=2wJ^{-2/3}\boldsymbol{Q}
    !! +K(J-1)J\boldsymbol{F}^{-T}.\]
    !! The exact material tangent is
    !! \(\mathbb{A}_{\mathrm{alg}}=
    !! \partial\boldsymbol{P}/\partial\boldsymbol{F}\)
    !! at fixed material parameters. No numerical-difference
    !! approximation is used. Analytical partial derivatives with
    !! respect to \(c_{10}\), \(c_{20}\), \(c_{30}\) and \(K\)
    !! are offered when selected during initialization. The output
    !! sensitivity axis follows the *requested* parameter order.
    !!
    !! Input `deformation.f` (dimensionless) and output `stress.piola`
    !! (stress units) are full \(3\times3\) tensors packed in Fortran column-major
    !! component-index order \((11,21,31,12,22,32,13,23,33)\), with the
    !! declared spatial-reference frame. The model is stateless,
    !! isothermal and rate-independent. Nonfinite trial fields or
    !! \(J\le0\) fail at the affected point.
    !!
    !! Reference: O. H. Yeoh, *Some Forms of the Strain Energy Function
    !! for Rubber*, Rubber Chemistry and Technology **66** (1993),
    !! 754--771, DOI: 10.5254/1.3538343. The isochoric/volumetric
    !! split and restricted coefficient domain are ConstituF adaptations.
    !! The compressible formulation is documented in Dassault Systèmes,
    !! *Abaqus 2026 Material Reference*, "Hyperelastic Behavior of Rubberlike
    !! Materials", "Yeoh Form":
    !! https://docs.software.vt.edu/abaqusv2026/English/SIMACAEMATRefMap/simamat-c-hyperelastic.htm.
    !! Here \(J^{\mathrm{el}}=J\), \(1/D_1=K/2\) and
    !! \(1/D_2=1/D_3=0\), retaining only the quadratic volume term.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constituf_author, only: wp, material_model_t, material_spec_t, cf_ok, cf_invalid_argument
    implicit none (type, external)
    private

    integer, parameter :: ndim=3
    integer, parameter :: nfull=9
    public :: yeoh_t

    type, extends(material_model_t) :: yeoh_t
        !! Three-term Yeoh parameters in the chosen stress unit. Construct as `yeoh_t(c10,c20,c30,K)`.
        real(wp) :: c10=0.0_wp
            !! First modified-invariant term, positive, in the chosen stress unit.
        real(wp) :: c20=0.0_wp
            !! Second modified-invariant term, nonnegative, in the chosen stress unit.
        real(wp) :: c30=0.0_wp
            !! Third modified-invariant term, nonnegative, in the chosen stress unit.
        real(wp) :: bulk_modulus=0.0_wp
            !! Quadratic volumetric penalty modulus \(K>0\), in the chosen stress unit.
    contains
        procedure :: describe => describe_yeoh
            !! Declare input and output fields, parameter names and constitutive capabilities.
        procedure :: eval_batch => evaluate_yeoh
            !! Evaluate the batched constitutive response and optional consistent tangent.
        procedure :: eval_sensitivity_batch => sensitivity_yeoh
            !! Evaluate exact scalar-parameter derivatives of the constitutive response.
    end type yeoh_t

contains

    subroutine describe_yeoh(self,spec,status)
        !! Validate constitutive parameters and register `deformation.f` input and `stress.piola` response fields.
        class(yeoh_t),intent(in) :: self
            !! Constitutive model instance containing the fixed material parameters.
        type(material_spec_t),intent(out) :: spec
            !! Packed field schema, material parameter metadata and declared capabilities.
        integer,intent(out) :: status
            !! Operation status; `cf_ok` indicates success.

        status=cf_invalid_argument
        if (.not. all(ieee_is_finite([self%c10,self%c20,self%c30,self%bulk_modulus]))) return
        if (self%c10<=0.0_wp.or.self%c20<0.0_wp.or.self%c30<0.0_wp.or.self%bulk_modulus<=0.0_wp) return
        call spec%add_input("deformation.f", nfull, "tensor_full", "deformation_gradient", "spatial_reference", [ndim,ndim])
        call spec%add_output("stress.piola", nfull, "tensor_full", "first_piola_kirchhoff", "spatial_reference", [ndim,ndim])
        call spec%add_parameter("c10")
        call spec%add_parameter("c20")
        call spec%add_parameter("c30")
        call spec%add_parameter("bulk_modulus")
        call spec%set_has_tangent(.true.)
        call spec%set_has_sensitivity(.true.)
        status=cf_ok
    end subroutine describe_yeoh

    subroutine evaluate_yeoh(self,x_old,x_new,state_old,state_trial,response,dt,tangent,point_status)
        !! Evaluate `stress.piola` from `deformation.f` over the current point batch.
        class(yeoh_t),intent(in) :: self
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
    end subroutine evaluate_yeoh

    subroutine sensitivity_yeoh(self,x_old,x_new,state_old,state_trial,selected_parameters,sensitivity_old,sensitivity_trial,response,dsdp, dt,tangent,point_status)
        !! Evaluate exact derivatives of `stress.piola` with respect to the selected scalar parameters.
        class(yeoh_t),intent(in) :: self
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
    end subroutine sensitivity_yeoh

    subroutine evaluate_points(self,x,response,point_status,tangent,selected_parameters,dsdp)
        !! Compute `stress.piola` and optional exact input and parameter derivatives for each point.
        class(yeoh_t),intent(in) :: self
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
        real(wp) :: cof(ndim,ndim)
        real(wp) :: a(ndim,ndim)
        real(wp) :: q(ndim,ndim)
        real(wp) :: p(ndim,ndim)
        real(wp) :: base(ndim,ndim)
        real(wp) :: jac
        real(wp) :: jm23
        real(wp) :: i1
        real(wp) :: s
        real(wp) :: w
        real(wp) :: wprime
        real(wp) :: delta
        real(wp) :: dq
        real(wp) :: val
        integer :: ipoint
        integer :: i
        integer :: j
        integer :: k
        integer :: l
        integer :: row
        integer :: col
        integer :: iparam

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
            jm23=jac**(-2.0_wp/3.0_wp)
            i1=sum(f*f)
            s=jm23*i1-3.0_wp
            w=self%c10+2.0_wp*self%c20*s+3.0_wp*self%c30*s**2
            wprime=2.0_wp*self%c20+6.0_wp*self%c30*s
            if (.not. all(ieee_is_finite([jm23,i1,s,w,wprime]))) cycle point_loop
            q=f-i1*a/3.0_wp
            base=2.0_wp*jm23*q
            p=w*base+self%bulk_modulus*(jac-1.0_wp)*jac*a
            if (.not. all(ieee_is_finite(p))) cycle point_loop

            if (present(tangent)) then
                do j=1,ndim
                    do i=1,ndim
                        row=i+ndim*(j-1)
                        do l=1,ndim
                            do k=1,ndim
                                col=k+ndim*(l-1)
                                delta=0.0_wp
                                if (i==k.and.j==l) delta=1.0_wp
                                dq=delta-(2.0_wp/3.0_wp)*f(k,l)*a(i,j)+ (i1/3.0_wp)*a(i,l)*a(k,j)
                                val=4.0_wp*wprime*jm23**2*q(i,j)*q(k,l)+ 2.0_wp*w*jm23*(dq-(2.0_wp/3.0_wp)*a(k,l)*q(i,j))+ self%bulk_modulus*jac*((2.0_wp*jac-1.0_wp)*a(k,l)*a(i,j)- (jac-1.0_wp)*a(i,l)*a(k,j))
                                tangent(row,col,ipoint)=val
                            end do
                        end do
                    end do
                end do
                if (.not. all(ieee_is_finite(tangent(:,:,ipoint)))) cycle point_loop
            end if

            if (present(selected_parameters)) then
                do iparam=1,size(selected_parameters)
                    select case(selected_parameters(iparam))
                    case(1)
                        dsdp(:,iparam,ipoint)=reshape(base,[nfull])
                    case(2)
                        dsdp(:,iparam,ipoint)=reshape(2.0_wp*s*base,[nfull])
                    case(3)
                        dsdp(:,iparam,ipoint)=reshape(3.0_wp*s*s*base,[nfull])
                    case(4)
                        dsdp(:,iparam,ipoint)=reshape((jac-1.0_wp)*jac*a,[nfull])
                    case default
                        cycle point_loop
                    end select
                end do
                if (.not. all(ieee_is_finite(dsdp(:,:,ipoint)))) cycle point_loop
            end if
            response(:,ipoint)=reshape(p,[nfull])
            point_status(ipoint)=cf_ok
        end do point_loop
    end subroutine evaluate_points

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

end module constituf_yeoh

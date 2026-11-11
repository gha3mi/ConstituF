! SPDX-License-Identifier: Apache-2.0

program test_ogden
    !! Test convention: stress and stress-like material parameters in Pa.
    !! Strain and deformation gradients are dimensionless. Explicit rescaling
    !! cases state their factors.
    !!
    !! Spectral identity, Ogden alpha=2 limit, exact Hessian/sensitivities.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
    use constituf, only: wp,material_t,verify_tangent,cf_ok,cf_invalid_argument,cf_evaluation_failed
    use constituf_ogden, only: ogden_t
    use constituf_yeoh, only: yeoh_t
    implicit none (type, external)

    integer, parameter :: nfull=9
    type(material_t) :: mat,plain,yeoh,plus_mat,minus_mat
    type(ogden_t) :: pos_model,neg_model
    real(wp) :: x(nfull,2),y(nfull,2),t(nfull,nfull,2),dp(nfull,3,2),f(3,3)
    real(wp) :: yp(nfull,2),ym(nfull,2),tp(nfull,nfull,2),xpert(nfull,2)
    real(wp) :: numerical(nfull,2),coarse(nfull,2),tolerance(nfull,2)
    real(wp) :: expected(nfull)
    real(wp) :: h,error
    real(wp), parameter :: reference_tolerance=512.0_wp*epsilon(1.0_wp)
    real(wp), parameter :: parameter_values(3)=[120.0_wp,1.4_wp,6.0_wp]
    integer :: status,i,j,k,point,refinement,refinements

    f=0.0_wp
    f(1,1)=1.10_wp
    f(2,2)=0.95_wp
    f(3,3)=1.02_wp
    f(1,2)=0.20_wp
    f(2,1)=0.02_wp
    f(2,3)=0.07_wp
    x(:,1)=reshape(f,[9])
    f=0.0_wp
    f(1,1)=1.0_wp
    f(2,2)=1.0_wp
    f(3,3)=1.0_wp
    f(1,2)=0.12_wp
    x(:,2)=reshape(f,[9])

    call mat%init(ogden_t(6.0_wp,1.4_wp,120.0_wp),2,status, &
        sensitivity_parameters=[character(len=13) :: "bulk_modulus","exponent","shear_modulus"])
    call expect(status,cf_ok)
    if (mat%input_size()/=9.or.mat%output_size()/=9) error stop "Ogden schema"
    if (mat%state_size()/=0.or.mat%state_storage_bytes()/=0) error stop "Ogden stateless"
    call mat%eval(x,y,t,status,dsdp=dp)
    call expect(status,cf_ok)
    do point=1,2
        do j=1,9
            do i=1,9
                call close(t(i,j,point),t(j,i,point),reference_tolerance)
            end do
        end do
    end do

    call plain%init(ogden_t(6.0_wp,1.4_wp,120.0_wp),2,status)
    call expect(status,cf_ok)
    call verify_tangent(plain,x,error,status)
    call expect(status,cf_ok)
    if (.not. ieee_is_finite(error).or.error>1.0_wp) error stop "Ogden tangent mismatch"

    ! a=2 recovers the isochoric first-invariant Yeoh specialization.
    call plain%init(ogden_t(6.0_wp,2.0_wp,120.0_wp),2,status)
    call expect(status,cf_ok)
    call yeoh%init(yeoh_t(3.0_wp,0.0_wp,0.0_wp,120.0_wp),2,status)
    call expect(status,cf_ok)
    call plain%eval(x,yp,tp,status)
    call expect(status,cf_ok)
    call yeoh%eval(x,ym,t,status)
    call expect(status,cf_ok)
    do k=1,2
        do i=1,9
            call close(yp(i,k),ym(i,k),reference_tolerance)
            do j=1,9
                call close(tp(i,j,k),t(i,j,k),reference_tolerance)
            end do
        end do
    end do

    do k=1,3
        ! Moduli enter affinely; only the exponent needs a truncation estimate.
        h=0.25_wp*parameter_values(k)
        refinements=1
        if (k==2) then
            h=epsilon(1.0_wp)**(1.0_wp/3.0_wp)*parameter_values(k)
            refinements=2
        end if
        do refinement=1,refinements
            pos_model=ogden_t(6.0_wp,1.4_wp,120.0_wp)
            neg_model=pos_model
            select case(k)
            case(1)
                pos_model%bulk_modulus=pos_model%bulk_modulus+h
                neg_model%bulk_modulus=neg_model%bulk_modulus-h
            case(2)
                pos_model%exponent=pos_model%exponent+h
                neg_model%exponent=neg_model%exponent-h
            case(3)
                pos_model%shear_modulus=pos_model%shear_modulus+h
                neg_model%shear_modulus=neg_model%shear_modulus-h
            case default
                error stop "Unexpected parameter index"
            end select
            call plus_mat%init(pos_model,2,status)
            call expect(status,cf_ok)
            call minus_mat%init(neg_model,2,status)
            call expect(status,cf_ok)
            call plus_mat%eval(x,yp,status=status)
            call expect(status,cf_ok)
            call minus_mat%eval(x,ym,status=status)
            call expect(status,cf_ok)
            numerical=(yp-ym)/(2.0_wp*h)
            if (refinement==1) coarse=numerical
            if (refinement<refinements) h=0.5_wp*h
        end do
        tolerance=64.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(yp),abs(ym))/h
        if (k==2) tolerance=tolerance+abs(numerical-coarse)/3.0_wp
        if (.not. all(ieee_is_finite(numerical)).or..not. all(ieee_is_finite(dp(:,k,:)))) then
            error stop "Nonfinite Ogden parameter derivative"
        end if
        if (any(abs(numerical-dp(:,k,:))>tolerance)) error stop "Ogden parameter derivative"
    end do

    call test_repeated_stretches()
    call test_spectral_resolution()
    call plain%init(ogden_t(6.0_wp,0.1_wp,120.0_wp),2,status)
    call expect(status,cf_invalid_argument)
    call plain%init(ogden_t(6.0_wp,-2.0_wp,120.0_wp),2,status, &
        sensitivity_parameters=[character(len=13) :: "bulk_modulus","exponent","shear_modulus"])
    call expect(status,cf_ok)
    xpert=x
    xpert(:,1)=0.0_wp
    xpert([1,5,9],1)=[2.0_wp,1.0_wp,0.5_wp]
    call plain%eval(xpert,yp,status=status,dsdp=dp)
    call expect(status,cf_ok)
    ! J=1 and lambda^alpha=[1/4,1,4] give an independent negative-exponent reference.
    expected=0.0_wp
    expected([1,5,9])=[4.5_wp,4.5_wp,-27.0_wp]
    do i=1,nfull
        call close(yp(i,1),expected(i),reference_tolerance)
        call close(dp(i,1,1),0.0_wp,reference_tolerance)
        call close(dp(i,3,1),expected(i)/6.0_wp,reference_tolerance)
    end do
    expected([1,5,9])=[2.25_wp-4.5_wp*log(2.0_wp),2.25_wp-7.5_wp*log(2.0_wp), &
        33.0_wp*log(2.0_wp)-13.5_wp]
    do i=1,nfull
        call close(dp(i,2,1),expected(i),reference_tolerance)
    end do
    xpert=x
    xpert(1,1)=ieee_value(0.0_wp,ieee_quiet_nan)
    call plain%eval(xpert,yp,status=status,dsdp=dp)
    call expect(status,cf_evaluation_failed)

contains

    subroutine expect(actual,expected)
        integer,intent(in) :: actual,expected
        if (actual/=expected) then
            print *, "Ogden status",actual,expected
            error stop "Ogden status mismatch"
        end if
    end subroutine expect

    subroutine close(actual,expected,relative)
        real(wp),intent(in) :: actual,expected,relative
        if (.not. ieee_is_finite(actual).or..not. ieee_is_finite(expected)) then
            error stop "Nonfinite Ogden reference comparison"
        end if
        if (abs(actual-expected)>relative*max(1.0_wp,abs(expected))) then
            print *, "Ogden values",actual,expected
            error stop "Ogden reference mismatch"
        end if
    end subroutine close

    subroutine test_repeated_stretches()
        type(material_t) :: repeated
        real(wp) :: driving(nfull,2),response(nfull,2),jacobian(nfull,nfull,2)
        real(wp), parameter :: exponents(2)=[1.4_wp,-2.0_wp]
        integer :: result,exponent_index,component

        driving=0.0_wp
        driving([1,5,9],1)=1.0_wp
        driving([1,5,9],2)=1.25_wp
        do exponent_index=1,size(exponents)
            call repeated%init(ogden_t(6.0_wp,exponents(exponent_index),120.0_wp),2,result)
            call expect(result,cf_ok)
            call repeated%eval(driving,response,jacobian,result)
            call expect(result,cf_ok)
            do component=1,nfull
                call close(response(component,1),0.0_wp,reference_tolerance)
            end do
            call close(response(1,2),22875.0_wp/128.0_wp,reference_tolerance)
            call close(response(5,2),22875.0_wp/128.0_wp,reference_tolerance)
            call close(response(9,2),22875.0_wp/128.0_wp,reference_tolerance)
            ! Triple principal stretches: the initial shear modulus is mu for either exponent.
            call close(jacobian(1,1,1),128.0_wp,reference_tolerance)
            call close(jacobian(1,5,1),116.0_wp,reference_tolerance)
            call close(jacobian(4,4,1),6.0_wp,reference_tolerance)
            call close(jacobian(4,2,1),6.0_wp,reference_tolerance)
            ! At F=aI the isochoric Hessian scales by 1/a^2; volume adds prestress terms.
            call close(jacobian(1,1,2),128.0_wp/25.0_wp+9375.0_wp/32.0_wp,reference_tolerance)
            call close(jacobian(1,5,2),-64.0_wp/25.0_wp+6975.0_wp/16.0_wp,reference_tolerance)
            call close(jacobian(4,4,2),96.0_wp/25.0_wp,reference_tolerance)
            call close(jacobian(4,2,2),96.0_wp/25.0_wp-4575.0_wp/32.0_wp,reference_tolerance)
        end do
    end subroutine test_repeated_stretches

    subroutine test_spectral_resolution()
        type(material_t) :: spectral
        real(wp) :: driving(nfull,1),response(nfull,1),jacobian(nfull,nfull,1)
        real(wp) :: gamma,c1,c2,determinant,weight,reference
        real(wp), parameter :: alpha=2.0e6_wp
        integer :: result

        ! At alpha=4 the shear law is the polynomial P12=mu*gamma*(1+gamma^2/2).
        gamma=32.0_wp*epsilon(1.0_wp)
        driving=0.0_wp
        driving([1,5,9],1)=1.0_wp
        driving(4,1)=gamma
        call spectral%init(ogden_t(1.0_wp,4.0_wp,10.0_wp),1,result)
        call expect(result,cf_ok)
        call spectral%eval(driving,response,status=result)
        call expect(result,cf_ok)
        ! Spectral reconstruction subtracts O(1) terms: retain its O(eps/gamma) error bound.
        call close(response(4,1)/gamma,1.0_wp+0.5_wp*gamma**2,2.0_wp*epsilon(1.0_wp)/gamma)

        ! Distinct nearby eigenvalues need their divided difference, even at large alpha.
        driving=0.0_wp
        driving(1,1)=1.0_wp
        driving(5,1)=sqrt(1.0_wp+8.0e-7_wp)
        driving(9,1)=1.0_wp/driving(5,1)
        call spectral%init(ogden_t(1.0_wp,alpha,10.0_wp),1,result)
        call expect(result,cf_ok)
        call spectral%eval(driving,response,jacobian,result)
        call expect(result,cf_ok)
        c1=driving(1,1)**2
        c2=driving(5,1)**2
        determinant=driving(1,1)*driving(5,1)*driving(9,1)
        weight=(2.0_wp/alpha)*determinant**(-alpha/3.0_wp)
        reference=weight*(c1**(alpha/2.0_wp)-c2**(alpha/2.0_wp))/(c1-c2)
        call close(jacobian(4,4,1),reference,reference_tolerance)
    end subroutine test_spectral_resolution

end program test_ogden

! SPDX-License-Identifier: Apache-2.0

program test_yeoh
    !! Test convention: stress and stress-like material parameters in Pa.
    !! Strain and deformation gradients are dimensionless. Explicit rescaling
    !! cases state their factors.
    !!
    !! Independent simple-shear reference, finite-difference audit, API paths.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use constituf, only: wp,material_t,verify_tangent,cf_ok,cf_invalid_argument,cf_evaluation_failed,cf_no_trial
    use constituf_yeoh, only: yeoh_t
    implicit none (type, external)

    integer, parameter :: nfull=9
    type(material_t) :: mat,plus_mat,minus_mat,plain
    type(yeoh_t) :: plus_model,minus_model
    real(wp) :: x(nfull,2),y(nfull,2),t(nfull,nfull,2),dp(nfull,4,2),xp(nfull,2),yp(nfull,2),ym(nfull,2)
    real(wp) :: base(4),expected(nfull)
    real(wp) :: h,ref,gamma,numerical,tolerance,error
    real(wp), parameter :: reference_tolerance=256.0_wp*epsilon(1.0_wp)
    integer :: status,k,i,point

    base=[4.0_wp,1.5_wp,0.3_wp,100.0_wp]
    call mat%init(yeoh_t(base(1),base(2),base(3),base(4)),2,status, &
        sensitivity_parameters=[character(len=12) :: "bulk_modulus","c30","c20","c10"])
    call expect(status,cf_ok)
    if (mat%input_size()/=9.or.mat%output_size()/=9) error stop "Yeoh shape"
    if (mat%state_size()/=0.or.mat%state_storage_bytes()/=0) error stop "Yeoh is stateless"
    if (.not. mat%field_matches("input","deformation.f","tensor_full", &
        "deformation_gradient","spatial_reference",[3,3])) error stop "Yeoh metadata"
    if (mat%sensitivity_size()/=4) error stop "Yeoh parameter selection"

    gamma=0.25_wp
    x=0.0_wp
    x(1,:)=1.0_wp
    x(5,:)=1.0_wp
    x(9,:)=1.0_wp
    x(4,1)=gamma
    x(2,2)=-0.1_wp
    x(5,2)=0.95_wp
    call mat%eval(x,y,t,status,dsdp=dp)
    call expect(status,cf_ok)
    ref=2.0_wp*(base(1)+2.0_wp*base(2)*gamma**2+ &
        3.0_wp*base(3)*gamma**4)*gamma
    ! The complete simple-shear Piola stress includes normal stresses and P21 /= P12.
    expected=0.0_wp
    expected([1,5,9])=-gamma*ref/3.0_wp
    expected(2)=ref*(1.0_wp+gamma**2/3.0_wp)
    expected(4)=ref
    do i=1,nfull
        call close(y(i,1),expected(i),reference_tolerance)
    end do
    ref=2.0_wp*base(1)+12.0_wp*base(2)*gamma**2+30.0_wp*base(3)*gamma**4
    call close(t(4,4,1),ref,reference_tolerance)

    ! Check both the isochoric shear and the second point with J = 0.95.
    call plain%init(yeoh_t(base(1),base(2),base(3),base(4)),2,status)
    call expect(status,cf_ok)
    call verify_tangent(plain,x,error,status)
    call expect(status,cf_ok)
    if (.not. ieee_is_finite(error).or.error>1.0_wp) error stop "Yeoh tangent mismatch"

    do k=1,4
        plus_model=yeoh_t(base(1),base(2),base(3),base(4))
        minus_model=plus_model
        ! Stress is affine in all four parameters: these secants have no truncation error.
        h=0.25_wp*base(5-k)
        select case(k)
        case(1)
            plus_model%bulk_modulus=plus_model%bulk_modulus+h
            minus_model%bulk_modulus=minus_model%bulk_modulus-h
        case(2)
            plus_model%c30=plus_model%c30+h
            minus_model%c30=minus_model%c30-h
        case(3)
            plus_model%c20=plus_model%c20+h
            minus_model%c20=minus_model%c20-h
        case(4)
            plus_model%c10=plus_model%c10+h
            minus_model%c10=minus_model%c10-h
        case default
            error stop "Unexpected parameter index"
        end select
        call plus_mat%init(plus_model,2,status)
        call expect(status,cf_ok)
        call minus_mat%init(minus_model,2,status)
        call expect(status,cf_ok)
        call plus_mat%eval(x,yp,status=status)
        call expect(status,cf_ok)
        call minus_mat%eval(x,ym,status=status)
        call expect(status,cf_ok)
        do point=1,2
            do i=1,nfull
                numerical=(yp(i,point)-ym(i,point))/(2.0_wp*h)
                tolerance=64.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(yp(i,point)),abs(ym(i,point)))/h
                if (.not. ieee_is_finite(numerical).or..not. ieee_is_finite(dp(i,k,point))) then
                    error stop "Nonfinite Yeoh parameter derivative"
                end if
                if (abs(numerical-dp(i,k,point))>tolerance) error stop "Yeoh parameter derivative"
            end do
        end do
    end do

    ! Identity gives mu0 = 2*c10; dilation F = (5/4)I isolates the volume penalty.
    xp=0.0_wp
    xp([1,5,9],1)=1.0_wp
    xp([1,5,9],2)=1.25_wp
    call mat%eval(xp,yp,t,status,dsdp=dp)
    call expect(status,cf_ok)
    expected=0.0_wp
    expected([1,5,9])=1525.0_wp/1024.0_wp
    do i=1,nfull
        call close(yp(i,1),0.0_wp,reference_tolerance)
        call close(yp(i,2),base(4)*expected(i),reference_tolerance)
        call close(dp(i,1,2),expected(i),reference_tolerance)
        do k=2,4
            call close(dp(i,k,2),0.0_wp,reference_tolerance)
        end do
    end do
    call close(t(1,1,1),332.0_wp/3.0_wp,reference_tolerance)
    call close(t(1,5,1),284.0_wp/3.0_wp,reference_tolerance)
    call close(t(4,4,1),8.0_wp,reference_tolerance)
    call close(t(4,2,1),8.0_wp,reference_tolerance)

    call mat%commit(status)
    call expect(status,cf_ok)
    call mat%revert()
    call mat%commit(status)
    call expect(status,cf_no_trial)

    call plain%init(yeoh_t(-1.0_wp,0.0_wp,0.0_wp,100.0_wp),2,status)
    call expect(status,cf_invalid_argument)
    call plain%init(yeoh_t(base(1),base(2),base(3),base(4)),2,status)
    call expect(status,cf_ok)
    xp=x
    xp(1,2)=-1.0_wp
    call plain%eval(xp,yp,status=status)
    call expect(status,cf_evaluation_failed)
    xp=x
    xp(1,2)=ieee_value(0.0_wp,ieee_quiet_nan)
    call plain%eval(xp,yp,status=status)
    call expect(status,cf_evaluation_failed)

contains

    subroutine expect(actual,expected)
        integer,intent(in) :: actual,expected
        if (actual/=expected) then
            print *, "Yeoh status", actual, "expected", expected
            error stop "Yeoh status mismatch"
        end if
    end subroutine expect

    subroutine close(actual,expected,relative)
        real(wp),intent(in) :: actual,expected,relative
        if (.not. ieee_is_finite(actual).or..not. ieee_is_finite(expected)) then
            error stop "Nonfinite Yeoh reference comparison"
        end if
        if (abs(actual-expected)>relative*max(1.0_wp,abs(expected))) then
            print *, "Yeoh values",actual,expected
            error stop "Yeoh reference mismatch"
        end if
    end subroutine close

end program test_yeoh

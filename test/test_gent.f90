! SPDX-License-Identifier: Apache-2.0

program test_gent
    !! Test convention: stress and stress-like material parameters in Pa.
    !! Strain and deformation gradients are dimensionless. Explicit rescaling
    !! cases state their factors.
    !!
    !! Shear limit, analytic Hessian, objectivity, sensitivities and failed steps.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use constituf, only: wp, material_t, cf_ok, cf_invalid_argument, &
                        cf_evaluation_failed, cf_invalid_shape, cf_no_trial
    use constituf_gent, only: gent_t
    use constituf_mooney_rivlin, only: mooney_rivlin_t
    implicit none (type, external)

    integer, parameter :: ndim = 3, nfull = ndim*ndim

    call test_reference()
    call test_tangent_and_rotation()
    call test_parameter_derivatives()
    call test_failures()

contains

    subroutine expect(actual, wanted)
        integer, intent(in) :: actual, wanted

        if (actual /= wanted) then
            print *, "status", actual, "wanted", wanted
            error stop "Gent status mismatch"
        end if
    end subroutine expect

    subroutine close(actual, wanted)
        real(wp), intent(in) :: actual, wanted

        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(wanted)) then
            error stop "Nonfinite Gent comparison"
        end if
        if (abs(actual-wanted) > &
            1024.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(wanted))) then
            print *, "actual", actual, "wanted", wanted
            error stop "Gent reference mismatch"
        end if
    end subroutine close

    pure function identity() result(f)
        real(wp) :: f(ndim,ndim)
        integer :: i

        f = 0.0_wp
        do i = 1, ndim
            f(i,i) = 1.0_wp
        end do
    end function identity

    subroutine test_reference()
        type(material_t) :: material, limiting
        real(wp) :: f(ndim,ndim), x(nfull,1), p(nfull,1), p0(nfull,1)
        real(wp) :: tangent(nfull,nfull,1)
        real(wp) :: dsdp(nfull,1,1)
        real(wp) :: value
        integer :: status, i, j, k, l, row, column

        call material%init(gent_t(12.0_wp,4.0_wp,80.0_wp),1,status, &
            sensitivity_parameters=[character(len=18) :: "limiting_parameter"])
        call expect(status,cf_ok)
        if (material%input_size() /= nfull .or. material%output_size() /= nfull) then
            error stop "Gent tensor sizes"
        end if
        if (material%state_storage_bytes() /= 0) error stop "Nonzero Gent state"
        if (.not. material%field_matches("input", "deformation.f", &
            "tensor_full", "deformation_gradient", "spatial_reference", [ndim,ndim])) then
            error stop "Gent input metadata"
        end if
        if (.not. material%field_matches("output", "stress.piola", &
            "tensor_full", "first_piola_kirchhoff", "spatial_reference", [ndim,ndim])) then
            error stop "Gent output metadata"
        end if
        x(:,1) = reshape(identity(),[nfull])
        call material%eval(x,p,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        if (maxval(abs(p)) > 1.0e-12_wp) error stop "Gent zero stress"
        do j = 1, ndim
            do i = 1, ndim
                row = i+ndim*(j-1)
                do l = 1, ndim
                    do k = 1, ndim
                        column = k+ndim*(l-1)
                        value = 0.0_wp
                        if (i == j .and. k == l) value = 80.0_wp-8.0_wp
                        if (i == k .and. j == l) value = value+12.0_wp
                        if (i == l .and. j == k) value = value+12.0_wp
                        call close(tangent(row,column,1),value)
                    end do
                end do
            end do
        end do

        f = identity()
        f(1,2) = 1.0_wp
        x(:,1) = reshape(f,[nfull])
        call material%eval(x,p,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)
        ! Gamma = 1 gives slack = 3/4 and independently known locking derivatives.
        call close(p(4,1),16.0_wp)
        call close(tangent(4,4,1),80.0_wp/3.0_wp)
        call close(dsdp(4,1,1),-4.0_wp/3.0_wp)
        call material%commit(status)
        call expect(status,cf_ok)

        call material%init(gent_t(12.0_wp,1.0e9_wp,80.0_wp),1,status)
        call expect(status,cf_ok)
        call limiting%init(mooney_rivlin_t(6.0_wp,0.0_wp,80.0_wp),1,status)
        call expect(status,cf_ok)
        call material%eval(x,p,status=status)
        call expect(status,cf_ok)
        call limiting%eval(x,p0,status=status)
        call expect(status,cf_ok)
        if (maxval(abs(p-p0)) > 1.0e-7_wp*max(1.0_wp,maxval(abs(p0)))) then
            error stop "Gent large-Jm Mooney-Rivlin limit"
        end if
    end subroutine test_reference

    subroutine test_tangent_and_rotation()
        type(material_t) :: batch, single
        real(wp) :: f(ndim,ndim), rotation(ndim,ndim), stress(ndim,ndim)
        real(wp) :: x(nfull,2), p(nfull,2), tangent(nfull,nfull,2)
        real(wp) :: trial(nfull,2), plus(nfull,2), minus(nfull,2), alone(nfull,1)
        real(wp) :: expected(ndim,ndim)
        real(wp) :: coarse(nfull)
        real(wp) :: derivative, h, tolerance
        integer :: status, i, j, refinement

        f = identity()
        f(1,1) = 1.13_wp
        f(2,2) = 0.94_wp
        f(3,3) = 1.08_wp
        f(1,2) = 0.12_wp
        f(2,3) = -0.09_wp
        f(3,1) = 0.05_wp
        rotation = reshape([2.0_wp,1.0_wp,-2.0_wp,1.0_wp,2.0_wp,2.0_wp, &
                            2.0_wp,-2.0_wp,1.0_wp],[ndim,ndim])/3.0_wp
        x(:,1) = reshape(f,[nfull])
        x(:,2) = reshape(matmul(rotation,f),[nfull])
        call batch%init(gent_t(8.0_wp,5.0_wp,80.0_wp),2,status)
        call expect(status,cf_ok)
        call single%init(gent_t(8.0_wp,5.0_wp,80.0_wp),1,status)
        call expect(status,cf_ok)
        call batch%eval(x,p,tangent,status)
        call expect(status,cf_ok)
        call single%eval(x(:,1:1),alone,status=status)
        call expect(status,cf_ok)
        do i = 1, nfull
            call close(p(i,1),alone(i,1))
        end do
        stress = reshape(p(:,1),[ndim,ndim])
        expected = matmul(rotation,stress)
        do i = 1, nfull
            call close(p(i,2),expected(mod(i-1,ndim)+1,(i-1)/ndim+1))
        end do

        do j = 1, nfull
            h = epsilon(1.0_wp)**(1.0_wp/3.0_wp)*max(1.0_wp,abs(x(j,1)))
            do refinement = 1, 2
                trial = x
                trial(j,1) = x(j,1)+h
                call batch%eval(trial,plus,status=status)
                call expect(status,cf_ok)
                trial(j,1) = x(j,1)-h
                call batch%eval(trial,minus,status=status)
                call expect(status,cf_ok)
                if (refinement == 1) then
                    coarse = (plus(:,1)-minus(:,1))/(2.0_wp*h)
                    h = 0.5_wp*h
                end if
            end do
            do i = 1, nfull
                derivative = (plus(i,1)-minus(i,1))/(2.0_wp*h)
                tolerance = abs(derivative-coarse(i))/3.0_wp + &
                    64.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(plus(i,1)),abs(minus(i,1)))/h
                if (.not. ieee_is_finite(derivative) .or. .not. ieee_is_finite(coarse(i)) .or. &
                    abs(derivative-tangent(i,j,1)) > tolerance) then
                    error stop "Gent tangent finite difference"
                end if
                call close(tangent(i,j,1),tangent(j,i,1))
            end do
        end do
        call batch%revert()
        call batch%commit(status)
        call expect(status,cf_no_trial)
        call batch%eval(x,p,status=status)
        call expect(status,cf_ok)
        call batch%commit(status)
        call expect(status,cf_ok)
    end subroutine test_tangent_and_rotation

    subroutine test_parameter_derivatives()
        type(material_t) :: tracked, upper, lower
        real(wp) :: f(ndim,ndim), x(nfull,1), p(nfull,1), dsdp(nfull,3,1)
        real(wp) :: positive(nfull,1), negative(nfull,1)
        real(wp) :: coarse(nfull)
        real(wp) :: derivative, h, tolerance
        integer :: status, iparam, i, refinement, refinements

        f = identity()
        f(1,1) = 1.10_wp
        f(2,2) = 0.93_wp
        f(3,3) = 1.04_wp
        f(1,2) = 0.23_wp
        x(:,1) = reshape(f,[nfull])
        call tracked%init(gent_t(7.0_wp,8.0_wp,90.0_wp),1,status, &
            sensitivity_parameters=[character(len=18) :: &
            "bulk_modulus","limiting_parameter","shear_modulus"])
        call expect(status,cf_ok)
        if (tracked%sensitivity_size() /= 3) error stop "Gent derivative selection"
        if (tracked%state_storage_bytes() /= 0) error stop "Gent derivative history"
        call tracked%eval(x,p,status=status,dsdp=dsdp)
        call expect(status,cf_ok)

        do iparam = 1, 3
            refinements = 1
            select case (iparam)
            case (1)
                h = 0.25_wp*90.0_wp
            case (2)
                h = epsilon(1.0_wp)**(1.0_wp/3.0_wp)*8.0_wp
                refinements = 2
            case (3)
                h = 0.25_wp*7.0_wp
            case default
                error stop "Invalid Gent parameter index"
            end select
            ! K and mu enter affinely. Only Jm needs a truncation estimate.
            do refinement = 1, refinements
                select case (iparam)
                case (1)
                    call upper%init(gent_t(7.0_wp,8.0_wp,90.0_wp+h),1,status)
                    call expect(status,cf_ok)
                    call lower%init(gent_t(7.0_wp,8.0_wp,90.0_wp-h),1,status)
                case (2)
                    call upper%init(gent_t(7.0_wp,8.0_wp+h,90.0_wp),1,status)
                    call expect(status,cf_ok)
                    call lower%init(gent_t(7.0_wp,8.0_wp-h,90.0_wp),1,status)
                case (3)
                    call upper%init(gent_t(7.0_wp+h,8.0_wp,90.0_wp),1,status)
                    call expect(status,cf_ok)
                    call lower%init(gent_t(7.0_wp-h,8.0_wp,90.0_wp),1,status)
                case default
                    error stop "Invalid Gent parameter index"
                end select
                call expect(status,cf_ok)
                call upper%eval(x,positive,status=status)
                call expect(status,cf_ok)
                call lower%eval(x,negative,status=status)
                call expect(status,cf_ok)
                if (refinement == 1) then
                    coarse = (positive(:,1)-negative(:,1))/(2.0_wp*h)
                    if (refinements == 2) h = 0.5_wp*h
                end if
            end do
            do i = 1, nfull
                derivative = (positive(i,1)-negative(i,1))/(2.0_wp*h)
                tolerance = abs(derivative-coarse(i))/3.0_wp + &
                    64.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(positive(i,1)),abs(negative(i,1)))/h
                if (.not. ieee_is_finite(derivative) .or. .not. ieee_is_finite(coarse(i)) .or. &
                    .not. ieee_is_finite(dsdp(i,iparam,1)) .or. abs(derivative-dsdp(i,iparam,1)) > tolerance) then
                    error stop "Gent parameter derivative mismatch"
                end if
            end do
        end do
        call tracked%commit(status)
        call expect(status,cf_ok)
        call tracked%eval(x,p,status=status)
        call expect(status,cf_invalid_argument)
        call tracked%eval(x,p,status=status,dsdp=dsdp)
        call expect(status,cf_ok)
    end subroutine test_parameter_derivatives

    subroutine test_failures()
        type(material_t) :: material
        real(wp) :: x(nfull,2), p(nfull,2), wrong(nfull-1,2)
        integer :: status

        call material%init(gent_t(0.0_wp,4.0_wp,80.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(gent_t(8.0_wp,0.0_wp,80.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(gent_t(8.0_wp,4.0_wp,-1.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(gent_t(8.0_wp,4.0_wp,80.0_wp),2,status)
        call expect(status,cf_ok)
        x(:,1) = reshape(identity(),[nfull])
        x(:,2) = x(:,1)
        x(4,2) = 2.0_wp
        call material%eval(x,p,status=status)
        call expect(status,cf_evaluation_failed)
        if (material%last_failure_point() /= 2) error stop "Gent lock point"
        call material%commit(status)
        call expect(status,cf_no_trial)
        x(4,2) = 3.0_wp
        call material%eval(x,p,status=status)
        call expect(status,cf_evaluation_failed)
        x(:,2) = x(:,1)
        x(1,2) = -1.0_wp
        call material%eval(x,p,status=status)
        call expect(status,cf_evaluation_failed)
        x(1,2) = ieee_value(0.0_wp,ieee_quiet_nan)
        call material%eval(x,p,status=status)
        call expect(status,cf_evaluation_failed)
        x(:,2) = x(:,1)
        x(4,2) = 1.95_wp
        call material%eval(x,p,status=status)
        call expect(status,cf_ok)
        if (p(4,2) <= 8.0_wp*1.95_wp) error stop "Gent stiffening"
        call material%eval(wrong,p,status=status)
        call expect(status,cf_invalid_shape)
    end subroutine test_failures

end program test_gent

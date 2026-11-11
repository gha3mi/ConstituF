! SPDX-License-Identifier: Apache-2.0

program test_mooney_rivlin
    !! Test convention: stress and stress-like material parameters in Pa.
    !! Strain and deformation gradients are dimensionless. Explicit rescaling
    !! cases state their factors.
    !!
    !! Energy-derived reference values, objectivity, tangents and sensitivities.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
    use constituf, only: wp, material_t, cf_ok, cf_invalid_argument, &
                        cf_evaluation_failed, cf_invalid_shape, cf_no_trial
    use constituf_mooney_rivlin, only: mooney_rivlin_t
    implicit none (type, external)

    integer, parameter :: ndim = 3, nfull = ndim*ndim

    call test_reference_and_metadata()
    call test_tangent_and_objectivity()
    call test_parameter_sensitivities()
    call test_invalid_and_recovery()

contains

    subroutine expect(actual, expected)
        integer, intent(in) :: actual, expected

        if (actual /= expected) then
            print *, "Status", actual, "expected", expected
            error stop "Mooney-Rivlin status failure"
        end if
    end subroutine expect

    subroutine close(actual, expected)
        real(wp), intent(in) :: actual, expected

        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) then
            error stop "Nonfinite Mooney-Rivlin comparison"
        end if
        if (abs(actual-expected) > &
            1024.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(expected))) then
            print *, "Value", actual, "expected", expected
            error stop "Mooney-Rivlin exact reference failure"
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

    subroutine test_reference_and_metadata()
        type(material_t) :: material
        real(wp) :: x(nfull,1), stress(nfull,1), tangent(nfull,nfull,1)
        real(wp) :: f(ndim,ndim)
        real(wp) :: expected(nfull), dsdp(nfull,3,1)
        real(wp) :: shear, lambda, reference
        integer :: status, i, j, k, l, row, column

        call material%init(mooney_rivlin_t(7.0_wp,3.0_wp,90.0_wp),1,status)
        call expect(status,cf_ok)
        if (material%input_size() /= nfull .or. material%output_size() /= nfull) then
            error stop "Mooney-Rivlin full tensor sizes"
        end if
        if (material%state_storage_bytes() /= 0) error stop "Nonzero stateless history"
        if (.not. material%field_matches("input", "deformation.f", &
            "tensor_full", "deformation_gradient", "spatial_reference", [ndim,ndim])) then
            error stop "Incorrect deformation frame"
        end if
        if (.not. material%field_matches("output", "stress.piola", &
            "tensor_full", "first_piola_kirchhoff", "spatial_reference", [ndim,ndim])) then
            error stop "Incorrect stress frame"
        end if
        x(:,1) = reshape(identity(),[nfull])
        call material%eval(x,stress,tangent,status)
        call expect(status,cf_ok)
        if (maxval(abs(stress)) > 1.0e-12_wp) error stop "Reference stress nonzero"
        shear = 2.0_wp*(7.0_wp+3.0_wp)
        lambda = 90.0_wp-(2.0_wp/3.0_wp)*shear
        do j = 1, ndim
            do i = 1, ndim
                row = i+ndim*(j-1)
                do l = 1, ndim
                    do k = 1, ndim
                        column = k+ndim*(l-1)
                        reference = 0.0_wp
                        if (i == j .and. k == l) reference = reference+lambda
                        if (i == k .and. j == l) reference = reference+shear
                        if (i == l .and. j == k) reference = reference+shear
                        call close(tangent(row,column,1),reference)
                    end do
                end do
            end do
        end do

        ! Unit-volume simple shear also has Piola normal stresses from the isochoric split.
        f = identity()
        f(1,2) = 0.25_wp
        x(:,1) = reshape(f,[nfull])
        call material%eval(x,stress,status=status)
        call expect(status,cf_ok)
        expected = 0.0_wp
        expected([1,5]) = -13.0_wp/24.0_wp
        expected(9) = -1.0_wp/6.0_wp
        expected(4) = 5.0_wp
        expected(2) = 5.0_wp+13.0_wp/96.0_wp
        do i = 1, nfull
            call close(stress(i,1),expected(i))
        end do

        ! F = (5/4) I has no isochoric stress. Pii/K = (a**3-1)*a**2 = 1525/1024.
        call material%init(mooney_rivlin_t(7.0_wp,3.0_wp,90.0_wp),1,status, &
            sensitivity_parameters=[character(len=12) :: "c10","c01","bulk_modulus"])
        call expect(status,cf_ok)
        x(:,1) = reshape(1.25_wp*identity(),[nfull])
        call material%eval(x,stress,status=status,dsdp=dsdp)
        call expect(status,cf_ok)
        expected = 0.0_wp
        expected([1,5,9]) = 1525.0_wp/1024.0_wp
        do i = 1, nfull
            call close(stress(i,1),90.0_wp*expected(i))
            call close(dsdp(i,1,1),0.0_wp)
            call close(dsdp(i,2,1),0.0_wp)
            call close(dsdp(i,3,1),expected(i))
        end do
    end subroutine test_reference_and_metadata

    subroutine test_tangent_and_objectivity()
        type(material_t) :: material, single
        real(wp) :: f(ndim,ndim), rotation(ndim,ndim), p(ndim,ndim)
        real(wp) :: x(nfull,2), stress(nfull,2), tangent(nfull,nfull,2)
        real(wp) :: perturb(nfull,2), plus(nfull,2), minus(nfull,2)
        real(wp) :: single_stress(nfull,1), rotated(ndim,ndim)
        real(wp) :: coarse(nfull)
        real(wp) :: h, derivative, tolerance
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
        call material%init(mooney_rivlin_t(8.0_wp,5.0_wp,80.0_wp),2,status)
        call expect(status,cf_ok)
        call single%init(mooney_rivlin_t(8.0_wp,5.0_wp,80.0_wp),1,status)
        call expect(status,cf_ok)
        call material%eval(x,stress,tangent,status)
        call expect(status,cf_ok)
        call single%eval(x(:,1:1),single_stress,status=status)
        call expect(status,cf_ok)
        do i = 1, nfull
            call close(stress(i,1),single_stress(i,1))
        end do
        p = reshape(stress(:,1),[ndim,ndim])
        rotated = matmul(rotation,p)
        do i = 1, nfull
            call close(stress(i,2),rotated(mod(i-1,ndim)+1,(i-1)/ndim+1))
        end do

        do j = 1, nfull
            h = epsilon(1.0_wp)**(1.0_wp/3.0_wp)*max(1.0_wp,abs(x(j,1)))
            do refinement = 1, 2
                perturb = x
                perturb(j,1) = x(j,1)+h
                call material%eval(perturb,plus,status=status)
                call expect(status,cf_ok)
                perturb(j,1) = x(j,1)-h
                call material%eval(perturb,minus,status=status)
                call expect(status,cf_ok)
                if (refinement == 1) then
                    coarse = (plus(:,1)-minus(:,1))/(2.0_wp*h)
                    h = 0.5_wp*h
                end if
            end do
            do i = 1, nfull
                derivative = (plus(i,1)-minus(i,1))/(2.0_wp*h)
                ! The step pair estimates centered-difference truncation error.
                tolerance = abs(derivative-coarse(i))/3.0_wp + &
                    64.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(plus(i,1)),abs(minus(i,1)))/h
                if (.not. ieee_is_finite(derivative) .or. .not. ieee_is_finite(coarse(i)) .or. &
                    abs(derivative-tangent(i,j,1)) > tolerance) then
                    error stop "Inconsistent Mooney-Rivlin tangent"
                end if
                call close(tangent(i,j,1),tangent(j,i,1))
            end do
        end do
        call material%revert()
        call material%commit(status)
        call expect(status,cf_no_trial)
        call material%eval(x,stress,status=status)
        call expect(status,cf_ok)
        call material%commit(status)
        call expect(status,cf_ok)
    end subroutine test_tangent_and_objectivity

    subroutine test_parameter_sensitivities()
        type(material_t) :: tracked, plus_material, minus_material
        real(wp) :: f(ndim,ndim), x(nfull,1), stress(nfull,1)
        real(wp) :: dsdp(nfull,3,1), tangent(nfull,nfull,1)
        real(wp) :: plus(nfull,1), minus(nfull,1)
        real(wp) :: h, derivative, tolerance
        integer :: status, iparam, i

        f = identity()
        f(1,1) = 1.10_wp
        f(2,2) = 0.92_wp
        f(3,3) = 1.04_wp
        f(1,2) = 0.11_wp
        x(:,1) = reshape(f,[nfull])
        call tracked%init(mooney_rivlin_t(7.0_wp,3.0_wp,90.0_wp),1,status, &
            sensitivity_parameters=[character(len=12) :: "bulk_modulus","c01","c10"])
        call expect(status,cf_ok)
        if (tracked%sensitivity_size() /= 3) error stop "Missing sensitivities"
        if (tracked%state_storage_bytes() /= 0) error stop "Derivative history allocated"
        call tracked%eval(x,stress,tangent,status,dsdp=dsdp)
        call expect(status,cf_ok)

        do iparam = 1, 3
            ! At fixed F the stress is affine in each parameter.
            select case (iparam)
            case (1)
                h = 0.25_wp*90.0_wp
                call plus_material%init(mooney_rivlin_t(7.0_wp,3.0_wp,90.0_wp+h),1,status)
                call expect(status,cf_ok)
                call minus_material%init(mooney_rivlin_t(7.0_wp,3.0_wp,90.0_wp-h),1,status)
            case (2)
                h = 0.25_wp*3.0_wp
                call plus_material%init(mooney_rivlin_t(7.0_wp,3.0_wp+h,90.0_wp),1,status)
                call expect(status,cf_ok)
                call minus_material%init(mooney_rivlin_t(7.0_wp,3.0_wp-h,90.0_wp),1,status)
            case (3)
                h = 0.25_wp*7.0_wp
                call plus_material%init(mooney_rivlin_t(7.0_wp+h,3.0_wp,90.0_wp),1,status)
                call expect(status,cf_ok)
                call minus_material%init(mooney_rivlin_t(7.0_wp-h,3.0_wp,90.0_wp),1,status)
            case default
                error stop "Invalid test parameter index"
            end select
            call expect(status,cf_ok)
            call plus_material%eval(x,plus,status=status)
            call expect(status,cf_ok)
            call minus_material%eval(x,minus,status=status)
            call expect(status,cf_ok)
            do i = 1, nfull
                derivative = (plus(i,1)-minus(i,1))/(2.0_wp*h)
                tolerance = 64.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(plus(i,1)),abs(minus(i,1)))/h
                if (.not. ieee_is_finite(derivative) .or. .not. ieee_is_finite(dsdp(i,iparam,1)) .or. &
                    abs(derivative-dsdp(i,iparam,1)) > tolerance) then
                    error stop "Mooney-Rivlin parameter sensitivity mismatch"
                end if
            end do
        end do
        call tracked%commit(status)
        call expect(status,cf_ok)
        call tracked%eval(x,stress,status=status)
        call expect(status,cf_invalid_argument)
        call tracked%eval(x,stress,status=status,dsdp=dsdp)
        call expect(status,cf_ok)
    end subroutine test_parameter_sensitivities

    subroutine test_invalid_and_recovery()
        type(material_t) :: material
        real(wp) :: x(nfull,2), stress(nfull,2), wrong(nfull-1,2)
        integer :: status

        call material%init(mooney_rivlin_t(0.0_wp,1.0_wp,10.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(mooney_rivlin_t(1.0_wp,-1.0_wp,10.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(mooney_rivlin_t(1.0_wp,0.0_wp,0.0_wp),2,status)
        call expect(status,cf_invalid_argument)
        call material%init(mooney_rivlin_t(1.0_wp,0.0_wp,10.0_wp),2,status)
        call expect(status,cf_ok)
        x(:,1) = reshape(identity(),[nfull])
        x(:,2) = x(:,1)
        x(1,2) = -1.0_wp
        call material%eval(x,stress,status=status)
        call expect(status,cf_evaluation_failed)
        if (material%last_failure_point() /= 2) error stop "Bad point not flagged"
        call material%commit(status)
        call expect(status,cf_no_trial)
        x(1,2) = 0.0_wp
        call material%eval(x,stress,status=status)
        call expect(status,cf_evaluation_failed)
        x(1,2) = ieee_value(0.0_wp,ieee_quiet_nan)
        call material%eval(x,stress,status=status)
        call expect(status,cf_evaluation_failed)
        x(1,2) = 1.0_wp
        call material%eval(x,stress,status=status)
        call expect(status,cf_ok)
        if (maxval(abs(stress)) > 1.0e-12_wp) error stop "Trial recovery failed"
        call material%eval(wrong,stress,status=status)
        call expect(status,cf_invalid_shape)
    end subroutine test_invalid_and_recovery

end program test_mooney_rivlin

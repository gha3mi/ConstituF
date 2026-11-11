! SPDX-License-Identifier: Apache-2.0

program test_neo_hookean
    !! Test convention: stress and stress-like material parameters in Pa.
    !! Strain and deformation gradients are dimensionless. Explicit rescaling
    !! cases state their factors.
    !!
    !! Analytical hyperelastic references, objectivity and complete API tests.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use constituf, only: wp, material_t, cf_ok, cf_invalid_argument, &
                        cf_evaluation_failed, cf_invalid_shape, cf_no_trial
    use constituf_neo_hookean, only: neo_hookean_t
    implicit none (type, external)

    integer, parameter :: ndim = 3
    integer, parameter :: nfull = ndim*ndim

    call test_metadata_and_identity()
    call test_simple_deformations()
    call test_tangent_and_objectivity()
    call test_sensitivities()
    call test_failure_recovery()

contains

    subroutine expect(actual, expected)
        integer, intent(in) :: actual, expected

        if (actual /= expected) then
            print *, "Expected status", expected, "received", actual
            error stop "Neo-Hookean operation failed"
        end if
    end subroutine expect

    subroutine close(actual, expected)
        real(wp), intent(in) :: actual, expected

        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) then
            error stop "Nonfinite neo-Hookean comparison"
        end if
        if (abs(actual-expected) > &
            512.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(expected))) then
            print *, "Expected", expected, "received", actual
            error stop "Neo-Hookean reference mismatch"
        end if
    end subroutine close

    pure function identity() result(value)
        real(wp) :: value(ndim,ndim)
        integer :: i

        value = 0.0_wp
        do i = 1, ndim
            value(i,i) = 1.0_wp
        end do
    end function identity

    subroutine test_metadata_and_identity()
        type(material_t) :: material
        real(wp) :: x(nfull,1), y(nfull,1), tangent(nfull,nfull,1)
        real(wp) :: f(ndim,ndim)
        real(wp) :: reference, lambda, mu
        integer :: status, i, j, k, l, row, column

        lambda = 50.0_wp
        mu = 20.0_wp
        call material%init(neo_hookean_t(lambda,mu), 1, status)
        call expect(status, cf_ok)
        if (material%input_size() /= nfull .or. material%output_size() /= nfull) then
            error stop "Full tensor sizes"
        end if
        if (material%state_storage_bytes() /= 0) error stop "Stateless model memory"
        if (material%input_offset("deformation.f") /= 1) error stop "Input offset"
        if (material%output_offset("stress.piola") /= 1) error stop "Output offset"
        if (.not. material%field_matches("input", "deformation.f", &
            "tensor_full", "deformation_gradient", "spatial_reference", [ndim,ndim])) then
            error stop "Incorrect deformation metadata"
        end if
        if (.not. material%field_matches("output", "stress.piola", &
            "tensor_full", "first_piola_kirchhoff", "spatial_reference", [ndim,ndim])) then
            error stop "Incorrect first Piola metadata"
        end if
        if (material%field_matches("output", "stress.piola", &
            "mandel", "cauchy_stress", "spatial", [ndim,ndim])) then
            error stop "Stress convention incorrectly accepted"
        end if

        x(:,1) = reshape(identity(),[nfull])
        call material%eval(x, y, tangent, status)
        call expect(status, cf_ok)
        if (maxval(abs(y)) > 1.0e-14_wp) error stop "Nonzero reference stress"

        do j = 1, ndim
            do i = 1, ndim
                row = i+ndim*(j-1)
                do l = 1, ndim
                    do k = 1, ndim
                        column = k+ndim*(l-1)
                        reference = 0.0_wp
                        if (i == k .and. j == l) reference = reference+mu
                        if (i == j .and. k == l) reference = reference+lambda
                        if (i == l .and. k == j) reference = reference+mu
                        call close(tangent(row,column,1), reference)
                    end do
                end do
            end do
        end do

        ! Diagonal deformation has an independent closed-form Piola stress.
        f = 0.0_wp
        f(1,1) = 1.2_wp
        f(2,2) = 0.8_wp
        f(3,3) = 1.1_wp
        x(:,1) = reshape(f,[nfull])
        call material%eval(x, y, tangent, status)
        call expect(status, cf_ok)
        do i = 1, ndim
            reference = mu*(f(i,i)-1.0_wp/f(i,i)) + &
                lambda*log(f(1,1)*f(2,2)*f(3,3))/f(i,i)
            call close(y(i+ndim*(i-1),1), reference)
        end do
        do j = 1, ndim
            do i = 1, ndim
                if (i /= j) call close(y(i+ndim*(j-1),1), 0.0_wp)
            end do
        end do
    end subroutine test_metadata_and_identity

    subroutine test_simple_deformations()
        type(material_t) :: material
        real(wp) :: x(nfull,2), stress(nfull,2), expected(nfull,2)
        integer :: status, i, ipoint

        call material%init(neo_hookean_t(50.0_wp,20.0_wp),2,status)
        call expect(status,cf_ok)
        ! F = I + (1/4) e1 tensor e2 has J = 1 and P12 = P21 = 5 Pa.
        x(:,1) = reshape(identity(),[nfull])
        x(4,1) = 0.25_wp
        x(:,2) = reshape(1.25_wp*identity(),[nfull])
        expected = 0.0_wp
        expected(2,1) = 5.0_wp
        expected(4,1) = 5.0_wp
        ! Uniform dilation: Pii = mu*(a-1/a) + 3*lambda*log(a)/a.
        expected([1,5,9],2) = 9.0_wp+120.0_wp*log(1.25_wp)
        call material%eval(x,stress,status=status)
        call expect(status,cf_ok)
        do ipoint = 1, 2
            do i = 1, nfull
                call close(stress(i,ipoint),expected(i,ipoint))
            end do
        end do
    end subroutine test_simple_deformations

    subroutine test_tangent_and_objectivity()
        type(material_t) :: batch, single
        real(wp) :: f(ndim,ndim), rotation(ndim,ndim)
        real(wp) :: x(nfull,2), stress(nfull,2), tangent(nfull,nfull,2)
        real(wp) :: plus(nfull,2), minus(nfull,2), perturb(nfull,2)
        real(wp) :: single_stress(nfull,1), p(ndim,ndim), rotated_p(ndim,ndim)
        real(wp) :: coarse(nfull)
        real(wp) :: h, derivative, tolerance
        integer :: status, i, j, refinement

        call batch%init(neo_hookean_t(70.0_wp,30.0_wp), 2, status)
        call expect(status, cf_ok)
        call single%init(neo_hookean_t(70.0_wp,30.0_wp), 1, status)
        call expect(status, cf_ok)

        f = identity()
        f(1,1) = 1.11_wp
        f(2,2) = 0.92_wp
        f(3,3) = 1.04_wp
        f(1,2) = 0.12_wp
        f(2,3) = -0.07_wp
        f(3,1) = 0.05_wp
        ! A proper rotation that mixes all three axes.
        rotation = reshape([2.0_wp,1.0_wp,-2.0_wp,1.0_wp,2.0_wp,2.0_wp, &
                            2.0_wp,-2.0_wp,1.0_wp],[ndim,ndim])/3.0_wp
        x(:,1) = reshape(f,[nfull])
        x(:,2) = reshape(matmul(rotation,f),[nfull])

        call batch%eval(x, stress, tangent, status)
        call expect(status, cf_ok)
        call single%eval(x(:,1:1), single_stress, status=status)
        call expect(status, cf_ok)
        do i = 1, nfull
            call close(stress(i,1), single_stress(i,1))
        end do
        p = reshape(stress(:,1),[ndim,ndim])
        rotated_p = matmul(rotation,p)
        do i = 1, nfull
            call close(stress(i,2), rotated_p(mod(i-1,ndim)+1,(i-1)/ndim+1))
        end do
        call single%eval(reshape(rotation,[nfull,1]),single_stress,status=status)
        call expect(status,cf_ok)
        do i = 1, nfull
            call close(single_stress(i,1),0.0_wp)
        end do

        ! Exact material tangent is the Hessian of the hyperelastic energy.
        do j = 1, nfull
            do i = 1, nfull
                call close(tangent(i,j,1), tangent(j,i,1))
            end do
        end do
        do j = 1, nfull
            h = epsilon(1.0_wp)**(1.0_wp/3.0_wp)*max(1.0_wp,abs(x(j,1)))
            do refinement = 1, 2
                perturb = x
                perturb(j,1) = x(j,1)+h
                call batch%eval(perturb,plus,status=status)
                call expect(status,cf_ok)
                perturb(j,1) = x(j,1)-h
                call batch%eval(perturb,minus,status=status)
                call expect(status,cf_ok)
                if (refinement == 1) then
                    coarse = (plus(:,1)-minus(:,1))/(2.0_wp*h)
                    h = 0.5_wp*h
                end if
            end do
            do i = 1, nfull
                derivative = (plus(i,1)-minus(i,1))/(2.0_wp*h)
                ! Estimate centered-difference truncation from the step pair.
                tolerance = abs(derivative-coarse(i))/3.0_wp + &
                    64.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(plus(i,1)),abs(minus(i,1)))/h
                if (.not. ieee_is_finite(derivative) .or. .not. ieee_is_finite(coarse(i)) .or. &
                    abs(derivative-tangent(i,j,1)) > tolerance) then
                    error stop "Finite-strain material tangent mismatch"
                end if
            end do
        end do
        call batch%revert()
        call batch%commit(status)
        call expect(status, cf_no_trial)
        call batch%eval(x, stress, status=status)
        call expect(status, cf_ok)
        call batch%commit(status)
        call expect(status, cf_ok)
        call batch%eval(x, stress, status=status)
        call expect(status, cf_ok)
    end subroutine test_tangent_and_objectivity

    subroutine test_sensitivities()
        type(material_t) :: tracked, plus, minus
        real(wp) :: deformation(ndim,ndim), x(nfull,1), stress(nfull,1)
        real(wp) :: dsdp(nfull,2,1), plus_stress(nfull,1), minus_stress(nfull,1)
        real(wp) :: tangent(nfull,nfull,1)
        real(wp) :: h, derivative, tolerance
        integer :: status, i, iparam

        deformation = identity()
        deformation(1,1) = 1.1_wp
        deformation(2,2) = 0.87_wp
        deformation(3,3) = 1.04_wp
        deformation(1,2) = 0.15_wp
        x(:,1) = reshape(deformation,[nfull])
        call tracked%init(neo_hookean_t(60.0_wp,25.0_wp), 1, status, &
            sensitivity_parameters=[character(len=16) :: "shear_modulus", "lame_lambda"])
        call expect(status, cf_ok)
        if (tracked%state_storage_bytes() /= 0) error stop "Derivative state allocated"
        if (tracked%sensitivity_size() /= 2) error stop "Selected parameter count"
        call tracked%eval(x, stress, tangent, status, dsdp=dsdp)
        call expect(status, cf_ok)

        do iparam = 1, 2
            ! Piola stress is affine in each Lame modulus at fixed F.
            if (iparam == 1) then
                h = 0.25_wp*25.0_wp
                call plus%init(neo_hookean_t(60.0_wp,25.0_wp+h), 1, status)
                call expect(status, cf_ok)
                call minus%init(neo_hookean_t(60.0_wp,25.0_wp-h), 1, status)
            else
                h = 0.25_wp*60.0_wp
                call plus%init(neo_hookean_t(60.0_wp+h,25.0_wp), 1, status)
                call expect(status, cf_ok)
                call minus%init(neo_hookean_t(60.0_wp-h,25.0_wp), 1, status)
            end if
            call expect(status, cf_ok)
            call plus%eval(x, plus_stress, status=status)
            call expect(status, cf_ok)
            call minus%eval(x, minus_stress, status=status)
            call expect(status, cf_ok)
            do i = 1, nfull
                derivative = (plus_stress(i,1)-minus_stress(i,1))/(2.0_wp*h)
                tolerance = 64.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(plus_stress(i,1)),abs(minus_stress(i,1)))/h
                if (.not. ieee_is_finite(derivative) .or. .not. ieee_is_finite(dsdp(i,iparam,1)) .or. &
                    abs(derivative-dsdp(i,iparam,1)) > tolerance) then
                    error stop "Neo-Hookean analytical sensitivity mismatch"
                end if
            end do
        end do
        call tracked%commit(status)
        call expect(status, cf_ok)
        call tracked%eval(x, stress, status=status)
        call expect(status, cf_invalid_argument)
        call tracked%eval(x, stress, status=status, dsdp=dsdp)
        call expect(status, cf_ok)
    end subroutine test_sensitivities

    subroutine test_failure_recovery()
        type(material_t) :: material
        real(wp) :: x(nfull,2), y(nfull,2), bad_input(nfull-1,2)
        integer :: status

        call material%init(neo_hookean_t(-1.0_wp,20.0_wp), 2, status)
        call expect(status, cf_invalid_argument)
        call material%init(neo_hookean_t(10.0_wp,0.0_wp), 2, status)
        call expect(status, cf_invalid_argument)
        call material%init(neo_hookean_t(0.0_wp,20.0_wp), 2, status)
        call expect(status, cf_ok)

        x(:,1) = reshape(identity(),[nfull])
        x(:,2) = reshape(identity(),[nfull])
        x(1,2) = -1.0_wp
        call material%eval(x, y, status=status)
        call expect(status, cf_evaluation_failed)
        if (material%last_failure_point() /= 2) error stop "Inverted deformation accepted"
        call material%commit(status)
        call expect(status, cf_no_trial)
        x(1,2) = 0.0_wp
        call material%eval(x, y, status=status)
        call expect(status, cf_evaluation_failed)
        x(1,2) = ieee_value(0.0_wp,ieee_quiet_nan)
        call material%eval(x, y, status=status)
        call expect(status, cf_evaluation_failed)
        x(1,2) = 1.0_wp
        call material%eval(x, y, status=status)
        call expect(status, cf_ok)
        if (maxval(abs(y)) > 1.0e-14_wp) error stop "Failed trial corrupted material"

        bad_input = 0.0_wp
        call material%eval(bad_input, y, status=status)
        call expect(status, cf_invalid_shape)
        call material%init(neo_hookean_t(-1.0_wp,20.0_wp), 2, status)
        call expect(status, cf_invalid_argument)
        call material%eval(x, y, status=status)
        call expect(status, cf_ok)
    end subroutine test_failure_recovery
end program test_neo_hookean

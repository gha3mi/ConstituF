! SPDX-License-Identifier: Apache-2.0

program test_orthotropic_elastic
    !! Test convention: stress and stress-like material parameters in Pa.
    !! Strain and deformation gradients are dimensionless. Explicit rescaling
    !! cases state their factors.
    !!
    !! Reference compliance, Mandel stiffness and all nine parameter derivatives.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use constituf, only: wp, material_t, cf_ok, cf_invalid_argument, &
                        cf_invalid_shape, cf_evaluation_failed, cf_no_trial
    use constituf_orthotropic_elastic, only: orthotropic_elastic_t
    implicit none (type, external)

    integer, parameter :: ndim = 3, nmandel = 6, nparameters = 9
    integer, parameter :: npoints = 2

    call test_uniaxial_compliance()
    call test_isotropic_limit_and_tangent()
    call test_all_sensitivities()
    call test_modulus_scales()
    call test_invalid_and_trial()

contains

    function reference_model() result(model)
        type(orthotropic_elastic_t) :: model

        model = orthotropic_elastic_t(100.0_wp,80.0_wp,60.0_wp, &
                                     0.20_wp,0.13_wp,0.12_wp, &
                                     25.0_wp,30.0_wp,35.0_wp)
    end function reference_model

    subroutine expect_status(actual, expected)
        integer, intent(in) :: actual, expected

        if (actual /= expected) then
            print *, "Expected status", expected, "got", actual
            error stop "Orthotropic API status mismatch"
        end if
    end subroutine expect_status

    subroutine close(actual, expected)
        real(wp), intent(in) :: actual, expected
        real(wp) :: tolerance

        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) then
            error stop "Nonfinite orthotropic comparison"
        end if
        tolerance = 512.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(expected))
        if (abs(actual-expected) > tolerance) then
            print *, "Expected", expected, "got", actual
            error stop "Orthotropic numerical mismatch"
        end if
    end subroutine close

    subroutine test_uniaxial_compliance()
        type(material_t) :: mat
        type(orthotropic_elastic_t) :: model
        real(wp) :: x(nmandel,npoints), y(nmandel,npoints)
        real(wp) :: tangent(nmandel,nmandel,npoints)
        integer :: status, i, j

        model = reference_model()
        call mat%init(model,npoints,status)
        call expect_status(status,cf_ok)
        if (mat%state_size() /= 0 .or. mat%state_storage_bytes() /= 0) then
            error stop "Orthotropic elasticity must have no history"
        end if
        if (mat%input_size() /= nmandel .or. mat%output_size() /= nmandel) then
            error stop "Mandel field dimension mismatch"
        end if
        if (.not. mat%field_matches("input","strain.small","mandel", &
                                    "small_strain","spatial",[ndim,ndim])) then
            error stop "Orthotropic strain metadata mismatch"
        end if
        if (.not. mat%field_matches("output","stress.cauchy","mandel", &
                                    "cauchy_stress","spatial",[ndim,ndim])) then
            error stop "Orthotropic stress metadata mismatch"
        end if

        ! A unit axial stress in axis 1 produces the known compliance column.
        x = 0.0_wp
        x(1,1) = 1.0_wp/model%young_1
        x(2,1) = -model%poisson_12/model%young_1
        x(3,1) = -model%poisson_13/model%young_1
        ! A unit axial stress in axis 2 gives the second compliance column.
        x(1,2) = -model%poisson_12/model%young_1
        x(2,2) = 1.0_wp/model%young_2
        x(3,2) = -model%poisson_23/model%young_2

        call mat%eval(x,y,tangent,status)
        call expect_status(status,cf_ok)
        call close(y(1,1),1.0_wp)
        call close(y(2,2),1.0_wp)
        do i = 1, nmandel
            if (i /= 1) call close(y(i,1),0.0_wp)
            if (i /= 2) call close(y(i,2),0.0_wp)
        end do
        do i = 1, nmandel
            do j = 1, nmandel
                call close(tangent(i,j,1),tangent(j,i,1))
            end do
        end do
        if (dot_product(x(:,1),y(:,1)) <= 0.0_wp) then
            error stop "Nonpositive orthotropic elastic energy"
        end if

        ! Complete the reciprocal compliance columns and distinguish all shear planes.
        x = 0.0_wp
        x(1,1) = -model%poisson_13/model%young_1
        x(2,1) = -model%poisson_23/model%young_2
        x(3,1) = 1.0_wp/model%young_3
        x(4:6,2) = sqrt(2.0_wp)*[0.001_wp,0.002_wp,0.003_wp]
        call mat%eval(x,y,status=status)
        call expect_status(status,cf_ok)
        call close(y(3,1),1.0_wp)
        do i = 1, nmandel
            if (i /= 3) call close(y(i,1),0.0_wp)
            if (i <= ndim) call close(y(i,2),0.0_wp)
        end do
        call close(y(4,2),sqrt(2.0_wp)*0.05_wp)
        call close(y(5,2),sqrt(2.0_wp)*0.12_wp)
        call close(y(6,2),sqrt(2.0_wp)*0.21_wp)
        call mat%commit(status)
        call expect_status(status,cf_ok)
    end subroutine test_uniaxial_compliance

    subroutine test_isotropic_limit_and_tangent()
        type(material_t) :: mat, single
        real(wp) :: x(nmandel,npoints), y(nmandel,npoints)
        real(wp) :: tangent(nmandel,nmandel,npoints)
        real(wp) :: changed(nmandel,npoints), yp(nmandel,npoints), ym(nmandel,npoints)
        real(wp) :: one(nmandel,1)
        real(wp) :: h, numerical, reference, tolerance
        integer :: status, i, j

        call mat%init(orthotropic_elastic_t(100.0_wp,100.0_wp,100.0_wp, &
                       0.25_wp,0.25_wp,0.25_wp,40.0_wp,40.0_wp,40.0_wp), &
                       npoints,status)
        call expect_status(status,cf_ok)
        call single%init(orthotropic_elastic_t(100.0_wp,100.0_wp,100.0_wp, &
                         0.25_wp,0.25_wp,0.25_wp,40.0_wp,40.0_wp,40.0_wp), &
                         1,status)
        call expect_status(status,cf_ok)
        x(:,1) = [0.01_wp,-0.004_wp,0.003_wp,0.006_wp,0.002_wp,-0.005_wp]
        x(:,2) = [-0.008_wp,0.009_wp,0.004_wp,-0.002_wp,0.005_wp,0.003_wp]
        call mat%eval(x,y,tangent,status)
        call expect_status(status,cf_ok)
        call single%eval(x(:,1:1),one,status=status)
        call expect_status(status,cf_ok)
        do i = 1, nmandel
            call close(y(i,1),one(i,1))
        end do

        do i = 1, nmandel
            do j = 1, nmandel
                reference = 0.0_wp
                if (i <= ndim .and. j <= ndim) reference = 40.0_wp
                if (i == j) reference = reference+80.0_wp
                call close(tangent(i,j,1),reference)
            end do
        end do

        ! The strain map is linear, so a finite secant avoids small-step cancellation.
        h = 0.001_wp
        do j = 1, nmandel
            changed = x
            changed(j,1) = changed(j,1)+h
            call mat%eval(changed,yp,status=status)
            call expect_status(status,cf_ok)
            changed(j,1) = changed(j,1)-2.0_wp*h
            call mat%eval(changed,ym,status=status)
            call expect_status(status,cf_ok)
            do i = 1, nmandel
                numerical = (yp(i,1)-ym(i,1))/(2.0_wp*h)
                tolerance = 64.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(yp(i,1)),abs(ym(i,1)))/h
                if (.not. ieee_is_finite(numerical) .or. abs(numerical-tangent(i,j,1)) > tolerance) then
                    error stop "Orthotropic Mandel consistent tangent mismatch"
                end if
            end do
        end do
    end subroutine test_isotropic_limit_and_tangent

    subroutine perturb(model, parameter_id, delta)
        type(orthotropic_elastic_t), intent(inout) :: model
        integer, intent(in) :: parameter_id
        real(wp), intent(in) :: delta

        select case (parameter_id)
        case (1)
            model%young_1 = model%young_1+delta
        case (2)
            model%young_2 = model%young_2+delta
        case (3)
            model%young_3 = model%young_3+delta
        case (4)
            model%poisson_12 = model%poisson_12+delta
        case (5)
            model%poisson_13 = model%poisson_13+delta
        case (6)
            model%poisson_23 = model%poisson_23+delta
        case (7)
            model%shear_23 = model%shear_23+delta
        case (8)
            model%shear_13 = model%shear_13+delta
        case (9)
            model%shear_12 = model%shear_12+delta
        case default
            error stop "Invalid orthotropic parameter ID"
        end select
    end subroutine perturb

    subroutine test_all_sensitivities()
        type(material_t) :: tracked, plus, minus
        type(orthotropic_elastic_t) :: positive_model, negative_model
        real(wp) :: x(nmandel,npoints), y(nmandel,npoints)
        real(wp) :: dsdp(nmandel,nparameters,npoints)
        real(wp) :: yp(nmandel,npoints), ym(nmandel,npoints)
        real(wp) :: coarse(nmandel,npoints)
        real(wp), parameter :: parameter_scale(nparameters) = &
            [100.0_wp,80.0_wp,60.0_wp,1.0_wp,1.0_wp,1.0_wp,25.0_wp,30.0_wp,35.0_wp]
        real(wp) :: numerical, h, tolerance
        integer :: status, id, k, ipoint, i, refinement

        call tracked%init(reference_model(),npoints,status, &
            sensitivity_parameters=[character(len=12) :: &
            "shear_12","shear_13","shear_23","poisson_23","poisson_13", &
            "poisson_12","young_3","young_2","young_1"])
        call expect_status(status,cf_ok)
        if (tracked%sensitivity_size() /= nparameters) error stop "Selected parameter size"
        if (tracked%state_storage_bytes() /= 0) error stop "Unexpected stateless buffers"
        x(:,1) = [0.015_wp,-0.009_wp,0.004_wp,0.002_wp,0.007_wp,-0.003_wp]
        x(:,2) = [-0.011_wp,0.006_wp,0.008_wp,-0.005_wp,0.003_wp,0.004_wp]
        call tracked%eval(x,y,status=status,dsdp=dsdp)
        call expect_status(status,cf_ok)

        do k = 1, nparameters
            id = nparameters+1-k
            h = epsilon(1.0_wp)**(1.0_wp/3.0_wp)*parameter_scale(id)
            do refinement = 1, 2
                positive_model = reference_model()
                negative_model = reference_model()
                call perturb(positive_model,id,h)
                call perturb(negative_model,id,-h)
                call plus%init(positive_model,npoints,status)
                call expect_status(status,cf_ok)
                call minus%init(negative_model,npoints,status)
                call expect_status(status,cf_ok)
                call plus%eval(x,yp,status=status)
                call expect_status(status,cf_ok)
                call minus%eval(x,ym,status=status)
                call expect_status(status,cf_ok)
                if (refinement == 1) then
                    coarse = (yp-ym)/(2.0_wp*h)
                    h = 0.5_wp*h
                end if
            end do
            do ipoint = 1, npoints
                do i = 1, nmandel
                    numerical = (yp(i,ipoint)-ym(i,ipoint))/(2.0_wp*h)
                    ! The step pair estimates O(h**2) truncation, with a roundoff floor.
                    tolerance = abs(numerical-coarse(i,ipoint))/3.0_wp + &
                        64.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(yp(i,ipoint)),abs(ym(i,ipoint)))/h
                    if (.not. ieee_is_finite(numerical) .or. .not. ieee_is_finite(coarse(i,ipoint)) .or. &
                        abs(numerical-dsdp(i,k,ipoint)) > tolerance) then
                        error stop "Orthotropic nine-parameter sensitivity mismatch"
                    end if
                end do
            end do
        end do
        call tracked%commit(status)
        call expect_status(status,cf_ok)
        call tracked%eval(x,y,status=status)
        call expect_status(status,cf_invalid_argument)
        call tracked%eval(x,y,status=status,dsdp=dsdp)
        call expect_status(status,cf_ok)
    end subroutine test_all_sensitivities

    subroutine test_modulus_scales()
        type(material_t) :: mat
        real(wp) :: x(nmandel,1), y(nmandel,1), dsdp(nmandel,1,1)
        real(wp) :: modulus
        integer :: status, case_id, i

        ! The compliance determinant may underflow or overflow while C remains finite.
        x = 0.0_wp
        x(1,1) = 1.0_wp/1024.0_wp
        do case_id = 1, 2
            if (case_id == 1) then
                modulus = huge(1.0_wp)**0.4_wp
            else
                modulus = tiny(1.0_wp)**0.4_wp
            end if
            call mat%init(orthotropic_elastic_t(modulus,modulus,modulus, &
                0.0_wp,0.0_wp,0.0_wp,1.0_wp,1.0_wp,1.0_wp),1,status)
            call expect_status(status,cf_ok)
            call mat%eval(x,y,status=status)
            call expect_status(status,cf_ok)
            call close(y(1,1)/modulus,x(1,1))
            do i = 2, nmandel
                call close(y(i,1),0.0_wp)
            end do
        end do

        ! d sigma_11 / d E_1 = epsilon_11 even when 1/E_1**2 is unrepresentable.
        modulus = huge(1.0_wp)**0.75_wp
        call mat%init(orthotropic_elastic_t(modulus,1.0_wp,1.0_wp, &
            0.0_wp,0.0_wp,0.0_wp,1.0_wp,1.0_wp,1.0_wp),1,status, &
            sensitivity_parameters=[character(len=7) :: "young_1"])
        call expect_status(status,cf_ok)
        call mat%eval(x,y,status=status,dsdp=dsdp)
        call expect_status(status,cf_ok)
        call close(y(1,1)/modulus,x(1,1))
        call close(dsdp(1,1,1),x(1,1))
        do i = 2, nmandel
            call close(dsdp(i,1,1),0.0_wp)
        end do
    end subroutine test_modulus_scales

    subroutine test_invalid_and_trial()
        type(material_t) :: mat
        type(orthotropic_elastic_t) :: model
        real(wp) :: x(nmandel,npoints), y(nmandel,npoints)
        real(wp) :: bad(nmandel-1,npoints)
        integer :: status

        ! Anisotropy permits a Poisson ratio above 0.5 and an auxetic coupling.
        model = orthotropic_elastic_t(100.0_wp,25.0_wp,50.0_wp, &
            0.75_wp,-0.1_wp,0.05_wp,10.0_wp,15.0_wp,20.0_wp)
        call mat%init(model,npoints,status)
        call expect_status(status,cf_ok)
        ! Positive first and second leading minors do not suffice for stability.
        model = orthotropic_elastic_t(100.0_wp,100.0_wp,100.0_wp, &
            0.6_wp,0.6_wp,0.6_wp,40.0_wp,40.0_wp,40.0_wp)
        call mat%init(model,npoints,status)
        call expect_status(status,cf_invalid_argument)
        model = reference_model()
        model%young_1 = -1.0_wp
        call mat%init(model,npoints,status)
        call expect_status(status,cf_invalid_argument)
        model = reference_model()
        model%shear_23 = 0.0_wp
        call mat%init(model,npoints,status)
        call expect_status(status,cf_invalid_argument)
        model = reference_model()
        model%poisson_12 = 3.0_wp
        call mat%init(model,npoints,status)
        call expect_status(status,cf_invalid_argument)
        model = reference_model()
        model%poisson_13 = ieee_value(0.0_wp,ieee_quiet_nan)
        call mat%init(model,npoints,status)
        call expect_status(status,cf_invalid_argument)
        call mat%init(reference_model(),npoints,status)
        call expect_status(status,cf_ok)

        x = 0.0_wp
        x(1,1) = 0.01_wp
        x(2,2) = ieee_value(0.0_wp,ieee_quiet_nan)
        call mat%eval(x,y,status=status)
        call expect_status(status,cf_evaluation_failed)
        if (mat%last_failure_point() /= 2) error stop "Wrong failing orthotropic point"
        call mat%commit(status)
        call expect_status(status,cf_no_trial)
        x(2,2) = 0.02_wp
        call mat%eval(x,y,status=status)
        call expect_status(status,cf_ok)
        call mat%revert()
        call mat%commit(status)
        call expect_status(status,cf_no_trial)
        bad = 0.0_wp
        call mat%eval(bad,y,status=status)
        call expect_status(status,cf_invalid_shape)
        model = reference_model()
        model%shear_13 = -1.0_wp
        call mat%init(model,npoints,status)
        call expect_status(status,cf_invalid_argument)
        call mat%eval(x,y,status=status)
        call expect_status(status,cf_ok)
        call mat%commit(status)
        call expect_status(status,cf_ok)
    end subroutine test_invalid_and_trial

end program test_orthotropic_elastic

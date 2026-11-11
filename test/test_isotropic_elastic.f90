! SPDX-License-Identifier: Apache-2.0

program test_isotropic_elastic
    !! Test convention: stress and stress-like material parameters in Pa.
    !! Strain and deformation gradients are dimensionless. Explicit rescaling
    !! cases state their factors.
    !!
    !! Independent numerical and API checks for the first library material.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
    use constituf, only: wp, material_t, cf_ok, cf_invalid_argument, &
                        cf_invalid_shape, cf_evaluation_failed, cf_no_trial
    use constituf_isotropic_elastic, only: isotropic_elastic_t
    implicit none (type, external)

    integer, parameter :: nmandel = 6

    call test_metadata()
    call test_physical_limits()
    call test_reference_response()
    call test_unit_conventions()
    call test_rotation_and_batch()
    call test_sensitivities()
    call test_invalid_input()

contains

    subroutine expect_status(actual, expected, context)
        integer, intent(in) :: actual, expected
        character(*), intent(in) :: context

        if (actual /= expected) then
            print *, trim(context), " expected status ", expected, " got ", actual
            error stop "Isotropic material status mismatch"
        end if
    end subroutine expect_status

    subroutine expect_close(actual, expected, context)
        real(wp), intent(in) :: actual, expected
        character(*), intent(in) :: context
        real(wp) :: tolerance

        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) then
            print *, trim(context), " expected ", expected, " got ", actual
            error stop "Nonfinite isotropic material comparison"
        end if
        tolerance = 256.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(expected))
        if (abs(actual-expected) > tolerance) then
            print *, trim(context), " expected ", expected, " got ", actual
            error stop "Isotropic material numerical mismatch"
        end if
    end subroutine expect_close

    subroutine test_metadata()
        type(material_t) :: mat
        character(:), allocatable :: convention, measure, frame
        integer, allocatable :: shape(:)
        integer :: status
        logical :: found

        call mat%init(isotropic_elastic_t(100.0_wp,40.0_wp), 2, status)
        call expect_status(status, cf_ok, "metadata init")
        if (mat%input_size() /= 6 .or. mat%output_size() /= 6) then
            error stop "Mandel packed dimensions"
        end if
        if (mat%state_size() /= 0 .or. mat%integer_state_size() /= 0 .or. &
            mat%logical_state_size() /= 0) error stop "Material must be stateless"
        if (mat%state_storage_bytes() /= 0) error stop "Stateless retained memory"
        if (mat%input_offset("strain.small") /= 1) error stop "Strain offset"
        if (mat%output_offset("stress.cauchy") /= 1) error stop "Stress offset"
        if (.not. mat%field_matches("input", "strain.small", "mandel", &
                                   "small_strain", "spatial", [3,3])) then
            error stop "Strain schema"
        end if
        if (.not. mat%field_matches("output", "stress.cauchy", "mandel", &
                                   "cauchy_stress", "spatial", [3,3])) then
            error stop "Stress schema"
        end if
        if (mat%field_matches("output", "stress.cauchy", "mandel", &
                              "piola_stress", "spatial", [3,3])) then
            error stop "Incorrect stress measure accepted"
        end if
        if (mat%field_matches("input", "strain.small", "mandel", &
                              "small_strain", "reference", [3,3])) then
            error stop "Incorrect frame accepted"
        end if
        call mat%field_info("output", "stress.cauchy", convention, &
                            measure, frame, shape, found)
        if (.not. found) error stop "Stress field introspection"
        if (convention /= "mandel" .or. &
            measure /= "cauchy_stress" .or. frame /= "spatial") then
            error stop "Unexpected field metadata"
        end if
        if (any(shape /= [3,3])) error stop "Incorrect physical tensor shape"
    end subroutine test_metadata

    subroutine test_physical_limits()
        !! Independent hydrostatic, physical shear and auxetic uniaxial-stress cases.
        type(material_t) :: mat
        real(wp), parameter :: young = 120.0_wp, poisson = -0.25_wp
        real(wp) :: x(nmandel,2), y(nmandel,2), reference(nmandel,2)
        integer :: status, i, j

        call mat%init(isotropic_elastic_t(100.0_wp,40.0_wp), 2, status)
        call expect_status(status, cf_ok, "physical limits init")
        x = 0.0_wp
        x(1:3,1) = 0.001_wp
        ! Tensor shear epsilon_12 = 0.001 gives sigma_12 = 0.08 Pa.
        x(6,2) = sqrt(2.0_wp)*0.001_wp
        reference(:,1) = [0.3_wp,0.3_wp,0.3_wp,0.0_wp,0.0_wp,0.0_wp]
        reference(:,2) = [0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,sqrt(2.0_wp)*0.08_wp]
        call mat%eval(x, y, status=status)
        call expect_status(status, cf_ok, "hydrostatic and physical shear")
        do j = 1, 2
            do i = 1, nmandel
                call expect_close(y(i,j), reference(i,j), "independent physical response")
            end do
        end do

        ! E and nu prescribe the free lateral expansion independently of the K/G split.
        call mat%init(isotropic_elastic_t(young/(3.0_wp*(1.0_wp-2.0_wp*poisson)), &
                                         young/(2.0_wp*(1.0_wp+poisson))), 1, status)
        call expect_status(status, cf_ok, "auxetic init")
        x(:,1) = [0.001_wp,-poisson*0.001_wp,-poisson*0.001_wp,0.0_wp,0.0_wp,0.0_wp]
        reference(:,1) = [young*0.001_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp]
        call mat%eval(x(:,1:1), y(:,1:1), status=status)
        call expect_status(status, cf_ok, "auxetic uniaxial stress")
        do i = 1, nmandel
            call expect_close(y(i,1), reference(i,1), "auxetic free lateral expansion")
        end do
    end subroutine test_physical_limits

    subroutine test_reference_response()
        !! Independently evaluate the closed-form isotropic tangent and stress.
        type(material_t) :: mat
        real(wp) :: x(nmandel,1), y(nmandel,1), tangent(nmandel,nmandel,1)
        real(wp) :: yp(nmandel,1), ym(nmandel,1), trial(nmandel,1)
        real(wp) :: reference(nmandel)
        real(wp) :: lambda, h, derivative, trace_strain
        integer :: status, i, j

        call mat%init(isotropic_elastic_t(100.0_wp,40.0_wp), 1, status)
        call expect_status(status, cf_ok, "reference init")
        x(:,1) = [0.01_wp,0.02_wp,-0.01_wp,0.005_wp,-0.004_wp,0.003_wp]
        trace_strain = sum(x(1:3,1))
        reference(1:3) = 100.0_wp*trace_strain + &
                         80.0_wp*(x(1:3,1)-trace_strain/3.0_wp)
        reference(4:6) = 80.0_wp*x(4:6,1)
        call mat%eval(x, y, tangent, status)
        call expect_status(status, cf_ok, "reference evaluation")
        if (.not. all(ieee_is_finite(tangent))) error stop "Nonfinite Mandel tangent"
        do i = 1, 6
            call expect_close(y(i,1), reference(i), "closed-form stress")
        end do

        lambda = 100.0_wp - 80.0_wp/3.0_wp
        do j = 1, 6
            do i = 1, 6
                if (i <= 3 .and. j <= 3) then
                    if (i == j) then
                        call expect_close(tangent(i,j,1), lambda+80.0_wp, "normal diagonal")
                    else
                        call expect_close(tangent(i,j,1), lambda, "normal coupling")
                    end if
                else if (i == j) then
                    call expect_close(tangent(i,j,1), 80.0_wp, "Mandel shear")
                else
                    call expect_close(tangent(i,j,1), 0.0_wp, "zero coupling")
                end if
            end do
        end do
        if (0.5_wp*dot_product(x(:,1),y(:,1)) <= 0.0_wp) then
            error stop "Positive strain must store positive energy"
        end if

        h = max(1.0e-6_wp, sqrt(epsilon(1.0_wp)))
        do j = 1, 6
            trial = x
            trial(j,1) = trial(j,1) + h
            call mat%eval(trial, yp, status=status)
            call expect_status(status, cf_ok, "positive strain perturbation")
            trial(j,1) = trial(j,1) - 2.0_wp*h
            call mat%eval(trial, ym, status=status)
            call expect_status(status, cf_ok, "negative strain perturbation")
            do i = 1, 6
                derivative = (yp(i,1)-ym(i,1))/(2.0_wp*h)
                if (.not. ieee_is_finite(derivative)) error stop "Nonfinite strain difference"
                if (abs(derivative-tangent(i,j,1)) > &
                    max(1.0e-7_wp*max(1.0_wp,abs(tangent(i,j,1))), &
                        16.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(yp(i,1)),abs(ym(i,1)))/h)) then
                    error stop "Inconsistent Mandel algorithmic tangent"
                end if
            end do
        end do
        call mat%revert()
        call mat%commit(status)
        call expect_status(status, cf_no_trial, "rejected trial cannot commit")
        call mat%eval(x, y, status=status)
        call expect_status(status, cf_ok, "repeated stateless evaluation")
        call mat%commit(status)
        call expect_status(status, cf_ok, "stateless commit")
        call mat%eval(x, y, status=status)
        call expect_status(status, cf_ok, "stateless repeated loading")
        do i = 1, 6
            call expect_close(y(i,1), reference(i), "stateless reload")
        end do
    end subroutine test_reference_response

    subroutine test_unit_conventions()
        !! The same E=150 GPa, nu=1/4 specimen is supplied in Pa and MPa.
        !! Strain is dimensionless. Stress and tangent scale together, while
        !! derivatives of stress with respect to the moduli are dimensionless.
        real(wp), parameter :: stress_factors(2) = [1.0_wp,1.0e-6_wp]
        real(wp), parameter :: expected_stress(nmandel) = &
            [180.0e6_wp,60.0e6_wp,60.0e6_wp,0.0_wp,0.0_wp,0.0_wp]
        real(wp), parameter :: expected_shear_derivative(nmandel) = &
            [4.0_wp/3000.0_wp,-2.0_wp/3000.0_wp,-2.0_wp/3000.0_wp,0.0_wp,0.0_wp,0.0_wp]
        type(material_t) :: material
        real(wp) :: x(nmandel,1), stress(nmandel,1), tangent(nmandel,nmandel,1)
        real(wp) :: dsdp(nmandel,2,1), expected_tangent(nmandel,nmandel)
        real(wp) :: factor
        integer :: status, system, row, column

        x = 0.0_wp
        x(1,1) = 0.001_wp
        expected_tangent = 0.0_wp
        expected_tangent(1:3,1:3) = 60.0e9_wp
        do row = 1,3
            expected_tangent(row,row) = 180.0e9_wp
        end do
        do row = 4,nmandel
            expected_tangent(row,row) = 120.0e9_wp
        end do
        do system = 1,size(stress_factors)
            factor = stress_factors(system)
            call material%init(isotropic_elastic_t(100.0e9_wp*factor,60.0e9_wp*factor),1,status, &
                sensitivity_parameters=[character(len=16) :: "bulk_modulus","shear_modulus"])
            call expect_status(status,cf_ok,"unit convention setup")
            call material%eval(x,stress,tangent,status,dsdp=dsdp)
            call expect_status(status,cf_ok,"unit convention evaluation")
            do row = 1,nmandel
                call expect_close(stress(row,1)/factor,expected_stress(row),"Pa/MPa stress")
                call expect_close(dsdp(row,1,1),merge(0.001_wp,0.0_wp,row <= 3),"Pa/MPa bulk derivative")
                call expect_close(dsdp(row,2,1),expected_shear_derivative(row),"Pa/MPa shear derivative")
                do column = 1,nmandel
                    call expect_close(tangent(row,column,1)/factor,expected_tangent(row,column),"Pa/MPa tangent")
                end do
            end do
        end do
    end subroutine test_unit_conventions

    subroutine test_rotation_and_batch()
        !! A 45-degree rotation about axis 3 mixes normal and Mandel shear components.
        type(material_t) :: batched, single
        real(wp) :: x(nmandel,2), y(nmandel,2), one(nmandel,1)
        real(wp) :: root_two
        integer :: status, i

        call batched%init(isotropic_elastic_t(100.0_wp,40.0_wp), 2, status)
        call expect_status(status, cf_ok, "batch init")
        call single%init(isotropic_elastic_t(100.0_wp,40.0_wp), 1, status)
        call expect_status(status, cf_ok, "single init")

        x(:,1) = [0.02_wp,-0.01_wp,0.005_wp,0.006_wp,-0.004_wp,0.003_wp]
        root_two = sqrt(2.0_wp)
        x(:,2) = [0.5_wp*(x(1,1)+x(2,1))-x(6,1)/root_two, &
                  0.5_wp*(x(1,1)+x(2,1))+x(6,1)/root_two, x(3,1), &
                  (x(4,1)+x(5,1))/root_two, (x(5,1)-x(4,1))/root_two, &
                  (x(1,1)-x(2,1))/root_two]
        call batched%eval(x, y, status=status)
        call expect_status(status, cf_ok, "rotated batch")
        call single%eval(x(:,1:1), one, status=status)
        call expect_status(status, cf_ok, "single point")
        do i = 1, 6
            call expect_close(y(i,1), one(i,1), "single versus batch")
        end do
        call expect_close(y(1,2), 0.5_wp*(y(1,1)+y(2,1))-y(6,1)/root_two, "rotated 11")
        call expect_close(y(2,2), 0.5_wp*(y(1,1)+y(2,1))+y(6,1)/root_two, "rotated 22")
        call expect_close(y(3,2), y(3,1), "rotated 33")
        call expect_close(y(4,2), (y(4,1)+y(5,1))/root_two, "rotated 23")
        call expect_close(y(5,2), (y(5,1)-y(4,1))/root_two, "rotated 13")
        call expect_close(y(6,2), (y(1,1)-y(2,1))/root_two, "rotated 12")
        call expect_close(dot_product(x(:,1),y(:,1)), &
                          dot_product(x(:,2),y(:,2)), "rotational energy")
    end subroutine test_rotation_and_batch

    subroutine test_sensitivities()
        type(material_t) :: mat, plus, minus
        real(wp) :: x(nmandel,2), y(nmandel,2)
        real(wp) :: tangent(nmandel,nmandel,2), dsdp(nmandel,2,2)
        real(wp) :: wrong(nmandel-1,2,2)
        real(wp) :: yp(nmandel,2), ym(nmandel,2)
        real(wp) :: trace_strain, h, expected
        integer :: status, i, j, k

        call mat%init(isotropic_elastic_t(100.0_wp,40.0_wp), 2, status, &
            sensitivity_parameters=[character(len=16) :: "shear_modulus", "bulk_modulus"])
        call expect_status(status, cf_ok, "sensitivity init")
        if (mat%sensitivity_size() /= 2) error stop "Sensitivity selection count"
        if (mat%state_storage_bytes() /= 0) error stop "Stateless derivative storage"

        x(:,1) = [0.02_wp,0.01_wp,-0.01_wp,0.005_wp,-0.006_wp,0.004_wp]
        x(:,2) = [-0.01_wp,0.03_wp,0.02_wp,-0.005_wp,0.001_wp,0.008_wp]
        call mat%eval(x, y, tangent, status, dsdp=dsdp)
        call expect_status(status, cf_ok, "analytical derivatives")
        if (.not. all(ieee_is_finite(tangent))) error stop "Nonfinite sensitivity-path tangent"
        if (.not. all(ieee_is_finite(dsdp))) error stop "Nonfinite parameter sensitivity"
        do i = 1, 2
            trace_strain = sum(x(1:3,i))
            do k = 1, 3
                call expect_close(dsdp(k,1,i), &
                    2.0_wp*(x(k,i)-trace_strain/3.0_wp), "d sigma / d G")
                call expect_close(dsdp(k,2,i), trace_strain, "d sigma / d K")
            end do
            do k = 4, 6
                call expect_close(dsdp(k,1,i), 2.0_wp*x(k,i), "shear d sigma / d G")
                call expect_close(dsdp(k,2,i), 0.0_wp, "shear d sigma / d K")
            end do
        end do

        ! Central differences are verification only; the model's dsdp is exact.
        h = max(1.0e-4_wp, 100.0_wp*sqrt(epsilon(1.0_wp)))
        do j = 1, 2
            if (j == 1) then
                call plus%init(isotropic_elastic_t(100.0_wp,40.0_wp+h), 2, status)
                call expect_status(status, cf_ok, "shear plus")
                call minus%init(isotropic_elastic_t(100.0_wp,40.0_wp-h), 2, status)
            else
                call plus%init(isotropic_elastic_t(100.0_wp+h,40.0_wp), 2, status)
                call expect_status(status, cf_ok, "bulk plus")
                call minus%init(isotropic_elastic_t(100.0_wp-h,40.0_wp), 2, status)
            end if
            call expect_status(status, cf_ok, "perturbed parameter init")
            call plus%eval(x, yp, status=status)
            call expect_status(status, cf_ok, "positive parameter perturbation")
            call minus%eval(x, ym, status=status)
            call expect_status(status, cf_ok, "negative parameter perturbation")
            do i = 1, 2
                do k = 1, 6
                    expected = (yp(k,i)-ym(k,i))/(2.0_wp*h)
                    if (.not. ieee_is_finite(expected)) error stop "Nonfinite parameter difference"
                    if (abs(expected-dsdp(k,j,i)) > &
                        max(1.0e-7_wp, 16.0_wp*epsilon(1.0_wp)* &
                        max(1.0_wp,abs(yp(k,i)),abs(ym(k,i)))/h)) then
                        error stop "Analytical parameter sensitivity mismatch"
                    end if
                end do
            end do
        end do

        call mat%commit(status)
        call expect_status(status, cf_ok, "sensitized commit")
        call mat%eval(x, y, status=status)
        call expect_status(status, cf_invalid_argument, "missing dsdp")
        call mat%eval(x, y, status=status, dsdp=wrong)
        call expect_status(status, cf_invalid_shape, "wrong dsdp shape")
        call mat%eval(x, y, status=status, dsdp=dsdp)
        call expect_status(status, cf_ok, "sensitized retry")
        call mat%revert()

        call mat%init(isotropic_elastic_t(100.0_wp,40.0_wp), 2, status, &
            sensitivity_parameters=[character(len=16) :: "bulk_modulus", "bulk_modulus"])
        call expect_status(status, cf_invalid_argument, "duplicate sensitivity")
        call mat%init(isotropic_elastic_t(100.0_wp,40.0_wp), 2, status, &
            sensitivity_parameters=[character(len=16) :: "missing"])
        call expect_status(status, cf_invalid_argument, "unknown sensitivity")
        if (mat%sensitivity_size() /= 2) error stop "Failed init mutated existing model"
    end subroutine test_sensitivities

    subroutine test_invalid_input()
        type(material_t) :: mat
        real(wp) :: x(nmandel,2), y(nmandel,2), invalid_x(nmandel-1,2)
        real(wp) :: nan, infinity
        integer :: status

        nan = ieee_value(0.0_wp, ieee_quiet_nan)
        infinity = ieee_value(0.0_wp, ieee_positive_inf)
        call mat%init(isotropic_elastic_t(0.0_wp,40.0_wp), 2, status)
        call expect_status(status, cf_invalid_argument, "zero bulk modulus")
        call mat%init(isotropic_elastic_t(100.0_wp,-1.0_wp), 2, status)
        call expect_status(status, cf_invalid_argument, "negative shear modulus")
        call mat%init(isotropic_elastic_t(nan,40.0_wp), 2, status)
        call expect_status(status, cf_invalid_argument, "nonfinite bulk modulus")
        call mat%init(isotropic_elastic_t(100.0_wp,infinity), 2, status)
        call expect_status(status, cf_invalid_argument, "nonfinite shear modulus")
        call mat%init(isotropic_elastic_t(100.0_wp,40.0_wp), 2, status)
        call expect_status(status, cf_ok, "valid model following invalid")

        x = 0.0_wp
        x(1,:) = 0.01_wp
        x(1,2) = nan
        call mat%eval(x, y, status=status)
        call expect_status(status, cf_evaluation_failed, "invalid integration point")
        if (mat%last_failure_point() /= 2) error stop "Wrong failing integration point"
        call mat%commit(status)
        call expect_status(status, cf_no_trial, "invalid trial was committed")
        x(1,2) = 0.01_wp
        call mat%eval(x, y, status=status)
        call expect_status(status, cf_ok, "valid retry")
        call expect_close(y(1,2), 100.0_wp*0.01_wp + &
            80.0_wp*(0.01_wp-0.01_wp/3.0_wp), "recovered response")

        invalid_x = 0.0_wp
        call mat%eval(invalid_x, y, status=status)
        call expect_status(status, cf_invalid_shape, "invalid input shape")
        call mat%init(isotropic_elastic_t(-1.0_wp,40.0_wp), 2, status)
        call expect_status(status, cf_invalid_argument, "failed reinit")
        call mat%eval(x, y, status=status)
        call expect_status(status, cf_ok, "failed reinit preserved model")
    end subroutine test_invalid_input

end program test_isotropic_elastic

! SPDX-License-Identifier: Apache-2.0

program test_isotropic_thermoelastic
    !! Test convention: moduli and stress in Pa, absolute temperatures in K
    !! and expansion coefficient in 1/K. Strain is dimensionless.
    !!
    !! Thermoelastic analytical references and coupled rectangular Jacobian.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use constituf, only: wp, material_t, cf_ok, cf_invalid_argument, &
                        cf_invalid_shape, cf_evaluation_failed, cf_no_trial
    use constituf_isotropic_thermoelastic, only: isotropic_thermoelastic_t
    implicit none (type, external)

    integer, parameter :: ndim = 3, nmandel = 6, ninput = nmandel+1
    integer, parameter :: npoints = 2, nparameters = 4
    real(wp), parameter :: bulk = 100.0_wp, shear = 40.0_wp
    real(wp), parameter :: expansion = 1.0e-4_wp, reference_t = 300.0_wp

    call test_schema_and_thermal_references()
    call test_coupled_tangent()
    call test_parameter_sensitivities()
    call test_unit_conventions()
    call test_large_elastic_limit()
    call test_rejected_trials()

contains

    subroutine expect_status(value, expected)
        integer, intent(in) :: value, expected

        if (value /= expected) then
            print *, "Expected status", expected, "got", value
            error stop "Thermoelastic API status mismatch"
        end if
    end subroutine expect_status

    subroutine expect_close(actual, expected)
        real(wp), intent(in) :: actual, expected
        real(wp) :: tolerance

        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) then
            error stop "Nonfinite thermoelastic comparison"
        end if
        tolerance = 512.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(expected))
        if (abs(actual-expected) > tolerance) then
            print *, "Expected", expected, "got", actual
            error stop "Thermoelastic reference mismatch"
        end if
    end subroutine expect_close

    subroutine test_unit_conventions()
        !! One physical case in Pa,K and MPa,degR. Rankine is absolute:
        !! T_R=(9/5)*T_K and alpha_R=(5/9)*alpha_K, with no offset.
        !! Each tangent and sensitivity column follows its own input or parameter.
        real(wp), parameter :: stress_factor = 1.0e-6_wp, temperature_factor = 9.0_wp/5.0_wp
        real(wp), parameter :: parameter_factors(nparameters) = &
            [stress_factor,stress_factor,1.0_wp/temperature_factor,temperature_factor]
        type(material_t) :: reference, scaled
        real(wp) :: x(ninput,1), scaled_input(ninput,1), input_factors(ninput)
        real(wp) :: stress(nmandel,1), scaled_stress(nmandel,1)
        real(wp) :: tangent(nmandel,ninput,1), scaled_tangent(nmandel,ninput,1)
        real(wp) :: dsdp(nmandel,nparameters,1), scaled_dsdp(nmandel,nparameters,1)
        integer :: status, row, column

        call reference%init(isotropic_thermoelastic_t(100.0e9_wp,40.0e9_wp,1.0e-5_wp,300.0_wp), &
            1,status,sensitivity_parameters=[character(len=24) :: &
                "bulk_modulus","shear_modulus","thermal_expansion","reference_temperature"])
        call expect_status(status,cf_ok)
        call scaled%init(isotropic_thermoelastic_t(100.0e9_wp*stress_factor,40.0e9_wp*stress_factor, &
            1.0e-5_wp/temperature_factor,300.0_wp*temperature_factor), &
            1,status,sensitivity_parameters=[character(len=24) :: &
                "bulk_modulus","shear_modulus","thermal_expansion","reference_temperature"])
        call expect_status(status,cf_ok)
        x(:,1) = [0.002_wp,-0.001_wp,0.0005_wp,0.0004_wp,-0.0002_wp,0.0003_wp,315.0_wp]
        input_factors = 1.0_wp
        input_factors(ninput) = temperature_factor
        scaled_input(:,1) = input_factors*x(:,1)
        call reference%eval(x,stress,tangent,status,dsdp=dsdp)
        call expect_status(status,cf_ok)
        call scaled%eval(scaled_input,scaled_stress,scaled_tangent,status,dsdp=scaled_dsdp)
        call expect_status(status,cf_ok)
        do row = 1,nmandel
            call expect_close(scaled_stress(row,1)/stress_factor,stress(row,1))
            do column = 1,ninput
                call expect_close(scaled_tangent(row,column,1)*input_factors(column)/stress_factor, &
                    tangent(row,column,1))
            end do
            do column = 1,nparameters
                call expect_close(scaled_dsdp(row,column,1)*parameter_factors(column)/stress_factor, &
                    dsdp(row,column,1))
            end do
        end do
    end subroutine test_unit_conventions

    subroutine test_schema_and_thermal_references()
        type(material_t) :: material
        real(wp) :: x(ninput,npoints), stress(nmandel,npoints)
        real(wp) :: tangent(nmandel,ninput,npoints)
        real(wp) :: temperature_change
        integer :: status, i, j
        integer :: empty_shape(0)

        call material%init(isotropic_thermoelastic_t(bulk,shear,expansion,reference_t), &
                           npoints,status)
        call expect_status(status,cf_ok)
        if (material%input_size() /= ninput .or. material%output_size() /= nmandel) then
            error stop "Coupled input/output sizes"
        end if
        if (material%input_offset("strain.small") /= 1) error stop "Strain offset"
        if (material%input_offset("temperature.absolute") /= ninput) then
            error stop "Temperature offset"
        end if
        if (material%state_size() /= 0 .or. material%state_storage_bytes() /= 0) then
            error stop "Thermoelasticity must have zero retained history"
        end if
        if (.not. material%field_matches("input","strain.small","mandel", &
                                         "small_strain","spatial",[ndim,ndim])) then
            error stop "Mandel strain field mismatch"
        end if
        if (.not. material%field_matches("input","temperature.absolute","scalar", &
                                         "absolute_temperature","material",empty_shape)) then
            error stop "Temperature metadata mismatch"
        end if
        if (material%field_matches("input","temperature.absolute","scalar", &
                                    "temperature_difference","material",empty_shape)) then
            error stop "Temperature difference accepted as absolute temperature"
        end if
        if (.not. material%field_matches("output","stress.cauchy","mandel", &
                                         "cauchy_stress","spatial",[ndim,ndim])) then
            error stop "Cauchy stress metadata mismatch"
        end if

        temperature_change = 15.0_wp
        x = 0.0_wp
        x(ninput,:) = reference_t+temperature_change
        x(1:ndim,1) = expansion*temperature_change
        call material%eval(x,stress,tangent,status)
        call expect_status(status,cf_ok)
        if (maxval(abs(stress(:,1))) > &
            64.0_wp*epsilon(1.0_wp)*(3.0_wp*bulk+2.0_wp*shear)*abs(expansion*temperature_change)) then
            error stop "Free thermal expansion produced stress"
        end if
        do i = 1, ndim
            call expect_close(stress(i,2),-3.0_wp*bulk*expansion*temperature_change)
            call expect_close(tangent(i,ninput,1),-3.0_wp*bulk*expansion)
            do j = 1, ndim
                if (i == j) then
                    call expect_close(tangent(i,j,1),bulk+4.0_wp*shear/3.0_wp)
                else
                    call expect_close(tangent(i,j,1),bulk-2.0_wp*shear/3.0_wp)
                end if
            end do
        end do
        call expect_close(tangent(4,4,1),2.0_wp*shear)
        call expect_close(tangent(5,ninput,1),0.0_wp)
        call material%commit(status)
        call expect_status(status,cf_ok)
        call material%eval(x,stress,status=status)
        call expect_status(status,cf_ok)
        call material%revert()
        call material%commit(status)
        call expect_status(status,cf_no_trial)
    end subroutine test_schema_and_thermal_references

    subroutine test_coupled_tangent()
        type(material_t) :: material, single
        real(wp) :: x(ninput,npoints), stress(nmandel,npoints)
        real(wp) :: tangent(nmandel,ninput,npoints)
        real(wp) :: perturbed(ninput,npoints), plus(nmandel,npoints)
        real(wp) :: minus(nmandel,npoints), alone(nmandel,1)
        real(wp) :: h, numerical, reference, tolerance
        integer :: status, j, i

        call material%init(isotropic_thermoelastic_t(bulk,shear,expansion,reference_t), &
                           npoints,status)
        call expect_status(status,cf_ok)
        call single%init(isotropic_thermoelastic_t(bulk,shear,expansion,reference_t), &
                         1,status)
        call expect_status(status,cf_ok)
        x(:,1) = [0.003_wp,-0.001_wp,0.004_wp,0.005_wp,-0.002_wp,0.007_wp,310.0_wp]
        x(:,2) = [-0.005_wp,0.008_wp,0.002_wp,-0.003_wp,0.006_wp,0.002_wp,290.0_wp]
        call material%eval(x,stress,tangent,status)
        call expect_status(status,cf_ok)
        call single%eval(x(:,1:1),alone,status=status)
        call expect_status(status,cf_ok)
        do i = 1, nmandel
            call expect_close(stress(i,1),alone(i,1))
        end do

        ! Independent direct isotropic stress including its thermal eigenstrain.
        reference = bulk*(sum(x(1:ndim,1))-3.0_wp*expansion*10.0_wp) + &
                    2.0_wp*shear*(x(1,1)-sum(x(1:ndim,1))/3.0_wp)
        call expect_close(stress(1,1),reference)
        call expect_close(stress(4,1),2.0_wp*shear*x(4,1))
        do j = 1, ninput
            ! Each input enters affinely. Use a resolved secant in its physical units.
            h = 0.001_wp
            if (j == ninput) h = 1.0_wp
            perturbed = x
            perturbed(j,1) = perturbed(j,1)+h
            call material%eval(perturbed,plus,status=status)
            call expect_status(status,cf_ok)
            perturbed(j,1) = perturbed(j,1)-2.0_wp*h
            call material%eval(perturbed,minus,status=status)
            call expect_status(status,cf_ok)
            do i = 1, nmandel
                numerical = (plus(i,1)-minus(i,1))/(2.0_wp*h)
                tolerance = 64.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(plus(i,1)),abs(minus(i,1)))/h
                if (.not. ieee_is_finite(numerical) .or. .not. ieee_is_finite(tangent(i,j,1)) .or. &
                    abs(numerical-tangent(i,j,1)) > tolerance) then
                    error stop "Incorrect coupled thermoelastic Jacobian"
                end if
            end do
        end do
    end subroutine test_coupled_tangent

    subroutine perturb_parameter(model,id,h)
        type(isotropic_thermoelastic_t), intent(inout) :: model
        integer, intent(in) :: id
        real(wp), intent(in) :: h

        select case (id)
        case (1)
            model%bulk_modulus = model%bulk_modulus+h
        case (2)
            model%shear_modulus = model%shear_modulus+h
        case (3)
            model%thermal_expansion = model%thermal_expansion+h
        case (4)
            model%reference_temperature = model%reference_temperature+h
        case default
            error stop "Unknown thermoelastic parameter"
        end select
    end subroutine perturb_parameter

    subroutine test_parameter_sensitivities()
        type(material_t) :: tracked, plus_mat, minus_mat
        type(isotropic_thermoelastic_t) :: plus_model, minus_model
        real(wp) :: x(ninput,npoints), stress(nmandel,npoints)
        real(wp) :: dsdp(nmandel,nparameters,npoints)
        real(wp) :: tangent(nmandel,ninput,npoints)
        real(wp) :: yp(nmandel,npoints), ym(nmandel,npoints)
        real(wp), parameter :: parameters(nparameters) = [bulk,shear,expansion,reference_t]
        real(wp) :: h, numerical, tolerance
        integer :: status, k, id, i, ipoint

        call tracked%init(isotropic_thermoelastic_t(bulk,shear,expansion,reference_t), &
            npoints,status,sensitivity_parameters=[character(len=24) :: &
            "reference_temperature","thermal_expansion","shear_modulus","bulk_modulus"])
        call expect_status(status,cf_ok)
        if (tracked%sensitivity_size() /= nparameters) error stop "Sensitivity count"
        if (tracked%state_storage_bytes() /= 0) error stop "Derivative history allocated"
        x(:,1) = [0.006_wp,-0.003_wp,0.004_wp,0.002_wp,0.005_wp,-0.003_wp,311.0_wp]
        x(:,2) = [-0.004_wp,0.008_wp,0.002_wp,0.005_wp,-0.001_wp,0.004_wp,282.0_wp]
        call tracked%eval(x,stress,tangent,status,dsdp=dsdp)
        call expect_status(status,cf_ok)
        do k = 1, nparameters
            id = nparameters+1-k
            ! The law is affine in each parameter at fixed remaining parameters.
            h = 0.25_wp*parameters(id)
            plus_model = isotropic_thermoelastic_t(bulk,shear,expansion,reference_t)
            minus_model = plus_model
            call perturb_parameter(plus_model,id,h)
            call perturb_parameter(minus_model,id,-h)
            call plus_mat%init(plus_model,npoints,status)
            call expect_status(status,cf_ok)
            call minus_mat%init(minus_model,npoints,status)
            call expect_status(status,cf_ok)
            call plus_mat%eval(x,yp,status=status)
            call expect_status(status,cf_ok)
            call minus_mat%eval(x,ym,status=status)
            call expect_status(status,cf_ok)
            do ipoint = 1, npoints
                do i = 1, nmandel
                    numerical = (yp(i,ipoint)-ym(i,ipoint))/(2.0_wp*h)
                    tolerance = 64.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(yp(i,ipoint)),abs(ym(i,ipoint)))/h
                    if (.not. ieee_is_finite(numerical) .or. .not. ieee_is_finite(dsdp(i,k,ipoint)) .or. &
                        abs(numerical-dsdp(i,k,ipoint)) > tolerance) then
                        error stop "Thermoelastic parameter derivative mismatch"
                    end if
                end do
            end do
        end do
        call tracked%commit(status)
        call expect_status(status,cf_ok)
        call tracked%eval(x,stress,status=status)
        call expect_status(status,cf_invalid_argument)
        call tracked%eval(x,stress,status=status,dsdp=dsdp)
        call expect_status(status,cf_ok)
    end subroutine test_parameter_sensitivities

    subroutine test_large_elastic_limit()
        type(material_t) :: material
        real(wp) :: x(ninput,npoints), stress(nmandel,npoints)
        real(wp) :: tangent(nmandel,ninput,npoints), dsdp(nmandel,nparameters,npoints)
        integer :: status, i

        ! Alpha = 0 is an elastic limit even when the intermediate 3*K overflows.
        call material%init(isotropic_thermoelastic_t(0.4_wp*huge(1.0_wp), &
            0.1_wp*huge(1.0_wp),0.0_wp,reference_t),npoints,status, &
            sensitivity_parameters=[character(len=24) :: &
            "reference_temperature","thermal_expansion","shear_modulus","bulk_modulus"])
        call expect_status(status,cf_ok)
        x = 0.0_wp
        x(ninput,:) = reference_t
        x(ninput,2) = reference_t+0.125_wp
        call material%eval(x,stress,tangent,status,dsdp=dsdp)
        call expect_status(status,cf_ok)
        if (.not. all(ieee_is_finite(tangent)) .or. .not. all(ieee_is_finite(dsdp))) then
            error stop "Nonfinite thermoelastic elastic limit"
        end if
        do i = 1, nmandel
            call expect_close(stress(i,1),0.0_wp)
            call expect_close(stress(i,2),0.0_wp)
            call expect_close(tangent(i,ninput,1),0.0_wp)
            call expect_close(dsdp(i,2,1),0.0_wp)
            if (i <= ndim) then
                call expect_close(dsdp(i,2,2)/huge(1.0_wp),-0.15_wp)
            else
                call expect_close(dsdp(i,2,2),0.0_wp)
            end if
        end do
    end subroutine test_large_elastic_limit

    subroutine test_rejected_trials()
        type(material_t) :: material
        real(wp) :: x(ninput,npoints), stress(nmandel,npoints)
        real(wp) :: wrong(ninput-1,npoints)
        integer :: status, i

        call material%init(isotropic_thermoelastic_t(0.0_wp,shear,expansion,reference_t), &
                           npoints,status)
        call expect_status(status,cf_invalid_argument)
        call material%init(isotropic_thermoelastic_t(bulk,0.0_wp,expansion,reference_t), &
                           npoints,status)
        call expect_status(status,cf_invalid_argument)
        call material%init(isotropic_thermoelastic_t(bulk,shear,expansion,-1.0_wp), &
                           npoints,status)
        call expect_status(status,cf_invalid_argument)
        call material%init(isotropic_thermoelastic_t(bulk,shear,-expansion,reference_t), &
                           npoints,status)
        call expect_status(status,cf_ok)
        x = 0.0_wp
        x(ninput,:) = reference_t
        x(ninput,1) = reference_t+15.0_wp
        call material%eval(x,stress,status=status)
        call expect_status(status,cf_ok)
        do i = 1, ndim
            call expect_close(stress(i,1),0.45_wp)
        end do
        x(ninput,2) = -1.0_wp
        call material%eval(x,stress,status=status)
        call expect_status(status,cf_evaluation_failed)
        if (material%last_failure_point() /= 2) error stop "Negative kelvin accepted"
        call material%commit(status)
        call expect_status(status,cf_no_trial)
        x(ninput,2) = ieee_value(0.0_wp,ieee_quiet_nan)
        call material%eval(x,stress,status=status)
        call expect_status(status,cf_evaluation_failed)
        x(ninput,2) = reference_t
        call material%eval(x,stress,status=status)
        call expect_status(status,cf_ok)
        wrong = 0.0_wp
        call material%eval(wrong,stress,status=status)
        call expect_status(status,cf_invalid_shape)
        call material%init(isotropic_thermoelastic_t(bulk,0.0_wp,expansion,reference_t), &
                           npoints,status)
        call expect_status(status,cf_invalid_argument)
        call material%eval(x,stress,status=status)
        call expect_status(status,cf_ok)
        call material%commit(status)
        call expect_status(status,cf_ok)
    end subroutine test_rejected_trials

end program test_isotropic_thermoelastic

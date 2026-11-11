! SPDX-License-Identifier: Apache-2.0

program test_kelvin_voigt
    !! Test convention: moduli and stress in Pa, viscosities in Pa*s and
    !! dt in s. Strain is dimensionless.
    !!
    !! Independent Kelvin--Voigt elastic/viscous references and path checks.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use constituf, only: wp, material_t, cf_ok, cf_invalid_argument, &
                        cf_invalid_time, cf_invalid_shape, cf_evaluation_failed, &
                        cf_no_trial, cf_unsupported
    use constituf_kelvin_voigt, only: kelvin_voigt_t
    implicit none (type, external)

    integer, parameter :: ndim = 3, nmandel = 6, npoints = 2, nparameters = 4
    real(wp), parameter :: bulk = 100.0_wp, shear = 40.0_wp
    real(wp), parameter :: bulk_viscosity = 25.0_wp, shear_viscosity = 10.0_wp

    call test_reference_response()
    call test_time_scale_references()
    call test_consistent_tangent()
    call test_history_sensitivities()
    call test_invalid_trials()

contains

    subroutine expect_status(actual, expected)
        integer, intent(in) :: actual, expected

        if (actual /= expected) then
            print *, "Expected status", expected, "got", actual
            error stop "Kelvin-Voigt API status mismatch"
        end if
    end subroutine expect_status

    subroutine close(actual, expected)
        real(wp), intent(in) :: actual, expected
        real(wp) :: tolerance

        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) then
            error stop "Nonfinite Kelvin-Voigt comparison"
        end if
        tolerance = 512.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(expected))
        if (abs(actual-expected) > tolerance) then
            print *, "Expected", expected, "got", actual
            error stop "Kelvin-Voigt numerical mismatch"
        end if
    end subroutine close

    subroutine test_reference_response()
        type(material_t) :: mat, single
        real(wp) :: initial(nmandel,npoints), x(nmandel,npoints)
        real(wp) :: y(nmandel,npoints), tangent(nmandel,nmandel,npoints)
        real(wp) :: one(nmandel,1), rate(nmandel)
        real(wp) :: dsdp(nmandel,nparameters,npoints), previous(nmandel,npoints)
        real(wp) :: state(0,npoints), sensitivity_state(0,nparameters,npoints)
        real(wp) :: elastic_stress(nmandel)
        real(wp) :: viscosity_power
        integer :: i, status

        initial = 0.0_wp
        call mat%init(kelvin_voigt_t(bulk,shear,bulk_viscosity,shear_viscosity), &
                      npoints,status)
        call expect_status(status,cf_invalid_argument)
        call mat%init(kelvin_voigt_t(bulk,shear,bulk_viscosity,shear_viscosity), &
                      npoints,status,initial, &
                      sensitivity_parameters=[character(len=20) :: "bulk_modulus", &
                          "shear_modulus","bulk_viscosity","shear_viscosity"])
        call expect_status(status,cf_ok)
        call single%init(kelvin_voigt_t(bulk,shear,bulk_viscosity,shear_viscosity), &
                         1,status,initial(:,1:1))
        call expect_status(status,cf_ok)
        if (mat%input_size() /= nmandel .or. mat%output_size() /= nmandel) then
            error stop "Kelvin-Voigt Mandel dimensions"
        end if
        if (mat%state_size() /= 0 .or. mat%state_storage_bytes() /= 0) then
            error stop "Kelvin-Voigt has no extra constitutive state"
        end if
        if (.not. mat%field_matches("input","strain.small","mandel", &
            "small_strain","spatial",[ndim,ndim])) then
            error stop "Kelvin-Voigt strain metadata mismatch"
        end if
        if (.not. mat%field_matches("output","stress.cauchy","mandel", &
            "cauchy_stress","spatial",[ndim,ndim])) then
            error stop "Kelvin-Voigt Cauchy stress metadata mismatch"
        end if

        x(:,1) = [0.02_wp,-0.01_wp,-0.01_wp,0.008_wp,-0.006_wp,0.005_wp]
        x(:,2) = [0.01_wp,0.01_wp,0.01_wp,0.0_wp,0.0_wp,0.0_wp]
        call mat%eval(x,y,tangent,status,dt=0.5_wp,dsdp=dsdp)
        call expect_status(status,cf_ok)
        call single%eval(x(:,1:1),one,status=status,dt=0.5_wp)
        call expect_status(status,cf_ok)
        do i = 1, nmandel
            call close(y(i,1),120.0_wp*x(i,1))
            call close(y(i,1),one(i,1))
            call close(dsdp(i,1,1),0.0_wp)
            call close(dsdp(i,2,1),2.0_wp*x(i,1))
            call close(dsdp(i,3,1),0.0_wp)
            call close(dsdp(i,4,1),4.0_wp*x(i,1))
            call close(dsdp(i,1,2),3.0_wp*x(i,2))
            call close(dsdp(i,2,2),0.0_wp)
            call close(dsdp(i,3,2),6.0_wp*x(i,2))
            call close(dsdp(i,4,2),0.0_wp)
        end do
        do i = 1, ndim
            call close(y(i,2),4.5_wp)
        end do
        call close(tangent(1,1,1),230.0_wp)
        call close(tangent(1,2,1),110.0_wp)
        call close(tangent(4,4,1),120.0_wp)

        rate = (x(:,1)-initial(:,1))/0.5_wp
        elastic_stress = 2.0_wp*shear*x(:,1)
        viscosity_power = dot_product(y(:,1)-elastic_stress,rate)
        if (.not. ieee_is_finite(viscosity_power) .or. viscosity_power < 0.0_wp) then
            error stop "Negative Kelvin-Voigt viscous dissipation"
        end if

        call mat%export_committed(state,status,previous_input=previous,sensitivity_state=sensitivity_state)
        call expect_status(status,cf_ok)
        if (.not. all(ieee_is_finite(previous)) .or. maxval(abs(previous-initial)) > 0.0_wp) then
            error stop "Kelvin-Voigt trial changed previous input"
        end if
        call mat%commit(status)
        call expect_status(status,cf_ok)
        call mat%export_committed(state,status,previous_input=previous,sensitivity_state=sensitivity_state)
        call expect_status(status,cf_ok)
        if (.not. all(ieee_is_finite(previous)) .or. maxval(abs(previous-x)) > 0.0_wp) then
            error stop "Kelvin-Voigt commit did not accept previous input"
        end if
        call mat%eval(x,y,tangent,status,dt=0.5_wp,dsdp=dsdp)
        call expect_status(status,cf_ok)
        do i = 1, nmandel
            call close(y(i,1),2.0_wp*shear*x(i,1))
            call close(dsdp(i,1,1),0.0_wp)
            call close(dsdp(i,2,1),2.0_wp*x(i,1))
            call close(dsdp(i,3,1),0.0_wp)
            call close(dsdp(i,4,1),0.0_wp)
            call close(dsdp(i,1,2),3.0_wp*x(i,2))
            call close(dsdp(i,2,2),0.0_wp)
            call close(dsdp(i,3,2),0.0_wp)
            call close(dsdp(i,4,2),0.0_wp)
        end do
        ! A zero rate removes viscous stress, but the endpoint tangent still
        ! includes viscosity because the committed previous strain is fixed.
        call close(tangent(1,1,1),230.0_wp)
        call close(tangent(1,2,1),110.0_wp)
        call close(tangent(4,4,1),120.0_wp)
        do i = 1, ndim
            call close(y(i,2),3.0_wp)
        end do
        call mat%eval(0.5_wp*x,y,status=status,dt=0.5_wp,dsdp=dsdp)
        call expect_status(status,cf_ok)
        call mat%revert()
        call mat%commit(status)
        call expect_status(status,cf_no_trial)
        call mat%export_committed(state,status,previous_input=previous,sensitivity_state=sensitivity_state)
        call expect_status(status,cf_ok)
        if (.not. all(ieee_is_finite(previous)) .or. maxval(abs(previous-x)) > 0.0_wp) then
            error stop "Kelvin-Voigt revert changed previous input"
        end if

        ! Zero viscosities recover elasticity for every positive step length.
        call mat%init(kelvin_voigt_t(bulk,shear,0.0_wp,0.0_wp), &
                      npoints,status,initial)
        call expect_status(status,cf_ok)
        call mat%eval(x,y,status=status,dt=0.2_wp)
        call expect_status(status,cf_ok)
        do i = 1, nmandel
            call close(y(i,1),2.0_wp*shear*x(i,1))
        end do
        call close(y(1,2),3.0_wp)
    end subroutine test_reference_response

    subroutine test_time_scale_references()
        !! These finite responses must not overflow from an intermediate 2*eta,
        !! (4/3)*eta or 2/dt. The model is affine, so both limits have exact references.
        type(material_t) :: mat
        real(wp) :: initial(nmandel,1), x(nmandel,1), y(nmandel,1)
        real(wp) :: tangent(nmandel,nmandel,1), dsdp(nmandel,nparameters,1)
        real(wp) :: time_scale, time_step
        integer :: status, i

        initial = 0.0_wp
        time_scale = huge(1.0_wp)/16.0_wp
        call mat%init(kelvin_voigt_t(bulk,shear,0.0_wp,14.0_wp*time_scale), &
            1,status,initial,sensitivity_parameters=[character(len=20) :: &
                "bulk_modulus","shear_modulus","bulk_viscosity","shear_viscosity"])
        call expect_status(status,cf_ok)
        x = 0.0_wp
        x(nmandel,1) = 0.125_wp
        call mat%eval(x,y,tangent,status,dt=0.5_wp*time_scale,dsdp=dsdp)
        call expect_status(status,cf_ok)
        do i = 1,nmandel
            call close(y(i,1),136.0_wp*x(i,1))
            call close(dsdp(i,1,1),0.0_wp)
            call close(dsdp(i,2,1),2.0_wp*x(i,1))
            call close(dsdp(i,3,1),0.0_wp)
            call close(dsdp(i,4,1)*time_scale,4.0_wp*x(i,1))
        end do
        call close(tangent(1,1,1),572.0_wp/3.0_wp)
        call close(tangent(1,2,1),164.0_wp/3.0_wp)
        call close(tangent(nmandel,nmandel,1),136.0_wp)

        call mat%init(kelvin_voigt_t(bulk,shear,0.0_wp,0.0_wp), &
            1,status,initial,sensitivity_parameters=[character(len=20) :: &
                "bulk_modulus","shear_modulus","bulk_viscosity","shear_viscosity"])
        call expect_status(status,cf_ok)
        time_step = 1.0_wp/huge(1.0_wp)
        x = 0.0_wp
        x(1:ndim,1) = time_step
        x(nmandel,1) = time_step
        call mat%eval(x,y,tangent,status,dt=time_step,dsdp=dsdp)
        call expect_status(status,cf_ok)
        do i = 1,ndim
            call close(y(i,1)/time_step,300.0_wp)
            call close(dsdp(i,1,1)/time_step,3.0_wp)
            call close(dsdp(i,2,1),0.0_wp)
            call close(dsdp(i,3,1),3.0_wp)
            call close(dsdp(i,4,1),0.0_wp)
        end do
        call close(y(nmandel,1)/time_step,80.0_wp)
        call close(dsdp(nmandel,1,1),0.0_wp)
        call close(dsdp(nmandel,2,1)/time_step,2.0_wp)
        call close(dsdp(nmandel,3,1),0.0_wp)
        call close(dsdp(nmandel,4,1),2.0_wp)
        call close(tangent(1,1,1),460.0_wp/3.0_wp)
        call close(tangent(1,2,1),220.0_wp/3.0_wp)
        call close(tangent(nmandel,nmandel,1),80.0_wp)
    end subroutine test_time_scale_references

    subroutine test_consistent_tangent()
        type(material_t) :: mat
        real(wp) :: initial(nmandel,1), x(nmandel,1)
        real(wp) :: y(nmandel,1), tangent(nmandel,nmandel,1)
        real(wp) :: changed(nmandel,1), yp(nmandel,1), ym(nmandel,1)
        real(wp) :: numerical, h, time_step, roundoff
        integer :: status, i, j

        initial = 0.0_wp
        call mat%init(kelvin_voigt_t(bulk,shear,bulk_viscosity,shear_viscosity), &
                      1,status,initial)
        call expect_status(status,cf_ok)
        x(:,1) = [0.013_wp,-0.004_wp,0.007_wp,0.002_wp,-0.003_wp,0.004_wp]
        call mat%eval(x,y,status=status,dt=0.2_wp)
        call expect_status(status,cf_ok)
        call mat%commit(status)
        call expect_status(status,cf_ok)
        x(:,1) = [0.020_wp,0.002_wp,-0.008_wp,0.006_wp,0.001_wp,-0.005_wp]
        time_step = 0.35_wp
        call mat%eval(x,y,tangent,status,dt=time_step)
        call expect_status(status,cf_ok)
        if (.not. all(ieee_is_finite(tangent))) error stop "Nonfinite Kelvin-Voigt tangent"
        ! This response is affine in strain: a finite central step has no
        ! truncation error, and avoids cancellation at every supported kind.
        h = 0.25_wp*maxval(abs(x))
        do j = 1, nmandel
            changed = x
            changed(j,1) = changed(j,1)+h
            call mat%eval(changed,yp,status=status,dt=time_step)
            call expect_status(status,cf_ok)
            changed(j,1) = changed(j,1)-2.0_wp*h
            call mat%eval(changed,ym,status=status,dt=time_step)
            call expect_status(status,cf_ok)
            roundoff = 32.0_wp*epsilon(1.0_wp)*maxval(abs(yp)+abs(ym))/(2.0_wp*h)
            if (.not. all(ieee_is_finite(yp)) .or. .not. all(ieee_is_finite(ym)) .or. &
                .not. ieee_is_finite(roundoff)) error stop "Nonfinite Kelvin-Voigt tangent secant"
            do i = 1, nmandel
                numerical = (yp(i,1)-ym(i,1))/(2.0_wp*h)
                if (.not. ieee_is_finite(numerical) .or. abs(numerical-tangent(i,j,1)) > roundoff) then
                    error stop "Inconsistent Kelvin-Voigt algorithmic tangent"
                end if
                call close(tangent(i,j,1),tangent(j,i,1))
            end do
        end do
        call mat%revert()
        call mat%eval(x,y,tangent,status,dt=time_step)
        call expect_status(status,cf_ok)
    end subroutine test_consistent_tangent

    subroutine test_history_sensitivities()
        !! Perturb all parameters over separately accepted prescribed histories.
        type(material_t) :: tracked, plain
        type(material_t) :: positive(nparameters), negative(nparameters)
        real(wp) :: initial(nmandel,npoints), x(nmandel,npoints)
        real(wp) :: y(nmandel,npoints), reference(nmandel,npoints)
        real(wp) :: tangent(nmandel,nmandel,npoints)
        real(wp) :: dsdp(nmandel,nparameters,npoints)
        real(wp) :: prior_dsdp(nmandel,nparameters,npoints)
        real(wp) :: yp(nmandel,npoints), ym(nmandel,npoints)
        real(wp) :: parameters(nparameters), plus_parameters(nparameters)
        real(wp) :: minus_parameters(nparameters), steps(nparameters)
        real(wp) :: time_step, numerical, roundoff
        integer :: status, k, id, ipoint, i, step

        initial = 0.0_wp
        initial(1,1) = 0.003_wp
        initial(2,2) = -0.002_wp
        parameters = [bulk,shear,bulk_viscosity,shear_viscosity]
        call tracked%init(kelvin_voigt_t(bulk,shear,bulk_viscosity,shear_viscosity), &
            npoints,status,initial,sensitivity_parameters=[character(len=20) :: &
                "shear_viscosity","bulk_viscosity","shear_modulus","bulk_modulus"])
        call expect_status(status,cf_ok)
        call plain%init(kelvin_voigt_t(bulk,shear,bulk_viscosity,shear_viscosity), &
                        npoints,status,initial)
        call expect_status(status,cf_ok)
        if (tracked%sensitivity_size() /= nparameters) error stop "Selected parameters"
        if (tracked%state_storage_bytes() /= 0 .or. plain%state_storage_bytes() /= 0) then
            error stop "Kelvin-Voigt unexpectedly allocated real history"
        end if

        do k = 1, nparameters
            id = nparameters+1-k
            ! Stress is also affine in each parameter at prescribed history.
            steps(k) = 0.25_wp*parameters(id)
            plus_parameters = parameters
            minus_parameters = parameters
            plus_parameters(id) = plus_parameters(id)+steps(k)
            minus_parameters(id) = minus_parameters(id)-steps(k)
            call positive(k)%init(kelvin_voigt_t(plus_parameters(1),plus_parameters(2), &
                plus_parameters(3),plus_parameters(4)),npoints,status,initial)
            call expect_status(status,cf_ok)
            call negative(k)%init(kelvin_voigt_t(minus_parameters(1),minus_parameters(2), &
                minus_parameters(3),minus_parameters(4)),npoints,status,initial)
            call expect_status(status,cf_ok)
        end do

        do step = 1, 4
            time_step = 0.2_wp+0.15_wp*real(step,wp)
            x = initial
            x(1,1) = x(1,1)+0.01_wp*real(step,wp)
            x(2,1) = x(2,1)-0.007_wp*real(step,wp)
            x(4,1) = 0.004_wp*real(3-step,wp)
            x(1,2) = x(1,2)+0.012_wp*real(3-step,wp)
            x(2,2) = x(2,2)+0.008_wp*real(step,wp)
            x(5,2) = 0.003_wp*real(step-2,wp)
            call tracked%eval(x,y,tangent,status,dt=time_step,dsdp=dsdp)
            call expect_status(status,cf_ok)
            if (.not. all(ieee_is_finite(dsdp))) error stop "Nonfinite Kelvin-Voigt parameter derivatives"
            call plain%eval(x,reference,status=status,dt=time_step)
            call expect_status(status,cf_ok)
            do ipoint = 1, npoints
                do i = 1, nmandel
                    call close(y(i,ipoint),reference(i,ipoint))
                end do
            end do

            do k = 1, nparameters
                call positive(k)%eval(x,yp,status=status,dt=time_step)
                call expect_status(status,cf_ok)
                call negative(k)%eval(x,ym,status=status,dt=time_step)
                call expect_status(status,cf_ok)
                if (.not. all(ieee_is_finite(yp)) .or. .not. all(ieee_is_finite(ym))) then
                    error stop "Nonfinite Kelvin-Voigt parameter secant response"
                end if
                do ipoint = 1, npoints
                    roundoff = 32.0_wp*epsilon(1.0_wp)* &
                        maxval(abs(yp(:,ipoint))+abs(ym(:,ipoint)))/(2.0_wp*steps(k))
                    if (.not. ieee_is_finite(roundoff)) error stop "Nonfinite Kelvin-Voigt sensitivity bound"
                    do i = 1, nmandel
                        numerical = (yp(i,ipoint)-ym(i,ipoint))/(2.0_wp*steps(k))
                        if (.not. ieee_is_finite(numerical) .or. abs(numerical-dsdp(i,k,ipoint)) > roundoff) then
                            error stop "Full-history Kelvin-Voigt sensitivity mismatch"
                        end if
                    end do
                end do
            end do

            if (step == 2) then
                prior_dsdp = dsdp
                call tracked%revert()
                call plain%revert()
                do k = 1, nparameters
                    call positive(k)%revert()
                    call negative(k)%revert()
                end do
                call tracked%eval(x,y,tangent,status,dt=time_step,dsdp=dsdp)
                call expect_status(status,cf_ok)
                do k = 1, nparameters
                    do ipoint = 1, npoints
                        do i = 1, nmandel
                            call close(dsdp(i,k,ipoint),prior_dsdp(i,k,ipoint))
                        end do
                    end do
                end do
                call plain%eval(x,reference,status=status,dt=time_step)
                call expect_status(status,cf_ok)
                do k = 1, nparameters
                    call positive(k)%eval(x,yp,status=status,dt=time_step)
                    call expect_status(status,cf_ok)
                    call negative(k)%eval(x,ym,status=status,dt=time_step)
                    call expect_status(status,cf_ok)
                end do
            end if

            call tracked%commit(status)
            call expect_status(status,cf_ok)
            call plain%commit(status)
            call expect_status(status,cf_ok)
            do k = 1, nparameters
                call positive(k)%commit(status)
                call expect_status(status,cf_ok)
                call negative(k)%commit(status)
                call expect_status(status,cf_ok)
            end do
        end do
        call tracked%commit(status)
        call expect_status(status,cf_no_trial)
    end subroutine test_history_sensitivities

    subroutine test_invalid_trials()
        type(material_t) :: mat, tracked
        real(wp) :: initial(nmandel,npoints), x(nmandel,npoints)
        real(wp) :: stress(nmandel,npoints), wrong(nmandel-1,npoints)
        real(wp) :: dsdp(nmandel,1,npoints)
        integer :: status

        initial = 0.0_wp
        call mat%init(kelvin_voigt_t(0.0_wp,shear,bulk_viscosity,shear_viscosity), &
                      npoints,status,initial)
        call expect_status(status,cf_invalid_argument)
        call mat%init(kelvin_voigt_t(bulk,shear,-1.0_wp,shear_viscosity), &
                      npoints,status,initial)
        call expect_status(status,cf_invalid_argument)
        call mat%init(kelvin_voigt_t(bulk,shear,bulk_viscosity,-1.0_wp), &
                      npoints,status,initial)
        call expect_status(status,cf_invalid_argument)
        call mat%init(kelvin_voigt_t(bulk,shear,bulk_viscosity,shear_viscosity), &
                      npoints,status,initial)
        call expect_status(status,cf_ok)

        x = 0.0_wp
        x(1,:) = 0.02_wp
        call mat%eval(x,stress,status=status)
        call expect_status(status,cf_invalid_time)
        call mat%eval(x,stress,status=status,dt=0.0_wp)
        call expect_status(status,cf_invalid_time)
        x(1,2) = ieee_value(0.0_wp,ieee_quiet_nan)
        call mat%eval(x,stress,status=status,dt=0.5_wp)
        call expect_status(status,cf_evaluation_failed)
        if (mat%last_failure_point() /= 2) error stop "Invalid Kelvin-Voigt point"
        call mat%commit(status)
        call expect_status(status,cf_no_trial)
        x(1,2) = 0.02_wp
        call mat%eval(x,stress,status=status,dt=0.5_wp)
        call expect_status(status,cf_ok)
        wrong = 0.0_wp
        call mat%eval(wrong,stress,status=status,dt=0.5_wp)
        call expect_status(status,cf_invalid_shape)
        call mat%init(kelvin_voigt_t(bulk,-1.0_wp,bulk_viscosity,shear_viscosity), &
                      npoints,status,initial)
        call expect_status(status,cf_invalid_argument)
        call mat%eval(x,stress,status=status,dt=0.5_wp)
        call expect_status(status,cf_ok)
        call mat%commit(status)
        call expect_status(status,cf_ok)
        call mat%eval(x,stress,status=status,dt=0.5_wp,dsdp=dsdp)
        call expect_status(status,cf_unsupported)

        call tracked%init(kelvin_voigt_t(bulk,shear,bulk_viscosity,shear_viscosity), &
            npoints,status,initial,sensitivity_parameters=[character(len=20) :: &
                "shear_viscosity"])
        call expect_status(status,cf_ok)
        call tracked%eval(x,stress,status=status,dt=0.5_wp)
        call expect_status(status,cf_invalid_argument)
        call tracked%eval(x,stress,status=status,dt=0.5_wp,dsdp=dsdp)
        call expect_status(status,cf_ok)
    end subroutine test_invalid_trials

end program test_kelvin_voigt

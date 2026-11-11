! SPDX-License-Identifier: Apache-2.0

program test_j2_plasticity
    !! Test convention: stress and stress-like material parameters in Pa.
    !! Strain and deformation gradients are dimensionless. Explicit rescaling
    !! cases state their factors.
    !!
    !! Independent J2 stress, tangent, dissipation and path-derivative tests.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use constituf, only: wp, lk, material_t, verify_tangent, cf_ok, cf_invalid_argument, &
                        cf_invalid_shape, cf_evaluation_failed, cf_no_trial, &
                        cf_nondifferentiable, cf_unsupported
    use constituf_j2_plasticity, only: j2_plasticity_t
    implicit none (type, external)

    integer, parameter :: nmandel = 6, nstate = nmandel+1
    integer, parameter :: npoints = 2, nparameters = 4
    real(wp), parameter :: bulk = 100.0_wp, shear = 40.0_wp
    real(wp), parameter :: initial_yield = 5.0_wp, hardening = 10.0_wp

    call test_schema_and_elasticity()
    call test_radial_return_and_tangent()
    call test_full_history_sensitivities()
    call test_reversed_plastic_state()
    call test_failure_and_yield_switch()

contains

    subroutine expect_status(value, expected)
        integer, intent(in) :: value, expected

        if (value /= expected) then
            print *, "Expected status", expected, "got", value
            error stop "J2 plasticity status mismatch"
        end if
    end subroutine expect_status

    subroutine expect_close(value, expected)
        real(wp), intent(in) :: value, expected

        if (.not. ieee_is_finite(value) .or. .not. ieee_is_finite(expected)) then
            error stop "Nonfinite J2 plasticity reference comparison"
        end if
        if (abs(value-expected) > &
            512.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(expected))) then
            print *, "Expected", expected, "got", value
            error stop "J2 plasticity numerical mismatch"
        end if
    end subroutine expect_close

    subroutine test_schema_and_elasticity()
        type(material_t) :: mat
        real(wp) :: x(nmandel,npoints), stress(nmandel,npoints)
        real(wp) :: tangent(nmandel,nmandel,npoints)
        integer(lk) :: expected_bytes
        integer :: status, ipoint
        integer :: empty_shape(0)

        call mat%init(j2_plasticity_t(0.0_wp,shear,initial_yield,hardening),npoints,status)
        call expect_status(status, cf_invalid_argument)
        call mat%init(j2_plasticity_t(bulk,0.0_wp,initial_yield,hardening),npoints,status)
        call expect_status(status, cf_invalid_argument)
        call mat%init(j2_plasticity_t(bulk,shear,0.0_wp,hardening),npoints,status)
        call expect_status(status, cf_invalid_argument)
        call mat%init(j2_plasticity_t(bulk,shear,initial_yield,-1.0_wp),npoints,status)
        call expect_status(status, cf_invalid_argument)
        call mat%init(j2_plasticity_t(bulk,shear,initial_yield,hardening),npoints,status)
        call expect_status(status, cf_ok)

        if (mat%input_size() /= nmandel .or. mat%output_size() /= nmandel) then
            error stop "J2 must have symmetric 3D Mandel fields"
        end if
        if (mat%state_size() /= nstate) error stop "J2 real state size mismatch"
        if (.not. mat%field_matches("state","plastic.strain","mandel", &
                                    "plastic_strain","spatial",[3,3])) then
            error stop "Plastic strain metadata mismatch"
        end if
        if (.not. mat%field_matches("state","plastic.accumulated","scalar", &
                                    "equivalent_plastic_strain","material", empty_shape)) then
            error stop "Equivalent plastic strain metadata mismatch"
        end if
        if (.not. mat%field_matches("output","stress.cauchy","mandel", &
                                    "cauchy_stress","spatial",[3,3])) then
            error stop "Cauchy stress metadata mismatch"
        end if
        expected_bytes = 2_lk*int(nstate,lk)*int(npoints,lk)* &
                         int(storage_size(0.0_wp)/8,lk)
        if (mat%state_storage_bytes() /= expected_bytes) then
            error stop "Incorrect retained plastic history storage"
        end if
        if (mat%sensitivity_size() /= 0) error stop "Derivatives allocated without request"

        x = 0.0_wp
        x(1,:) = 0.01_wp
        x(2,:) = 0.01_wp
        x(3,:) = 0.01_wp
        call mat%eval(x,stress,tangent,status)
        call expect_status(status,cf_ok)
        do ipoint = 1,npoints
            call expect_close(stress(1,ipoint),3.0_wp*bulk*0.01_wp)
            call expect_close(stress(2,ipoint),3.0_wp*bulk*0.01_wp)
            call expect_close(stress(3,ipoint),3.0_wp*bulk*0.01_wp)
            if (maxval(abs(stress(4:6,ipoint))) > 1.0e-14_wp) then
                error stop "Hydrostatic load created deviatoric stress"
            end if
            call expect_close(tangent(1,1,ipoint),bulk+4.0_wp*shear/3.0_wp)
            call expect_close(tangent(1,2,ipoint),bulk-2.0_wp*shear/3.0_wp)
            call expect_close(tangent(4,4,ipoint),2.0_wp*shear)
        end do
        call mat%commit(status)
        call expect_status(status,cf_ok)
    end subroutine test_schema_and_elasticity

    subroutine test_radial_return_and_tangent()
        type(material_t) :: mat
        real(wp) :: x(nmandel,1), stress(nmandel,1), tangent(nmandel,nmandel,1)
        real(wp) :: gamma, reference, error, yield_value
        real(wp) :: stress_deviator(nmandel)
        integer :: status, i, j

        call mat%init(j2_plasticity_t(bulk,shear,initial_yield,hardening),1,status)
        call expect_status(status,cf_ok)
        x = 0.0_wp
        x(1,1) = 0.2_wp
        gamma = (2.0_wp*shear*x(1,1)-initial_yield)/(3.0_wp*shear+hardening)
        call mat%eval(x,stress,tangent,status)
        call expect_status(status,cf_ok)
        reference = bulk*x(1,1)+4.0_wp*shear*x(1,1)/3.0_wp- &
                    2.0_wp*shear*gamma
        call expect_close(stress(1,1),reference)
        reference = bulk*x(1,1)-2.0_wp*shear*x(1,1)/3.0_wp+shear*gamma
        call expect_close(stress(2,1),reference)
        call expect_close(stress(3,1),reference)
        call mat%commit(status)
        call expect_status(status,cf_ok)

        ! Elastic unloading preserves plastic history and recovers elastic stiffness.
        x(1,1) = 0.15_wp
        call mat%eval(x,stress,tangent,status)
        call expect_status(status,cf_ok)
        reference = bulk*x(1,1)+2.0_wp*shear*(2.0_wp*x(1,1)/3.0_wp-gamma)
        call expect_close(stress(1,1),reference)
        call expect_close(tangent(1,1,1),bulk+4.0_wp*shear/3.0_wp)
        call mat%revert()
        call mat%eval(x,stress,tangent,status)
        call expect_status(status,cf_ok)
        call expect_close(stress(1,1),reference)
        call mat%commit(status)
        call expect_status(status,cf_ok)

        ! Reload into the plastic branch from a nonzero committed plastic state.
        x(1,1) = 0.25_wp
        x(4,1) = 0.02_wp
        call mat%eval(x,stress,tangent,status)
        call expect_status(status,cf_ok)
        do j = 1,nmandel
            do i = 1,nmandel
                call expect_close(tangent(i,j,1),tangent(j,i,1))
            end do
        end do
        call verify_tangent(mat,x,error,status)
        call expect_status(status,cf_ok)
        if (.not. ieee_is_finite(error) .or. error > 1.0_wp) then
            error stop "Inconsistent algorithmic J2 return tangent"
        end if

        ! Perfect plasticity limits the returned von Mises stress to yield.
        block
            type(material_t) :: perfect
            call perfect%init(j2_plasticity_t(bulk,shear,initial_yield,0.0_wp),1,status)
            call expect_status(status,cf_ok)
            x = 0.0_wp
            x(1,1) = 0.3_wp
            call perfect%eval(x,stress,status=status)
            call expect_status(status,cf_ok)
            stress_deviator = stress(:,1)
            stress_deviator(1:3) = stress_deviator(1:3) - sum(stress(1:3,1))/3.0_wp
            yield_value = sqrt(1.5_wp*dot_product(stress_deviator,stress_deviator))
            call expect_close(yield_value,initial_yield)
        end block
    end subroutine test_radial_return_and_tangent

    subroutine test_full_history_sensitivities()
        !! Exact derivatives compared against independently perturbed paths.
        type(material_t) :: tracked, plain, single
        type(material_t) :: positive(nparameters,2), negative(nparameters,2)
        real(wp) :: x(nmandel,npoints), stress(nmandel,npoints)
        real(wp) :: reference(nmandel,npoints), single_stress(nmandel,1)
        real(wp) :: tangent(nmandel,nmandel,npoints)
        real(wp) :: dsdp(nmandel,nparameters,npoints)
        real(wp) :: prior_dsdp(nmandel,nparameters,npoints)
        real(wp) :: plus_stress(nmandel,npoints), minus_stress(nmandel,npoints)
        real(wp) :: parameters(nparameters), plus_parameters(nparameters)
        real(wp) :: minus_parameters(nparameters), perturbations(nparameters,2)
        real(wp) :: numerical(nmandel,npoints,2), tolerance(nmandel,npoints)
        integer(lk) :: history_extra
        integer :: status, step, k, i, j, id, refinement

        parameters = [bulk,shear,initial_yield,hardening]
        call tracked%init(j2_plasticity_t(bulk,shear,initial_yield,hardening), &
            npoints,status,sensitivity_parameters=[character(len=20) :: &
                "hardening_modulus","yield_stress","shear_modulus","bulk_modulus"])
        call expect_status(status,cf_ok)
        call plain%init(j2_plasticity_t(bulk,shear,initial_yield,hardening),npoints,status)
        call expect_status(status,cf_ok)
        call single%init(j2_plasticity_t(bulk,shear,initial_yield,hardening),1,status)
        call expect_status(status,cf_ok)
        history_extra = 2_lk*int(nstate,lk)*int(npoints,lk)*int(nparameters,lk)* &
                        int(storage_size(0.0_wp)/8,lk)
        if (tracked%state_storage_bytes()-plain%state_storage_bytes() /= &
            history_extra) error stop "Incorrect plastic sensitivity state memory"
        if (plain%sensitivity_size() /= 0) error stop "Plain J2 tracked sensitivity"

        do k = 1,nparameters
            id = nparameters+1-k
            perturbations(k,1) = epsilon(1.0_wp)**(1.0_wp/3.0_wp)*max(1.0_wp,parameters(id))
            perturbations(k,2) = 0.5_wp*perturbations(k,1)
            do refinement = 1,2
                plus_parameters = parameters
                minus_parameters = parameters
                plus_parameters(id) = plus_parameters(id)+perturbations(k,refinement)
                minus_parameters(id) = minus_parameters(id)-perturbations(k,refinement)
                call positive(k,refinement)%init(j2_plasticity_t(plus_parameters(1),plus_parameters(2), &
                    plus_parameters(3),plus_parameters(4)),npoints,status)
                call expect_status(status,cf_ok)
                call negative(k,refinement)%init(j2_plasticity_t(minus_parameters(1),minus_parameters(2), &
                    minus_parameters(3),minus_parameters(4)),npoints,status)
                call expect_status(status,cf_ok)
            end do
        end do

        do step = 1,5
            x = 0.0_wp
            select case (step)
            case (1)
                x(1,1) = 0.02_wp
                x(1,2) = -0.01_wp
            case (2)
                x(1,1) = 0.14_wp
                x(1,2) = -0.10_wp
            case (3)
                x(1,1) = 0.09_wp
                x(1,2) = -0.045_wp
            case (4)
                x(1,1) = 0.21_wp
                x(1,2) = -0.14_wp
            case (5)
                x(1,1) = 0.16_wp
                x(1,2) = -0.08_wp
            case default
                error stop "Invalid path step"
            end select
            x(4,1) = 0.007_wp*real(step-2,wp)
            x(5,2) = 0.005_wp*real(step-2,wp)
            call tracked%eval(x,stress,tangent,status,dsdp=dsdp)
            call expect_status(status,cf_ok)
            call plain%eval(x,reference,status=status)
            call expect_status(status,cf_ok)
            do i = 1,npoints
                do j = 1,nmandel
                    call expect_close(stress(j,i),reference(j,i))
                end do
            end do
            if (step == 1) then
                call single%eval(x(:,1:1),single_stress,status=status)
                call expect_status(status,cf_ok)
                do j = 1,nmandel
                    call expect_close(stress(j,1),single_stress(j,1))
                end do
            end if
            do k = 1,nparameters
                ! Two independently committed paths estimate centered-difference truncation.
                do refinement = 1,2
                    call positive(k,refinement)%eval(x,plus_stress,status=status)
                    call expect_status(status,cf_ok)
                    call negative(k,refinement)%eval(x,minus_stress,status=status)
                    call expect_status(status,cf_ok)
                    numerical(:,:,refinement) = (plus_stress-minus_stress)/ &
                        (2.0_wp*perturbations(k,refinement))
                end do
                tolerance = abs(numerical(:,:,2)-numerical(:,:,1))/3.0_wp + &
                    64.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(plus_stress),abs(minus_stress))/perturbations(k,2)
                if (.not. all(ieee_is_finite(numerical)) .or. .not. all(ieee_is_finite(dsdp(:,k,:)))) then
                    error stop "Nonfinite accepted-history J2 parameter derivative"
                end if
                if (any(abs(numerical(:,:,2)-dsdp(:,k,:)) > tolerance)) then
                    error stop "Incorrect accepted-history J2 parameter derivative"
                end if
            end do

            if (step == 3) then
                prior_dsdp = dsdp
                call tracked%revert()
                call plain%revert()
                do k = 1,nparameters
                    do refinement = 1,2
                        call positive(k,refinement)%revert()
                        call negative(k,refinement)%revert()
                    end do
                end do
                call tracked%eval(x,stress,tangent,status,dsdp=dsdp)
                call expect_status(status,cf_ok)
                do k = 1,nparameters
                    do i = 1,npoints
                        do j = 1,nmandel
                            call expect_close(dsdp(j,k,i),prior_dsdp(j,k,i))
                        end do
                    end do
                end do
                call plain%eval(x,reference,status=status)
                call expect_status(status,cf_ok)
                do k = 1,nparameters
                    do refinement = 1,2
                        call positive(k,refinement)%eval(x,plus_stress,status=status)
                        call expect_status(status,cf_ok)
                        call negative(k,refinement)%eval(x,minus_stress,status=status)
                        call expect_status(status,cf_ok)
                    end do
                end do
            end if
            call tracked%commit(status)
            call expect_status(status,cf_ok)
            call plain%commit(status)
            call expect_status(status,cf_ok)
            do k = 1,nparameters
                do refinement = 1,2
                    call positive(k,refinement)%commit(status)
                    call expect_status(status,cf_ok)
                    call negative(k,refinement)%commit(status)
                    call expect_status(status,cf_ok)
                end do
            end do
        end do
        call tracked%commit(status)
        call expect_status(status,cf_no_trial)
    end subroutine test_full_history_sensitivities

    subroutine test_reversed_plastic_state()
        type(material_t) :: mat
        real(wp) :: x(nmandel,1), stress(nmandel,1), dsdp(nmandel,1,1)
        real(wp) :: state(nstate,1), state_dsdp(nstate,1,1)
        real(wp) :: expected(nstate), expected_dsdp(nstate)
        integer :: status, i

        call mat%init(j2_plasticity_t(bulk,shear,initial_yield,hardening),1,status, &
            sensitivity_parameters=[character(len=12) :: "yield_stress"])
        call expect_status(status,cf_ok)
        x = 0.0_wp
        x(1,1) = 0.2_wp
        call mat%eval(x,stress,status=status,dsdp=dsdp)
        call expect_status(status,cf_ok)
        call mat%commit(status)
        call expect_status(status,cf_ok)

        ! Reverse yield: gamma1=11/130, gamma2=16/169; accumulated strain increases.
        x(1,1) = -0.1_wp
        call mat%eval(x,stress,status=status,dsdp=dsdp)
        call expect_status(status,cf_ok)
        call expect_close(stress(1,1),-7366.0_wp/507.0_wp)
        call expect_close(stress(2,1),-3922.0_wp/507.0_wp)
        call expect_close(stress(3,1),-3922.0_wp/507.0_wp)
        call expect_close(dsdp(1,1,1),-88.0_wp/169.0_wp)
        call expect_close(dsdp(2,1,1),44.0_wp/169.0_wp)
        call expect_close(dsdp(3,1,1),44.0_wp/169.0_wp)
        call mat%commit(status)
        call expect_status(status,cf_ok)
        call mat%export_committed(state,status,sensitivity_state=state_dsdp)
        call expect_status(status,cf_ok)

        ! Differentiate both accepted plastic steps with respect to initial yield stress.
        expected = 0.0_wp
        expected([1,2,3,nstate]) = [-17.0_wp/1690.0_wp,17.0_wp/3380.0_wp, &
            17.0_wp/3380.0_wp,303.0_wp/1690.0_wp]
        expected_dsdp = 0.0_wp
        expected_dsdp([1,2,3,nstate]) = [11.0_wp/1690.0_wp,-11.0_wp/3380.0_wp, &
            -11.0_wp/3380.0_wp,-37.0_wp/1690.0_wp]
        do i = 1,nstate
            call expect_close(state(i,1),expected(i))
            call expect_close(state_dsdp(i,1,1),expected_dsdp(i))
        end do
    end subroutine test_reversed_plastic_state

    subroutine test_failure_and_yield_switch()
        type(material_t) :: tracked, plain
        real(wp) :: x(nmandel,npoints), stress(nmandel,npoints)
        real(wp) :: dsdp(nmandel,1,npoints), bad(nmandel-1,npoints)
        real(wp) :: unused_derivative(nmandel,1,npoints)
        integer :: status

        call tracked%init(j2_plasticity_t(bulk,shear,initial_yield,hardening), &
            npoints,status,sensitivity_parameters=[character(len=12) :: "yield_stress"])
        call expect_status(status,cf_ok)
        call plain%init(j2_plasticity_t(bulk,shear,initial_yield,hardening),npoints,status)
        call expect_status(status,cf_ok)
        x = 0.0_wp
        x(1,:) = initial_yield/(2.0_wp*shear)
        call tracked%eval(x,stress,status=status,dsdp=dsdp)
        call expect_status(status,cf_nondifferentiable)
        call tracked%commit(status)
        call expect_status(status,cf_no_trial)
        call plain%eval(x,stress,status=status)
        call expect_status(status,cf_ok)
        call plain%revert()

        x(1,:) = x(1,:)+0.01_wp
        x(1,2) = ieee_value(0.0_wp,ieee_quiet_nan)
        call tracked%eval(x,stress,status=status,dsdp=dsdp)
        call expect_status(status,cf_evaluation_failed)
        if (tracked%last_failure_point() /= 2) error stop "J2 invalid point mismatch"
        call tracked%commit(status)
        call expect_status(status,cf_no_trial)
        x(1,2) = x(1,1)
        call tracked%eval(x,stress,status=status)
        call expect_status(status,cf_invalid_argument)
        call tracked%eval(x,stress,status=status,dsdp=dsdp)
        call expect_status(status,cf_ok)
        call tracked%revert()

        bad = 0.0_wp
        call plain%eval(bad,stress,status=status)
        call expect_status(status,cf_invalid_shape)
        call plain%eval(x,stress,status=status,dsdp=unused_derivative)
        call expect_status(status,cf_unsupported)
        call plain%init(j2_plasticity_t(bulk,shear,-1.0_wp,hardening),npoints,status)
        call expect_status(status,cf_invalid_argument)
        call plain%eval(x,stress,status=status)
        call expect_status(status,cf_ok)
        call plain%commit(status)
        call expect_status(status,cf_ok)
    end subroutine test_failure_and_yield_switch

end program test_j2_plasticity

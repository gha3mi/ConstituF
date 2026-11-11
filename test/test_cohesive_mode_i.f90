! SPDX-License-Identifier: Apache-2.0

program test_cohesive_mode_i
    !! Test convention: opening in m, traction in Pa and stiffness in Pa/m.
    !! Interface energy per area is J/m2. Explicit rescaling cases state their factors.
    !!
    !! Bilinear Mode I energy, irreversible damage, tangents and path sensitivities.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use constituf, only: wp, lk, material_t, cf_ok, cf_invalid_argument, &
                        cf_invalid_shape, cf_evaluation_failed, cf_no_trial, &
                        cf_nondifferentiable, cf_unsupported
    use constituf_cohesive_mode_i, only: cohesive_mode_i_t
    implicit none (type, external)

    integer, parameter :: nscalar = 1, npoints = 2, nparameters = 3
    real(wp), parameter :: stiffness = 2000.0_wp
    real(wp), parameter :: onset = 0.01_wp, failure = 0.05_wp

    call test_bilinear_loading()
    call test_energy_and_history()
    call test_scaled_path(1.0_wp)
    call test_scaled_path(sqrt(tiny(1.0_wp))*epsilon(1.0_wp))
    call test_scaled_path(sqrt(huge(1.0_wp))/epsilon(1.0_wp))
    call test_separated_thresholds()
    call test_threshold_history()
    call test_complete_path_derivatives()
    call test_nondifferentiable_and_failure()

contains

    subroutine expect_status(actual, expected)
        integer, intent(in) :: actual, expected

        if (actual /= expected) then
            print *, "Expected status", expected, "got", actual
            error stop "Cohesive material API status mismatch"
        end if
    end subroutine expect_status

    subroutine close(actual, expected)
        real(wp), intent(in) :: actual, expected
        real(wp) :: tolerance

        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) then
            error stop "Nonfinite cohesive comparison"
        end if
        tolerance = 512.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(expected))
        if (abs(actual-expected) > tolerance) then
            print *, "Expected", expected, "got", actual
            error stop "Cohesive traction mismatch"
        end if
    end subroutine close

    subroutine test_bilinear_loading()
        type(material_t) :: mat, single
        real(wp) :: x(nscalar,npoints), t(nscalar,npoints)
        real(wp) :: tangent(nscalar,nscalar,npoints), one(nscalar,nscalar)
        real(wp) :: fraction
        integer(lk) :: expected_bytes
        integer :: status
        integer :: empty_shape(0)

        call mat%init(cohesive_mode_i_t(stiffness,onset,failure),npoints,status)
        call expect_status(status,cf_ok)
        call single%init(cohesive_mode_i_t(stiffness,onset,failure),1,status)
        call expect_status(status,cf_ok)
        if (mat%input_size() /= nscalar .or. mat%output_size() /= nscalar) then
            error stop "Cohesive scalar field dimension"
        end if
        if (mat%state_size() /= nscalar) error stop "Missing maximum-opening history"
        if (.not. mat%field_matches("input","opening.normal","scalar", &
                                    "normal_separation","interface",empty_shape)) then
            error stop "Normal opening metadata mismatch"
        end if
        if (.not. mat%field_matches("output","traction.normal","scalar", &
                                    "normal_traction","interface",empty_shape)) then
            error stop "Normal traction metadata mismatch"
        end if
        if (.not. mat%field_matches("state","opening.maximum","scalar", &
                                    "maximum_opening","interface",empty_shape)) then
            error stop "Maximum opening metadata mismatch"
        end if
        expected_bytes = 2_lk*int(npoints,lk)*int(storage_size(0.0_wp)/8,lk)
        if (mat%state_storage_bytes() /= expected_bytes) error stop "Wrong history bytes"

        x(1,:) = [0.005_wp,0.020_wp]
        call mat%eval(x,t,tangent,status)
        call expect_status(status,cf_ok)
        call close(t(1,1),stiffness*x(1,1))
        call close(tangent(1,1,1),stiffness)
        call close(t(1,2),stiffness*onset*(failure-x(1,2))/(failure-onset))
        call close(tangent(1,1,2),-stiffness*onset/(failure-onset))
        call check_tangent(mat,x,tangent)
        call single%eval(x(:,1:1),one,status=status)
        call expect_status(status,cf_ok)
        call close(one(1,1),t(1,1))
        call mat%commit(status)
        call expect_status(status,cf_ok)

        fraction = onset*(failure-0.020_wp)/(0.020_wp*(failure-onset))
        x(1,:) = [0.003_wp,0.012_wp]
        call mat%eval(x,t,tangent,status)
        call expect_status(status,cf_ok)
        call close(t(1,1),stiffness*x(1,1))
        call close(t(1,2),stiffness*x(1,2)*fraction)
        call close(tangent(1,1,2),stiffness*fraction)
        call check_tangent(mat,x,tangent)
        call mat%revert()
        call mat%eval(x,t,tangent,status)
        call expect_status(status,cf_ok)
        call close(t(1,2),stiffness*x(1,2)*fraction)
        call mat%commit(status)
        call expect_status(status,cf_ok)

        x(1,:) = [0.007_wp,0.034_wp]
        call mat%eval(x,t,tangent,status)
        call expect_status(status,cf_ok)
        call close(t(1,2),stiffness*onset*(failure-x(1,2))/(failure-onset))
        call close(tangent(1,1,2),-stiffness*onset/(failure-onset))
        call mat%commit(status)
        call expect_status(status,cf_ok)

        x(1,:) = [0.007_wp,0.055_wp]
        call mat%eval(x,t,tangent,status)
        call expect_status(status,cf_ok)
        call close(t(1,2),0.0_wp)
        call close(tangent(1,1,2),0.0_wp)
        call mat%commit(status)
        call expect_status(status,cf_ok)
        x(1,2) = 0.014_wp
        call mat%eval(x,t,tangent,status)
        call expect_status(status,cf_ok)
        call close(t(1,2),0.0_wp)
        call close(tangent(1,1,2),0.0_wp)
    end subroutine test_bilinear_loading

    subroutine check_tangent(mat,x,tangent)
        !! Independent centered differences on smooth loading/unloading branches.
        !! Restore the unperturbed trial before the caller commits it.
        type(material_t), intent(inout) :: mat
        real(wp), intent(in) :: x(nscalar,npoints), tangent(nscalar,nscalar,npoints)
        real(wp) :: changed(nscalar,npoints)
        real(wp) :: yplus(nscalar,npoints), yminus(nscalar,npoints)
        real(wp) :: restored(nscalar,npoints)
        real(wp) :: h, numerical, error_bound
        integer :: status, ipoint

        ! Traction is affine on each smooth branch. This finite step stays
        ! inside the chosen branches and has no truncation error.
        h = 0.125_wp*minval(x)
        do ipoint = 1, npoints
            changed = x
            changed(1,ipoint) = changed(1,ipoint)+h
            call mat%eval(changed,yplus,status=status)
            call expect_status(status,cf_ok)
            changed(1,ipoint) = changed(1,ipoint)-2.0_wp*h
            call mat%eval(changed,yminus,status=status)
            call expect_status(status,cf_ok)
            numerical = (yplus(1,ipoint)-yminus(1,ipoint))/(2.0_wp*h)
            error_bound = 32.0_wp*epsilon(1.0_wp)* &
                (abs(yplus(1,ipoint))+abs(yminus(1,ipoint)))/(2.0_wp*h)
            if (.not. all(ieee_is_finite([numerical,tangent(1,1,ipoint),error_bound]))) then
                error stop "Nonfinite cohesive tangent comparison"
            end if
            if (abs(numerical-tangent(1,1,ipoint)) > error_bound) then
                error stop "Cohesive consistent opening tangent mismatch"
            end if
        end do
        call mat%eval(x,restored,status=status)
        call expect_status(status,cf_ok)
    end subroutine check_tangent

    subroutine test_energy_and_history()
        !! Integrate monotone bilinear traction exactly via its two triangles.
        type(material_t) :: mat
        real(wp) :: x(nscalar,nscalar), traction(nscalar,nscalar)
        real(wp) :: peak, energy, reference
        integer :: status, i

        call mat%init(cohesive_mode_i_t(stiffness,onset,failure),1,status)
        call expect_status(status,cf_ok)
        x(1,1) = 0.0_wp
        call mat%eval(x,traction,status=status)
        call expect_status(status,cf_ok)
        call close(traction(1,1),0.0_wp)
        call mat%commit(status)
        call expect_status(status,cf_ok)

        x(1,1) = onset
        call mat%eval(x,traction,status=status)
        call expect_status(status,cf_ok)
        peak = traction(1,1)
        call close(peak,stiffness*onset)
        call mat%commit(status)
        call expect_status(status,cf_ok)
        energy = 0.5_wp*peak*onset

        x(1,1) = failure
        call mat%eval(x,traction,status=status)
        call expect_status(status,cf_ok)
        call close(traction(1,1),0.0_wp)
        call mat%commit(status)
        call expect_status(status,cf_ok)
        energy = energy + 0.5_wp*peak*(failure-onset)
        reference = 0.5_wp*stiffness*onset*failure
        call close(energy,reference)

        ! No healing: after complete decohesion, every admissible opening
        ! has zero tensile traction. Repeated accept/reject preserves this.
        do i = 1, 3
            x(1,1) = 0.005_wp*real(i,wp)
            call mat%eval(x,traction,status=status)
            call expect_status(status,cf_ok)
            call close(traction(1,1),0.0_wp)
            call mat%revert()
            call mat%eval(x,traction,status=status)
            call expect_status(status,cf_ok)
            call close(traction(1,1),0.0_wp)
        end do
    end subroutine test_energy_and_history

    subroutine test_scaled_path(length_scale)
        !! Changing length units leaves the traction fixed. Normalized tangent
        !! and parameter derivatives expose otherwise hidden product overflow.
        real(wp), intent(in) :: length_scale
        type(material_t) :: mat
        real(wp), parameter :: openings(2) = [0.02_wp,0.012_wp]
        real(wp), parameter :: tractions(2) = [15.0_wp,9.0_wp]
        real(wp), parameter :: slopes(2) = [-500.0_wp,750.0_wp]
        real(wp), parameter :: onset_derivative(2) = [1875.0_wp,1125.0_wp]
        real(wp), parameter :: failure_derivative(2) = [125.0_wp,75.0_wp]
        real(wp) :: x(nscalar,1), y(nscalar,1), tangent(nscalar,nscalar,1)
        real(wp) :: dsdp(nscalar,nparameters,1), state(nscalar,1)
        real(wp) :: state_derivative(nscalar,nparameters,1)
        integer :: status, step

        call mat%init(cohesive_mode_i_t(stiffness/length_scale,onset*length_scale, &
            failure*length_scale),1,status,sensitivity_parameters=[character(len=16) :: &
                "failure_opening","onset_opening","stiffness"])
        call expect_status(status,cf_ok)
        do step = 1,2
            x(1,1) = openings(step)*length_scale
            call mat%eval(x,y,tangent,status,dsdp=dsdp)
            call expect_status(status,cf_ok)
            call close(y(1,1),tractions(step))
            call close(tangent(1,1,1)*length_scale,slopes(step))
            call close(dsdp(1,1,1)*length_scale,failure_derivative(step))
            call close(dsdp(1,2,1)*length_scale,onset_derivative(step))
            call close(dsdp(1,3,1)/length_scale,tractions(step)/stiffness)
            call mat%commit(status)
            call expect_status(status,cf_ok)
            call mat%export_committed(state,status,sensitivity_state=state_derivative)
            call expect_status(status,cf_ok)
            call close(state(1,1)/length_scale,0.02_wp)
            if (.not. all(ieee_is_finite(state_derivative))) error stop "Nonfinite cohesive history derivative"
            call close(maxval(abs(state_derivative)),0.0_wp)
        end do
    end subroutine test_scaled_path

    subroutine test_separated_thresholds()
        !! A remote failure threshold must not broaden the onset switch band.
        type(material_t) :: mat
        real(wp) :: x(nscalar,1), y(nscalar,1), tangent(nscalar,nscalar,1)
        real(wp) :: dsdp(nscalar,nparameters,1)
        real(wp) :: scale
        integer :: status

        scale = epsilon(1.0_wp)**2
        call mat%init(cohesive_mode_i_t(stiffness,onset,onset/scale),1,status, &
            sensitivity_parameters=[character(len=16) :: &
                "failure_opening","onset_opening","stiffness"])
        call expect_status(status,cf_ok)
        x(1,1) = 2.0_wp*onset
        call mat%eval(x,y,tangent,status,dsdp=dsdp)
        call expect_status(status,cf_ok)
        call close(y(1,1),20.0_wp)
        call close(tangent(1,1,1)/scale,-stiffness)
        call close(dsdp(1,1,1)/scale**2,stiffness)
        call close(dsdp(1,2,1),stiffness)
        call close(dsdp(1,3,1),onset)
    end subroutine test_separated_thresholds

    subroutine test_threshold_history()
        !! Restart an accepted, prescribed history with its exact zero history
        !! derivatives. Positive unloading retains threshold-parameter kinks.
        type(material_t) :: plain, tracked
        real(wp) :: x(nscalar,npoints), y(nscalar,npoints)
        real(wp) :: tangent(nscalar,nscalar,npoints), state(nscalar,npoints)
        real(wp) :: dsdp(nscalar,nparameters,npoints)
        real(wp) :: state_derivative(nscalar,nparameters,npoints)
        integer :: status, offset

        call plain%init(cohesive_mode_i_t(stiffness,onset,failure),npoints,status)
        call expect_status(status,cf_ok)
        x(1,:) = [onset,failure]
        ! First arrival at either threshold has a kink in the opening response.
        call plain%eval(x,y,tangent,status)
        call expect_status(status,cf_nondifferentiable)
        if (plain%last_failure_point() /= 1) error stop "Missing cohesive onset tangent kink"
        call plain%commit(status)
        call expect_status(status,cf_no_trial)
        x(1,1) = 0.5_wp*onset
        call plain%eval(x,y,tangent,status)
        call expect_status(status,cf_nondifferentiable)
        if (plain%last_failure_point() /= 2) error stop "Missing cohesive first-failure tangent kink"
        call plain%commit(status)
        call expect_status(status,cf_no_trial)
        x(1,:) = [onset,failure]
        call plain%eval(x,y,status=status)
        call expect_status(status,cf_ok)
        call plain%commit(status)
        call expect_status(status,cf_ok)
        call plain%export_committed(state,status)
        call expect_status(status,cf_ok)

        ! Once failure is committed, traction is identically zero as opening varies.
        ! Keep the accepted maximum exactly at failure for all three tangent probes.
        do offset = -1,1
            x(1,:) = [0.5_wp*onset,failure]
            if (offset /= 0) x(1,2) = nearest(failure,real(offset,wp))
            call plain%eval(x,y,tangent,status)
            call expect_status(status,cf_ok)
            call close(y(1,1),0.5_wp*stiffness*onset)
            call close(y(1,2),0.0_wp)
            call close(tangent(1,1,1),stiffness)
            call close(tangent(1,1,2),0.0_wp)
        end do
        call plain%revert()
        call tracked%init(cohesive_mode_i_t(stiffness,onset,failure),npoints,status, &
            sensitivity_parameters=[character(len=16) :: &
                "failure_opening","onset_opening","stiffness"])
        call expect_status(status,cf_ok)
        state_derivative = 0.0_wp
        call tracked%import_committed(state,status,sensitivity_state=state_derivative)
        call expect_status(status,cf_ok)
        x = 0.5_wp*state
        call plain%eval(x,y,tangent,status)
        call expect_status(status,cf_ok)
        call close(tangent(1,1,1),stiffness)
        call close(tangent(1,1,2),0.0_wp)
        call tracked%eval(x,y,status=status,dsdp=dsdp)
        call expect_status(status,cf_nondifferentiable)
        call tracked%commit(status)
        call expect_status(status,cf_no_trial)
        ! Isolate the failed point: its failure-opening parameter derivative
        ! remains undefined at positive opening despite the zero opening tangent.
        x(1,1) = 0.0_wp
        call tracked%eval(x,y,status=status,dsdp=dsdp)
        call expect_status(status,cf_nondifferentiable)
        if (tracked%last_failure_point() /= 2) error stop "Missing cohesive failure-parameter kink"
        call tracked%commit(status)
        call expect_status(status,cf_no_trial)

        ! At exactly zero opening, traction is zero for every parameter even
        ! when the accepted maximum lies on a damage threshold.
        x = 0.0_wp
        call tracked%eval(x,y,tangent,status,dsdp=dsdp)
        call expect_status(status,cf_ok)
        if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(dsdp))) then
            error stop "Nonfinite cohesive threshold response"
        end if
        call close(maxval(abs(y)),0.0_wp)
        call close(maxval(abs(dsdp)),0.0_wp)
        call close(tangent(1,1,1),stiffness)
        call close(tangent(1,1,2),0.0_wp)
    end subroutine test_threshold_history

    subroutine perturb(model, id, delta)
        type(cohesive_mode_i_t), intent(inout) :: model
        integer, intent(in) :: id
        real(wp), intent(in) :: delta

        select case (id)
        case (1)
            model%stiffness = model%stiffness+delta
        case (2)
            model%onset_opening = model%onset_opening+delta
        case (3)
            model%failure_opening = model%failure_opening+delta
        case default
            error stop "Unknown cohesive parameter"
        end select
    end subroutine perturb

    subroutine test_complete_path_derivatives()
        !! Derivatives compared with entirely independent accepted histories.
        type(material_t) :: mat, plain
        type(material_t) :: plus(nparameters), minus(nparameters)
        type(cohesive_mode_i_t) :: positive_model, negative_model
        real(wp) :: x(nscalar,npoints), y(nscalar,npoints)
        real(wp) :: reference(nscalar,npoints), tangent(nscalar,nscalar,npoints)
        real(wp) :: dsdp(nscalar,nparameters,npoints)
        real(wp) :: prior_dsdp(nscalar,nparameters,npoints)
        real(wp) :: yplus(nscalar,npoints), yminus(nscalar,npoints)
        real(wp) :: parameters(nparameters), perturbations(nparameters)
        real(wp) :: h, numerical, error_bound
        integer(lk) :: history_difference
        integer :: status, step, k, id, ipoint

        parameters = [stiffness,onset,failure]
        call mat%init(cohesive_mode_i_t(stiffness,onset,failure),npoints,status, &
            sensitivity_parameters=[character(len=16) :: &
                "failure_opening","onset_opening","stiffness"])
        call expect_status(status,cf_ok)
        call plain%init(cohesive_mode_i_t(stiffness,onset,failure),npoints,status)
        call expect_status(status,cf_ok)
        history_difference = 2_lk*int(npoints,lk)*int(nparameters,lk)* &
                             int(storage_size(0.0_wp)/8,lk)
        if (mat%state_storage_bytes()-plain%state_storage_bytes() /= &
            history_difference) error stop "Cohesive sensitivity allocation mismatch"
        if (plain%sensitivity_size() /= 0) error stop "Unrequested derivative history"

        do k = 1, nparameters
            id = nparameters+1-k
            if (id == 1) then
                h = 0.25_wp*parameters(id)
            else
                h = epsilon(1.0_wp)**(1.0_wp/3.0_wp)*parameters(id)
            end if
            perturbations(k) = h
            positive_model = cohesive_mode_i_t(stiffness,onset,failure)
            negative_model = cohesive_mode_i_t(stiffness,onset,failure)
            call perturb(positive_model,id,h)
            call perturb(negative_model,id,-h)
            call plus(k)%init(positive_model,npoints,status)
            call expect_status(status,cf_ok)
            call minus(k)%init(negative_model,npoints,status)
            call expect_status(status,cf_ok)
        end do

        do step = 1, 5
            select case (step)
            case (1)
                x(1,:) = [0.006_wp,0.009_wp]
            case (2)
                x(1,:) = [0.022_wp,0.018_wp]
            case (3)
                x(1,:) = [0.012_wp,0.010_wp]
            case (4)
                x(1,:) = [0.038_wp,0.032_wp]
            case (5)
                x(1,:) = [0.014_wp,0.040_wp]
            case default
                error stop "Unexpected cohesive loading step"
            end select
            call mat%eval(x,y,tangent,status,dsdp=dsdp)
            call expect_status(status,cf_ok)
            call plain%eval(x,reference,status=status)
            call expect_status(status,cf_ok)
            do ipoint = 1, npoints
                call close(y(1,ipoint),reference(1,ipoint))
            end do

            do k = 1, nparameters
                call plus(k)%eval(x,yplus,status=status)
                call expect_status(status,cf_ok)
                call minus(k)%eval(x,yminus,status=status)
                call expect_status(status,cf_ok)
                do ipoint = 1, npoints
                    numerical = (yplus(1,ipoint)-yminus(1,ipoint))/ &
                                (2.0_wp*perturbations(k))
                    error_bound = 32.0_wp*epsilon(1.0_wp)* &
                        (abs(yplus(1,ipoint))+abs(yminus(1,ipoint)))/(2.0_wp*perturbations(k))
                    if (k < nparameters) then
                        ! Each threshold derivative has a squared denominator
                        ! (failure-onset); central-difference relative error is
                        ! (h/(failure-onset))**2 to leading order.
                        error_bound = error_bound+2.0_wp*(perturbations(k)/(failure-onset))**2* &
                            abs(dsdp(1,k,ipoint))
                    end if
                    if (.not. all(ieee_is_finite([numerical,dsdp(1,k,ipoint),error_bound]))) then
                        error stop "Nonfinite cohesive parameter comparison"
                    end if
                    if (abs(numerical-dsdp(1,k,ipoint)) > error_bound) then
                        error stop "Full-history cohesive parameter derivative mismatch"
                    end if
                end do
            end do

            if (step == 3) then
                prior_dsdp = dsdp
                call mat%revert()
                call plain%revert()
                do k = 1, nparameters
                    call plus(k)%revert()
                    call minus(k)%revert()
                end do
                call mat%eval(x,y,tangent,status,dsdp=dsdp)
                call expect_status(status,cf_ok)
                do k = 1, nparameters
                    do ipoint = 1, npoints
                        call close(dsdp(1,k,ipoint),prior_dsdp(1,k,ipoint))
                    end do
                end do
                call plain%eval(x,reference,status=status)
                call expect_status(status,cf_ok)
                do k = 1, nparameters
                    call plus(k)%eval(x,yplus,status=status)
                    call expect_status(status,cf_ok)
                    call minus(k)%eval(x,yminus,status=status)
                    call expect_status(status,cf_ok)
                end do
            end if
            call mat%commit(status)
            call expect_status(status,cf_ok)
            call plain%commit(status)
            call expect_status(status,cf_ok)
            do k = 1, nparameters
                call plus(k)%commit(status)
                call expect_status(status,cf_ok)
                call minus(k)%commit(status)
                call expect_status(status,cf_ok)
            end do
        end do
    end subroutine test_complete_path_derivatives

    subroutine test_nondifferentiable_and_failure()
        type(material_t) :: tracked, plain
        real(wp) :: x(nscalar,npoints), y(nscalar,npoints)
        real(wp) :: tangent(nscalar,nscalar,npoints)
        real(wp) :: dsdp(nscalar,1,npoints), bad(nscalar-1,npoints)
        real(wp) :: extra_dsdp(nscalar,1,npoints)
        integer :: status

        call tracked%init(cohesive_mode_i_t(stiffness,onset,failure),npoints,status, &
            sensitivity_parameters=[character(len=12) :: "stiffness"])
        call expect_status(status,cf_ok)
        call plain%init(cohesive_mode_i_t(stiffness,onset,failure),npoints,status)
        call expect_status(status,cf_ok)
        x(1,:) = onset
        call tracked%eval(x,y,status=status,dsdp=dsdp)
        call expect_status(status,cf_ok)
        if (.not. all(ieee_is_finite(dsdp))) error stop "Nonfinite cohesive stiffness derivative"
        call close(maxval(abs(dsdp-onset)),0.0_wp)
        call tracked%commit(status)
        call expect_status(status,cf_ok)
        call plain%eval(x,y,status=status)
        call expect_status(status,cf_ok)
        call plain%revert()

        x(1,:) = 0.025_wp
        call tracked%eval(x,y,tangent,status,dsdp=dsdp)
        call expect_status(status,cf_ok)
        call tracked%commit(status)
        call expect_status(status,cf_ok)
        call tracked%eval(x,y,tangent,status,dsdp=dsdp)
        call expect_status(status,cf_nondifferentiable)
        call tracked%commit(status)
        call expect_status(status,cf_no_trial)
        call tracked%eval(x,y,status=status,dsdp=dsdp)
        call expect_status(status,cf_ok)
        call tracked%revert()

        x(1,:) = [0.020_wp,-0.001_wp]
        call tracked%eval(x,y,status=status,dsdp=dsdp)
        call expect_status(status,cf_evaluation_failed)
        if (tracked%last_failure_point() /= 2) error stop "Wrong failed cohesive point"
        call tracked%commit(status)
        call expect_status(status,cf_no_trial)
        x(1,2) = 0.020_wp
        call tracked%eval(x,y,status=status)
        call expect_status(status,cf_invalid_argument)
        call tracked%eval(x,y,status=status,dsdp=dsdp)
        call expect_status(status,cf_ok)
        call tracked%revert()

        ! The stiffness derivative is also defined at complete separation.
        x = failure
        call tracked%eval(x,y,status=status,dsdp=dsdp)
        call expect_status(status,cf_ok)
        if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(dsdp))) then
            error stop "Nonfinite cohesive failure response"
        end if
        call close(maxval(abs(y)),0.0_wp)
        call close(maxval(abs(dsdp)),0.0_wp)
        call tracked%revert()

        ! Complete failure must return zero without first forming K*opening.
        x = huge(1.0_wp)
        call tracked%eval(x,y,tangent,status,dsdp=dsdp)
        call expect_status(status,cf_ok)
        if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(tangent)) .or. &
            .not. all(ieee_is_finite(dsdp))) error stop "Nonfinite fully failed cohesive response"
        call close(maxval(abs(y)),0.0_wp)
        call close(maxval(abs(tangent)),0.0_wp)
        call close(maxval(abs(dsdp)),0.0_wp)
        call tracked%revert()
        x = 0.02_wp

        bad = 0.0_wp
        call plain%eval(bad,y,status=status)
        call expect_status(status,cf_invalid_shape)
        call plain%eval(x,y,status=status,dsdp=extra_dsdp)
        call expect_status(status,cf_unsupported)
        call plain%init(cohesive_mode_i_t(-1.0_wp,onset,failure),npoints,status)
        call expect_status(status,cf_invalid_argument)
        call plain%init(cohesive_mode_i_t(stiffness,0.0_wp,failure),npoints,status)
        call expect_status(status,cf_invalid_argument)
        call plain%init(cohesive_mode_i_t(stiffness,onset,onset),npoints,status)
        call expect_status(status,cf_invalid_argument)
        call plain%init(cohesive_mode_i_t(stiffness,onset, &
            ieee_value(0.0_wp,ieee_quiet_nan)),npoints,status)
        call expect_status(status,cf_invalid_argument)
        call plain%eval(x,y,status=status)
        call expect_status(status,cf_ok)
        call plain%commit(status)
        call expect_status(status,cf_ok)
    end subroutine test_nondifferentiable_and_failure

end program test_cohesive_mode_i

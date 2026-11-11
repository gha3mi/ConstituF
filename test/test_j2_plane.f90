! SPDX-License-Identifier: Apache-2.0

program test_j2_plane
    !! Test convention: stress and stress-like material parameters in Pa.
    !! Strain and deformation gradients are dimensionless. Explicit rescaling
    !! cases state their factors.
    !!
    !! Plane constraints, 3D comparisons and exact plastic loading/unloading references.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use constituf, only: wp, material_t, verify_tangent, cf_ok, cf_invalid_argument, cf_evaluation_failed, cf_no_trial
    use constituf_j2_plane, only: j2_plane_t
    use constituf_j2_plasticity, only: j2_plasticity_t
    use constituf_plane_elastic, only: plane_elastic_t
    implicit none (type, external)

    integer, parameter :: nmandel = 6, nplane = 3, nstate3 = 7, nstate2 = 8
    real(wp), parameter :: bulk = 150000.0_wp, shear = 75000.0_wp
    real(wp), parameter :: yield_stress = 350.0_wp, hardening = 2000.0_wp
    integer, parameter :: inplane(nplane) = [1,2,6]
    type(material_t) :: plane, reference, elastic, backup
    real(wp) :: x(nplane,1), y(nplane,1), c(nplane,nplane,1)
    real(wp) :: x3(nmandel,1), stress3(nmandel,1), c3(nmandel,nmandel,1)
    real(wp) :: history2(nstate2,1), history3(nstate3,1), before(nstate2,1)
    real(wp) :: y_elastic(nplane,1)
    real(wp) :: tolerance, err, h, stress_scale, difference_tolerance
    real(wp), parameter :: reference_tolerance = 512.0_wp*epsilon(1.0_wp)
    real(wp), parameter :: factors(nmandel) = [0.05_wp,0.5_wp,1.0_wp,0.7_wp,-0.4_wp,0.3_wp]
    integer :: status, mode, step, i, j

    call plane%init(j2_plane_t(bulk,shear,yield_stress,hardening,.true.),1,status, &
        sensitivity_parameters=[character(len=12) :: "bulk_modulus"])
    if (status /= cf_ok) error stop "2D J2 analytical sensitivity setup"
    call plane%init(j2_plane_t(bulk,0.0_wp,yield_stress,hardening),1,status)
    if (status /= cf_invalid_argument) error stop "2D J2 accepted invalid shear modulus"
    call plane%init(j2_plane_t(bulk,shear,yield_stress,-1.0_wp),1,status)
    if (status /= cf_invalid_argument) error stop "2D J2 accepted negative hardening"

    do mode = 1, 2
        call plane%init(j2_plane_t(bulk,shear,yield_stress,hardening,mode==2),1,status)
        if (status /= cf_ok) error stop "2D J2 model setup"
        call reference%init(j2_plasticity_t(bulk,shear,yield_stress,hardening),1,status)
        if (status /= cf_ok) error stop "3D J2 model setup"
        if (.not. plane%field_matches("input","strain.small","mandel","small_strain","spatial",[2,2])) then
            error stop "Incorrect in-plane Mandel field declaration"
        end if
        if (plane%state_size() /= nstate2) error stop "2D J2 history count"

        do step = 1, size(factors)
            x(:,1) = factors(step)*[0.008_wp,-0.003_wp,0.0015_wp]
            call plane%eval(x,y,c,status)
            if (status /= cf_ok) error stop "2D J2 step failed"
            if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(c))) then
                error stop "Nonfinite plane J2 response"
            end if
            if (step == 1) then
                call elastic%init(plane_elastic_t(bulk,shear,mode==2),1,status)
                if (status /= cf_ok) error stop "Elastic 2D reference setup"
                call elastic%eval(x,y_elastic,status=status)
                if (status /= cf_ok) error stop "Elastic 2D reference response"
                if (.not. all(ieee_is_finite(y_elastic))) error stop "Nonfinite elastic reference"
                if (maxval(abs(y-y_elastic)) > reference_tolerance*max(1.0_wp,maxval(abs(y_elastic)))) then
                    error stop "2D J2 elastic limit"
                end if
            end if

            ! Store the accepted out-of-plane strain. The 3D reference
            ! retains its own independent committed material-point history.
            call plane%commit(status)
            if (status /= cf_ok) error stop "2D J2 accept step"
            call plane%export_committed(history2,status)
            if (status /= cf_ok) error stop "2D J2 state export"
            if (.not. all(ieee_is_finite(history2))) error stop "Nonfinite plane J2 history"
            x3 = 0.0_wp
            x3(inplane,1) = x(:,1)
            x3(3,1) = history2(nstate2,1)
            call reference%eval(x3,stress3,c3,status)
            if (status /= cf_ok) error stop "3D J2 independent reference failed"
            if (.not. all(ieee_is_finite(stress3)) .or. .not. all(ieee_is_finite(c3))) then
                error stop "Nonfinite 3D J2 reference"
            end if
            ! Include the 128*eps local plane-stress residual scale.
            tolerance = reference_tolerance*max(1.0_wp,yield_stress,maxval(abs(stress3)))
            if (maxval(abs(y(:,1)-stress3(inplane,1))) > tolerance) then
                error stop "2D J2 must match equilibrated 3D stress"
            end if
            if (mode == 2) then
                if (abs(stress3(3,1)) > tolerance) error stop "Plane stress sigma33 not zero"
            else
                if (abs(x3(3,1)) > epsilon(1.0_wp)) error stop "Plane strain epsilon33 not zero"
            end if
            do j = 1, nplane
                do i = 1, nplane
                    err = c3(inplane(i),inplane(j),1)
                    if (mode == 2) err = err - c3(inplane(i),3,1)*c3(3,inplane(j),1)/c3(3,3,1)
                    if (.not. ieee_is_finite(err)) error stop "Nonfinite plane J2 Schur reference"
                    if (abs(c(i,j,1)-err) > reference_tolerance*max(1.0_wp,abs(err))) then
                        error stop "Inconsistent plane J2 Schur tangent"
                    end if
                end do
            end do
            call reference%commit(status)
            if (status /= cf_ok) error stop "3D J2 reference acceptance"
            call reference%export_committed(history3,status)
            if (status /= cf_ok) error stop "3D J2 state export"
            if (.not. all(ieee_is_finite(history3))) error stop "Nonfinite 3D J2 history"
            if (maxval(abs(history2(1:nstate3,:)-history3)) > reference_tolerance) then
                error stop "2D/3D J2 committed plastic history differs"
            end if
        end do

        ! A verification trial is automatically reverted and must not change
        ! the accepted plastic history even after many perturbation evaluations.
        x(:,1) = [0.010_wp,-0.003_wp,0.002_wp]
        call plane%eval(x,y,c,status)
        if (status /= cf_ok) error stop "2D J2 tangent verification reference"
        if (.not. all(ieee_is_finite(y))) error stop "Nonfinite tangent verification reference"
        h = maxval(abs(x))*epsilon(1.0_wp)**(1.0_wp/3.0_wp)
        stress_scale = max(1.0_wp,yield_stress,maxval(abs(y)))
        difference_tolerance = 32.0_wp*epsilon(1.0_wp)*stress_scale/h
        ! Account separately for roundoff and the local Newton residual in each perturbed trial.
        if (mode == 2) difference_tolerance = difference_tolerance + &
            256.0_wp*epsilon(1.0_wp)*stress_scale/h
        call verify_tangent(plane,x,err,status,step=[(h,i=1,nplane)],absolute_tolerance=difference_tolerance)
        if (status /= cf_ok .or. .not. ieee_is_finite(err) .or. err > 1.0_wp) then
            error stop "2D J2 finite-difference tangent verification"
        end if
        call plane%export_committed(before,status)
        if (status /= cf_ok) error stop "2D J2 checkpoint after verification"
        if (.not. all(abs(before-history2) <= 0.0_wp)) error stop "Verification modified committed state"
        x(1,1) = ieee_value(0.0_wp,ieee_quiet_nan)
        call plane%eval(x,y,status=status)
        if (status /= cf_evaluation_failed) error stop "Invalid 2D J2 input not rejected"
        call plane%commit(status)
        if (status /= cf_no_trial) error stop "Failed 2D J2 trial was committed"
        call plane%export_committed(before,status)
        if (status /= cf_ok) error stop "2D J2 checkpoint after failed trial"
        if (.not. all(abs(before-history2) <= 0.0_wp)) error stop "Failed trial changed history"
    end do

    call backup%init(j2_plane_t(bulk,shear,yield_stress,hardening,.true.),1,status)
    if (status /= cf_ok) error stop "2D J2 reinitialized model"
    call backup%import_committed(history2,status)
    if (status /= cf_ok) error stop "2D J2 in-memory restart"
    call test_exact_plane_stress_history(1.0_wp)
    ! Changing stress units must preserve the constraint solve and smooth plastic sensitivities.
    call test_exact_plane_stress_history(epsilon(1.0_wp)**2)

contains

    subroutine test_exact_plane_stress_history(modulus_scale)
        real(wp), intent(in) :: modulus_scale
        type(material_t) :: material
        real(wp) :: driving(nplane,1), stress(nplane,1), jacobian(nplane,nplane,1)
        real(wp) :: dsdp(nplane,1,1), state(nstate2,1), state_dsdp(nstate2,1,1)
        real(wp) :: expected_stress(nplane), expected_tangent(nplane,nplane)
        real(wp) :: expected_state(nstate2), expected_state_dsdp(nstate2)
        real(wp), parameter :: expected_dsdp(nplane) = [16.0_wp/19.0_wp,-8.0_wp/57.0_wp,0.0_wp]
        integer :: result, load_step, row, col

        ! E=150, nu=1/4: the first plastic state has uniaxial stress 36/11.
        call material%init(j2_plane_t(100.0_wp*modulus_scale,60.0_wp*modulus_scale, &
            3.0_wp*modulus_scale,15.0_wp*modulus_scale,.true.),1,result, &
            sensitivity_parameters=[character(len=12) :: "yield_stress"])
        if (result /= cf_ok) error stop "Exact plane-stress history setup"
        driving(:,1) = [1.0_wp/25.0_wp,-4.0_wp/275.0_wp,0.0_wp]
        expected_stress = [36.0_wp/11.0_wp,0.0_wp,0.0_wp]
        expected_tangent = reshape([660.0_wp,840.0_wp,0.0_wp,840.0_wp,1760.0_wp, &
            0.0_wp,0.0_wp,0.0_wp,1140.0_wp],[nplane,nplane])/19.0_wp
        expected_state = 0.0_wp
        expected_state([1,2,3,7,8]) = [1.0_wp/55.0_wp,-1.0_wp/110.0_wp, &
            -1.0_wp/110.0_wp,1.0_wp/55.0_wp,-4.0_wp/275.0_wp]
        expected_state_dsdp = 0.0_wp
        expected_state_dsdp([1,2,3,7,8]) = [-1.0_wp/171.0_wp,2.0_wp/855.0_wp, &
            1.0_wp/285.0_wp,-1.0_wp/171.0_wp,2.0_wp/855.0_wp]

        do load_step = 1,2
            if (load_step == 2) then
                ! Elastic unloading retains plastic history and its transverse-corrected derivative.
                driving = 0.5_wp*driving
                expected_stress = [4.0_wp/11.0_wp,4.0_wp/11.0_wp,0.0_wp]
                expected_state(8) = -17.0_wp/1650.0_wp
                expected_tangent = reshape([160.0_wp,40.0_wp,0.0_wp,40.0_wp,160.0_wp, &
                    0.0_wp,0.0_wp,0.0_wp,120.0_wp],[nplane,nplane])
            end if
            call material%eval(driving,stress,jacobian,result,dsdp=dsdp)
            if (result /= cf_ok) error stop "Exact plane-stress history evaluation"
            do row = 1,nplane
                call check(stress(row,1)/modulus_scale,expected_stress(row))
                call check(dsdp(row,1,1),expected_dsdp(row))
                do col = 1,nplane
                    call check(jacobian(row,col,1)/modulus_scale,expected_tangent(row,col))
                end do
            end do
            call material%commit(result)
            if (result /= cf_ok) error stop "Exact plane-stress history acceptance"
            call material%export_committed(state,result,sensitivity_state=state_dsdp)
            if (result /= cf_ok) error stop "Exact plane-stress history export"
            do row = 1,nstate2
                call check(state(row,1),expected_state(row))
                call check(state_dsdp(row,1,1)*modulus_scale,expected_state_dsdp(row))
            end do
        end do
    end subroutine test_exact_plane_stress_history

    subroutine check(actual, expected)
        real(wp), intent(in) :: actual, expected

        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) then
            error stop "Nonfinite exact plane-stress comparison"
        end if
        if (abs(actual-expected) > reference_tolerance*max(1.0_wp,abs(expected))) then
            print *, "Expected",expected,"got",actual
            error stop "Exact plane-stress reference mismatch"
        end if
    end subroutine check

end program test_j2_plane

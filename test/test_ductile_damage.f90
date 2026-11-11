! SPDX-License-Identifier: Apache-2.0

program test_ductile_damage
    !! Test convention: moduli, stress and energy density in MPa (= N/mm2).
    !! Strain and damage are dimensionless.
    !!
    !! Ductile damage verification: J2 limit, exact discrete tangent,
    !! positive local dissipation, degradation, checkpoint and rollback.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
    use constituf, only: wp, material_t, verify_potential, verify_tangent, cf_ok, cf_invalid_argument, &
                        cf_evaluation_failed, cf_no_trial, cf_nondifferentiable
    use constituf_ductile_damage, only: ductile_damage_t
    use constituf_j2_plasticity, only: j2_plasticity_t
    implicit none (type, external)

    integer, parameter :: nmandel = 6, nstate = 8, noutput = 8, npoints = 2
    real(wp), parameter :: bulk = 150000.0_wp, shear = 75000.0_wp
    real(wp), parameter :: yield_stress = 350.0_wp, hardening = 2000.0_wp
    real(wp), parameter :: alpha0 = 0.001_wp, rate = 150.0_wp, damage_cap = 0.85_wp
    real(wp), parameter :: tolerance = 64.0_wp*epsilon(1.0_wp)
    type(material_t) :: material, j2, undamaged, restored, invalid
    real(wp) :: strain(nmandel,npoints), elastic_input(nmandel,npoints)
    real(wp) :: response(noutput,npoints), base_stress(nmandel,npoints), undamaged_y(noutput,npoints)
    real(wp) :: tangent(noutput,nmandel,npoints), j2_tangent(nmandel,nmandel,npoints)
    real(wp) :: history(nstate,npoints), previous_history(nstate,npoints), backup(nstate,npoints)
    real(wp) :: state_j2(nstate-1,npoints)
    real(wp) :: expected_damage, expected_dissipation, error, h
    real(wp) :: trace_strain, elastic_energy, hardening_energy
    real(wp) :: elastic_strain(nmandel), elastic_dev(nmandel)
    integer :: point_codes(npoints)
    integer :: status, step, i
    logical :: saw_damage

    call material%init(ductile_damage_t(bulk,shear,yield_stress,hardening,alpha0,rate,damage_cap),npoints,status, &
        sensitivity_parameters=[character(len=12) :: "bulk_modulus"])
    if (status /= cf_ok) error stop "Damage analytical sensitivity setup"
    call material%init(ductile_damage_t(bulk,shear,yield_stress,hardening,alpha0,rate,damage_cap),npoints,status)
    if (status /= cf_ok) error stop "Damage material setup"
    if (material%input_size() /= nmandel .or. material%output_size() /= noutput .or. &
        material%state_size() /= nstate) error stop "Damage declared field sizes"
    if (material%output_offset("energy.helmholtz") /= nmandel+1 .or. &
        material%output_offset("energy.dissipation_increment") /= noutput) then
        error stop "Damage output fields not declared"
    end if

    call j2%init(j2_plasticity_t(bulk,shear,yield_stress,hardening),npoints,status)
    if (status /= cf_ok) error stop "J2 effective stress reference"
    call undamaged%init(ductile_damage_t(bulk,shear,yield_stress,hardening,alpha0,rate,0.0_wp),npoints,status)
    if (status /= cf_ok) error stop "Zero-damage J2 specialization"

    strain = 0.0_wp
    strain(1,:) = [0.0007_wp,-0.0005_wp]
    strain(2,:) = [-0.00035_wp,0.00025_wp]
    strain(3,:) = [-0.00035_wp,0.00025_wp]
    call material%eval(strain,response,tangent,status)
    if (status /= cf_ok) error stop "Initial elastic damage trial"
    if (.not. all(ieee_is_finite(response))) error stop "Nonfinite initial damage response"
    if (maxval(abs(response(noutput,:))) > 0.0_wp) error stop "Damage dissipation before yielding"
    call material%commit(status)
    if (status /= cf_ok) error stop "Initial damage commit"
    call j2%eval(strain,base_stress,status=status)
    if (status /= cf_ok) error stop "Initial J2 response"
    call j2%commit(status)
    if (status /= cf_ok) error stop "Initial J2 commit"
    call undamaged%eval(strain,undamaged_y,status=status)
    if (status /= cf_ok) error stop "Initial undamaged response"
    call undamaged%commit(status)
    if (status /= cf_ok) error stop "Initial undamaged commit"
    call material%export_committed(previous_history,status)
    if (status /= cf_ok .or. .not. all(ieee_is_finite(previous_history)) .or. &
        maxval(abs(previous_history)) > tolerance) then
        error stop "Virgin damage state must remain zero"
    end if

    saw_damage = .false.
    loading_loop: do step = 1,5
        strain(:,1) = real(step,wp)*[0.004_wp,-0.002_wp,-0.002_wp,0.0008_wp,0.0_wp,0.0002_wp]
        strain(:,2) = real(step,wp)*[-0.003_wp,0.0015_wp,0.0015_wp,0.0004_wp,0.0_wp,-0.0001_wp]
        call material%eval(strain,response,tangent,status)
        if (status /= cf_ok) error stop "Damage active loading"
        call j2%eval(strain,base_stress,j2_tangent,status)
        if (status /= cf_ok) error stop "J2 reference loading"
        call undamaged%eval(strain,undamaged_y,status=status)
        if (status /= cf_ok) error stop "Undamaged loading"
        if (.not. all(ieee_is_finite(response)) .or. .not. all(ieee_is_finite(base_stress)) .or. &
            .not. all(ieee_is_finite(undamaged_y))) error stop "Nonfinite damage loading comparison"
        if (maxval(abs(undamaged_y(1:nmandel,:)-base_stress)) > tolerance*maxval(abs(base_stress))) then
            error stop "Zero damage cap does not reproduce J2 stress"
        end if
        call material%commit(status)
        if (status /= cf_ok) error stop "Damage accepted loading"
        call j2%commit(status)
        if (status /= cf_ok) error stop "J2 accepted loading"
        call undamaged%commit(status)
        if (status /= cf_ok) error stop "Undamaged accepted loading"
        call material%export_committed(history,status)
        if (status /= cf_ok) error stop "Damage state export"
        call j2%export_committed(state_j2,status)
        if (status /= cf_ok .or. .not. all(ieee_is_finite(history)) .or. .not. all(ieee_is_finite(state_j2)) .or. &
            maxval(abs(history(1:nstate-1,:)-state_j2)) > tolerance) then
            error stop "Damage changed the effective J2 plastic history"
        end if
        do i = 1,npoints
            expected_damage = 0.0_wp
            if (history(nstate-1,i) > alpha0) then
                expected_damage = damage_cap*(1.0_wp-exp(-rate*(history(nstate-1,i)-alpha0)))
            end if
            if (abs(history(nstate,i)-expected_damage) > tolerance) then
                error stop "Ductile damage law does not match committed alpha"
            end if
            if (history(nstate,i) < previous_history(nstate,i)) then
                error stop "Damage decreased after an accepted plastic step"
            end if
            if (history(nstate,i) > 0.0_wp) saw_damage = .true.
            if (maxval(abs(response(1:nmandel,i)-(1.0_wp-expected_damage)*base_stress(:,i))) > &
                tolerance*maxval(abs(base_stress(:,i)))) then
                error stop "Effective-stress damage degradation mismatch"
            end if

            elastic_strain = strain(:,i)-history(1:nmandel,i)
            trace_strain = sum(elastic_strain(1:3))
            elastic_dev = elastic_strain
            elastic_dev(1:3) = elastic_dev(1:3)-trace_strain/3.0_wp
            elastic_energy = 0.5_wp*bulk*trace_strain**2 + shear*dot_product(elastic_dev,elastic_dev)
            hardening_energy = 0.5_wp*hardening*history(nstate-1,i)**2
            if (abs(response(nmandel+1,i)-(1.0_wp-expected_damage)*(elastic_energy+hardening_energy)) > &
                tolerance*(elastic_energy+hardening_energy)) then
                error stop "Stored Helmholtz energy does not match damage potential"
            end if
            expected_dissipation = (1.0_wp-expected_damage)*yield_stress* &
                (history(nstate-1,i)-previous_history(nstate-1,i)) + &
                (elastic_energy+hardening_energy)*(expected_damage-previous_history(nstate,i))
            if (abs(response(noutput,i)-expected_dissipation) > tolerance*max(abs(response(noutput,i)), &
                abs(expected_dissipation))) then
                error stop "Discrete plastic/damage dissipation mismatch"
            end if
            if (response(noutput,i) < 0.0_wp) error stop "Negative ductile damage dissipation"
        end do
        previous_history = history
    end do loading_loop
    if (.not. saw_damage) error stop "Ductile damage never activated"
    if (maxval(history(nstate,:)) >= damage_cap) error stop "Damage incorrectly saturated"

    ! Verify *all eight* exact discrete output derivatives away from
    ! activation and yielding; committed damage and alpha remain fixed.
    strain(:,1) = [0.028_wp,-0.014_wp,-0.014_wp,0.0052_wp,0.0_wp,0.0014_wp]
    strain(:,2) = [-0.020_wp,0.010_wp,0.010_wp,0.0030_wp,0.0_wp,-0.0007_wp]
    call material%eval(strain,response,status=status)
    if (status /= cf_ok) error stop "Damage tangent reference response"
    h = maxval(abs(strain))*epsilon(1.0_wp)**(1.0_wp/3.0_wp)
    call verify_tangent(material,strain,error,status, &
        absolute_tolerance=64.0_wp*epsilon(1.0_wp)*maxval(abs(response))/h,step=[(h,i=1,nmandel)])
    if (status /= cf_ok .or. .not. ieee_is_finite(error) .or. error > 1.0_wp) then
        error stop "Damage algorithmic stress/energy/dissipation tangent"
    end if
    call material%export_committed(backup,status)
    if (status /= cf_ok .or. .not. all(ieee_is_finite(backup)) .or. maxval(abs(backup-history)) > 0.0_wp) then
        error stop "Numerical tangent verification altered committed damage"
    end if

    ! After unloading, no new plastic strain or damage occurs. At fixed
    ! committed history the Helmholtz gradient equals current stress.
    elastic_input(:,1) = 0.99_wp*[0.020_wp,-0.010_wp,-0.010_wp,0.004_wp,0.0_wp,0.001_wp]
    elastic_input(:,2) = 0.99_wp*[-0.015_wp,0.0075_wp,0.0075_wp,0.002_wp,0.0_wp,-0.0005_wp]
    call material%eval(elastic_input,response,status=status)
    if (status /= cf_ok) error stop "Damage unloaded response"
    if (.not. all(ieee_is_finite(response))) error stop "Nonfinite unloaded damage response"
    if (maxval(abs(response(noutput,:))) > tolerance) error stop "Elastic unloading produced plastic dissipation"
    call material%revert()
    h = maxval(abs(elastic_input))*epsilon(1.0_wp)**(1.0_wp/3.0_wp)
    call verify_potential(material,elastic_input,nmandel+1,[(i,i=1,nmandel)], &
        [(1.0_wp,i=1,nmandel)],error,status, &
        absolute_tolerance=64.0_wp*epsilon(1.0_wp)*maxval(abs(response(nmandel+1,:)))/h, &
        step=[(h,i=1,nmandel)])
    if (status /= cf_ok .or. .not. ieee_is_finite(error) .or. error > 1.0_wp) then
        error stop "Elastic damage branch violates Helmholtz stress conjugacy"
    end if

    call restored%init(ductile_damage_t(bulk,shear,yield_stress,hardening,alpha0,rate,damage_cap),npoints,status)
    if (status /= cf_ok) error stop "Damage restart model initialization"
    call restored%import_committed(history,status)
    if (status /= cf_ok) error stop "Damage committed state restart"
    call material%eval(strain,response,status=status)
    if (status /= cf_ok) error stop "Damage original continuation"
    call restored%eval(strain,undamaged_y,status=status)
    if (status /= cf_ok .or. .not. all(ieee_is_finite(response)) .or. &
        .not. all(ieee_is_finite(undamaged_y)) .or. maxval(abs(response-undamaged_y)) > tolerance) then
        error stop "Damage restarted response does not match original"
    end if
    call material%revert()
    call restored%revert()

    strain(1,2) = ieee_value(0.0_wp,ieee_quiet_nan)
    call material%eval(strain,response,status=status)
    if (status /= cf_evaluation_failed) error stop "Damage accepted nonfinite input"
    call material%commit(status)
    if (status /= cf_no_trial) error stop "Failed damage trial was committed"
    call material%point_statuses(point_codes,status)
    if (status /= cf_ok) error stop "Point diagnostics not retained for damage failure"
    if (point_codes(2) == cf_ok .or. point_codes(1) /= cf_ok) then
        error stop "Ductile damage point status is incorrect"
    end if
    call material%export_committed(backup,status)
    if (status /= cf_ok .or. .not. all(ieee_is_finite(backup)) .or. maxval(abs(backup-history)) > 0.0_wp) then
        error stop "Rejected damage trial modified accepted state"
    end if

    call invalid%init(ductile_damage_t(bulk,shear,yield_stress,hardening,alpha0,rate,1.0_wp),npoints,status)
    if (status /= cf_invalid_argument) error stop "Damage accepted unphysical Dmax=1"
    call invalid%init(ductile_damage_t(bulk,shear,yield_stress,hardening,alpha0,-rate,damage_cap),npoints,status)
    if (status /= cf_invalid_argument) error stop "Damage accepted negative rate"
    call invalid%init(ductile_damage_t(bulk,shear,yield_stress,hardening,-alpha0,rate,damage_cap),npoints,status)
    if (status /= cf_invalid_argument) error stop "Damage accepted negative initiation strain"
    call invalid%init(ductile_damage_t(bulk,shear,yield_stress,-hardening,alpha0,rate,damage_cap),npoints,status)
    if (status /= cf_invalid_argument) error stop "Damage accepted negative hardening"

    call test_exact_damage(1.0_wp)
    call test_exact_damage(epsilon(1.0_wp)**2)
    call test_small_damage_scales()
    call test_branch_switches()
    call test_path_sensitivities()

contains

    subroutine test_exact_damage(modulus_scale)
        real(wp), intent(in) :: modulus_scale
        type(material_t) :: model
        real(wp), parameter :: direction(nmandel) = [1.0_wp,-0.5_wp,-0.5_wp,0.0_wp,0.0_wp,0.0_wp]
        real(wp) :: x(nmandel,1), y(noutput,1), jacobian(noutput,nmandel,1), derivative(noutput,2,1)
        real(wp) :: state(nstate,1), state_derivative(nstate,2,1)
        real(wp) :: expected(noutput), expected_derivative(noutput,2)
        real(wp) :: logarithm
        integer :: result, load_step, component

        ! alpha=1/50, b*(alpha-alpha0)=log(2), D=2/5 and available energy=133/4000.
        logarithm = log(2.0_wp)
        call model%init(ductile_damage_t(100.0_wp*modulus_scale,60.0_wp*modulus_scale, &
            3.0_wp*modulus_scale,15.0_wp*modulus_scale,0.01_wp,100.0_wp*logarithm,0.8_wp),1,result, &
            sensitivity_parameters=[character(len=12) :: "yield_stress", "damage_rate"])
        call expect(result)
        x(:,1) = (23.0_wp/600.0_wp)*direction
        do load_step = 1,2
            expected(1:nmandel) = (33.0_wp/25.0_wp)*direction
            expected(nmandel+1:noutput) = [399.0_wp/20000.0_wp,493.0_wp/10000.0_wp]
            expected_derivative(1:nmandel,1) = ((72.0_wp+88.0_wp*logarithm)/195.0_wp)*direction
            expected_derivative(1:nmandel,2) = -(11.0_wp/1250.0_wp)*direction
            expected_derivative(nmandel+1,1) = 3.0_wp/325.0_wp+133.0_wp*logarithm/19500.0_wp
            expected_derivative(nmandel+1,2) = -133.0_wp/1000000.0_wp
            expected_derivative(noutput,1) = 29.0_wp/3250.0_wp+107.0_wp*logarithm/19500.0_wp
            expected_derivative(noutput,2) = -107.0_wp/1000000.0_wp
            if (load_step==2) then
                ! Elastic unloading preserves damage and all accepted plastic/damage sensitivities.
                x = 0.5_wp*x
                expected(1:nmandel) = -(3.0_wp/50.0_wp)*direction
                expected(nmandel+1:noutput) = [147.0_wp/80000.0_wp,0.0_wp]
                expected_derivative(1:nmandel,1) = ((72.0_wp-4.0_wp*logarithm)/195.0_wp)*direction
                expected_derivative(1:nmandel,2) = direction/2500.0_wp
                expected_derivative(nmandel+1,1) = -9.0_wp/6500.0_wp+49.0_wp*logarithm/78000.0_wp
                expected_derivative(nmandel+1,2) = -49.0_wp/4000000.0_wp
                expected_derivative(noutput,:) = 0.0_wp
            end if
            call model%eval(x,y,jacobian,result,dsdp=derivative)
            call expect(result)
            do component = 1,noutput
                call close(y(component,1)/modulus_scale,expected(component))
                call close(derivative(component,1,1),expected_derivative(component,1))
                call close(derivative(component,2,1)/modulus_scale,expected_derivative(component,2))
            end do
            if (load_step==1) then
                call close(jacobian(1,1,1)/modulus_scale,(828.0_wp-704.0_wp*logarithm)/13.0_wp)
                call close(jacobian(1,2,1)/modulus_scale,(756.0_wp+352.0_wp*logarithm)/13.0_wp)
                call close(jacobian(4,4,1)/modulus_scale,792.0_wp/23.0_wp)
                call close(jacobian(nmandel+1,1,1)/modulus_scale,(69.0_wp-266.0_wp*logarithm)/325.0_wp)
                call close(jacobian(noutput,1,1)/modulus_scale,(406.0_wp-214.0_wp*logarithm)/325.0_wp)
            else
                call close(jacobian(1,1,1)/modulus_scale,108.0_wp)
                call close(jacobian(1,2,1)/modulus_scale,36.0_wp)
                call close(jacobian(4,4,1)/modulus_scale,72.0_wp)
                call close(jacobian(nmandel+1,1,1)/modulus_scale,-3.0_wp/50.0_wp)
                call close(jacobian(noutput,1,1)/modulus_scale,0.0_wp)
            end if
            call model%commit(result)
            call expect(result)
            call model%export_committed(state,result,sensitivity_state=state_derivative)
            call expect(result)
            do component = 1,nmandel
                call close(state(component,1),direction(component)/50.0_wp)
                call close(state_derivative(component,1,1)*modulus_scale,-direction(component)/195.0_wp)
                call close(state_derivative(component,2,1),0.0_wp)
            end do
            call close(state(nstate-1,1),1.0_wp/50.0_wp)
            call close(state(nstate,1),2.0_wp/5.0_wp)
            call close(state_derivative(nstate-1,1,1)*modulus_scale,-1.0_wp/195.0_wp)
            call close(state_derivative(nstate,1,1)*modulus_scale,-8.0_wp*logarithm/39.0_wp)
            call close(state_derivative(nstate-1,2,1),0.0_wp)
            call close(state_derivative(nstate,2,1),1.0_wp/250.0_wp)
        end do
    end subroutine test_exact_damage

    subroutine test_small_damage_scales()
        type(material_t) :: model
        real(wp), parameter :: small_scale = epsilon(1.0_wp)**2
        real(wp), parameter :: direction(nmandel) = [1.0_wp,-0.5_wp,-0.5_wp,0.0_wp,0.0_wp,0.0_wp]
        real(wp) :: x(nmandel,1), y(noutput,1), derivative(noutput,1,1)
        real(wp) :: state(nstate,1), state_derivative(nstate,1,1)
        integer :: result, component

        ! For tiny b*(alpha-alpha0), D/[Dmax*b*(alpha-alpha0)] tends to one.
        call model%init(ductile_damage_t(100.0_wp,60.0_wp,3.0_wp,15.0_wp,0.01_wp,small_scale,0.8_wp),1,result, &
            sensitivity_parameters=[character(len=14) :: "maximum_damage"])
        call expect(result)
        x(:,1) = (23.0_wp/600.0_wp)*direction
        call model%eval(x,y,status=result,dsdp=derivative)
        call expect(result)
        call model%commit(result)
        call expect(result)
        call model%export_committed(state,result,sensitivity_state=state_derivative)
        call expect(result)
        call close(state(nstate,1)/(0.008_wp*small_scale),1.0_wp)
        call close(state_derivative(nstate,1,1)/(0.01_wp*small_scale),1.0_wp)
        call close(derivative(1,1,1)/(0.01_wp*small_scale),-11.0_wp/5.0_wp)

        ! A tiny positive alpha-alpha0 remains a smooth active branch after rescaling b.
        call model%init(ductile_damage_t(100.0_wp,60.0_wp,3.0_wp*small_scale,15.0_wp, &
            0.01_wp*small_scale,100.0_wp*log(2.0_wp)/small_scale,0.8_wp),1,result, &
            sensitivity_parameters=[character(len=14) :: "maximum_damage"])
        call expect(result)
        x(:,1) = small_scale*(23.0_wp/600.0_wp)*direction
        call model%eval(x,y,status=result,dsdp=derivative)
        call expect(result)
        do component = 1,nmandel
            call close(y(component,1)/small_scale,(33.0_wp/25.0_wp)*direction(component))
            call close(derivative(component,1,1)/small_scale,-(11.0_wp/10.0_wp)*direction(component))
        end do
        call close(y(nmandel+1,1)/small_scale**2,399.0_wp/20000.0_wp)
        call close(y(noutput,1)/small_scale**2,493.0_wp/10000.0_wp)
        call model%commit(result)
        call expect(result)
        call model%export_committed(state,result,sensitivity_state=state_derivative)
        call expect(result)
        call close(state(nstate,1),2.0_wp/5.0_wp)
        call close(state_derivative(nstate,1,1),0.5_wp)
    end subroutine test_small_damage_scales


    subroutine test_branch_switches()
        type(material_t) :: ordinary, tracked
        real(wp), parameter :: direction(nmandel) = [1.0_wp,-0.5_wp,-0.5_wp,0.0_wp,0.0_wp,0.0_wp]
        real(wp) :: x(nmandel,1), y(noutput,1), jacobian(noutput,nmandel,1), derivative(noutput,1,1)
        integer :: result, branch, component

        ! Binary-exact inputs distinguish initial yielding from damage initiation.
        call ordinary%init(ductile_damage_t(128.0_wp,64.0_wp,3.0_wp,0.0_wp, &
            1.0_wp/64.0_wp,32.0_wp,0.5_wp),1,result)
        call expect(result)
        call tracked%init(ductile_damage_t(128.0_wp,64.0_wp,3.0_wp,0.0_wp, &
            1.0_wp/64.0_wp,32.0_wp,0.5_wp),1,result, &
            sensitivity_parameters=[character(len=12) :: "yield_stress"])
        call expect(result)
        do branch = 1,2
            x(:,1) = real(branch,wp)*direction/64.0_wp
            call ordinary%eval(x,y,jacobian,result)
            call expect(result)
            do component = 1,nmandel
                call close(y(component,1),2.0_wp*direction(component))
            end do
            call close(y(nmandel+1,1),3.0_wp/128.0_wp)
            ! Ordinary tangents use the documented inactive-side convention.
            if (branch==1) then
                call close(y(noutput,1),0.0_wp)
                call close(jacobian(1,1,1),640.0_wp/3.0_wp)
                call close(jacobian(1,2,1),256.0_wp/3.0_wp)
                call close(jacobian(4,4,1),128.0_wp)
                call close(jacobian(nmandel+1,1,1),2.0_wp)
                call close(jacobian(noutput,1,1),0.0_wp)
            else
                call close(y(noutput,1),3.0_wp/64.0_wp)
                call close(jacobian(1,1,1),128.0_wp)
                call close(jacobian(1,2,1),128.0_wp)
                call close(jacobian(4,4,1),64.0_wp)
                call close(jacobian(nmandel+1,1,1),0.0_wp)
                call close(jacobian(noutput,1,1),2.0_wp)
            end if
            call tracked%eval(x,y,status=result,dsdp=derivative)
            if (result /= cf_nondifferentiable) error stop "Damage sensitivity accepted a branch switch"
        end do
    end subroutine test_branch_switches

    subroutine test_path_sensitivities()
        integer, parameter :: nparameters = 7, nsteps = 4
        character(len=18), parameter :: names(nparameters) = [character(len=18) :: "bulk_modulus", &
            "shear_modulus","yield_stress","hardening_modulus","damage_threshold","damage_rate","maximum_damage"]
        real(wp), parameter :: parameters(nparameters) = &
            [100.0_wp,60.0_wp,3.0_wp,15.0_wp,0.01_wp,100.0_wp*log(2.0_wp),0.8_wp]
        real(wp), parameter :: direction(nmandel) = [1.0_wp,-0.5_wp,-0.5_wp,0.0_wp,0.0_wp,0.0_wp]
        type(material_t) :: tracked
        type(material_t) :: plus(nparameters,2), minus(nparameters,2)
        real(wp) :: path(nmandel,nsteps), x(nmandel,1), y(noutput,1), derivative(noutput,nparameters,1)
        real(wp) :: state(nstate,1), state_derivative(nstate,nparameters,1)
        real(wp) :: pplus(nparameters), pminus(nparameters), steps(nparameters,2), elastic_dev(nmandel)
        real(wp) :: yp(noutput,1), ym(noutput,1), hp(nstate,1), hm(nstate,1)
        real(wp) :: output_difference(noutput,2), state_difference(nstate,2)
        real(wp) :: output_tolerance(noutput), state_tolerance(nstate)
        real(wp) :: output_scale(noutput), state_scale(nstate)
        real(wp) :: stress_scale, strain_scale, trace_e, available_energy, old_alpha
        integer :: result, load_step, parameter, refinement

        ! Plasticity below the threshold, active damage, elastic unloading, then reversal.
        ! Volumetric strain makes the bulk contribution to stress and energy nonzero.
        path(:,1) = (53.0_wp/2400.0_wp)*direction
        path(1:3,1) = path(1:3,1)+0.003_wp
        path(:,2) = (23.0_wp/600.0_wp)*direction
        path(1:3,2) = path(1:3,2)+0.004_wp
        path(:,3) = 0.75_wp*path(:,2)
        path(:,4) = [-0.025_wp,0.010_wp,0.006_wp,0.012_wp,-0.006_wp,0.009_wp]
        call tracked%init(ductile_damage_t(parameters(1),parameters(2),parameters(3),parameters(4), &
            parameters(5),parameters(6),parameters(7)),1,result,sensitivity_parameters=names)
        call expect(result)
        do parameter = 1,nparameters
            steps(parameter,1) = (32.0_wp*epsilon(1.0_wp))**(1.0_wp/3.0_wp)*parameters(parameter)
            ! Cap dependence is affine along a fixed effective-plastic history.
            if (parameter==nparameters) steps(parameter,1) = parameters(parameter)/8.0_wp
            steps(parameter,2) = 0.5_wp*steps(parameter,1)
            do refinement = 1,2
                pplus = parameters
                pminus = parameters
                pplus(parameter) = pplus(parameter)+steps(parameter,refinement)
                pminus(parameter) = pminus(parameter)-steps(parameter,refinement)
                call plus(parameter,refinement)%init(ductile_damage_t(pplus(1),pplus(2),pplus(3),pplus(4), &
                    pplus(5),pplus(6),pplus(7)),1,result)
                call expect(result)
                call minus(parameter,refinement)%init(ductile_damage_t(pminus(1),pminus(2),pminus(3),pminus(4), &
                    pminus(5),pminus(6),pminus(7)),1,result)
                call expect(result)
            end do
        end do
        state = 0.0_wp
        output_scale = 0.0_wp
        state_scale = 0.0_wp
        do load_step = 1,nsteps
            x(:,1) = path(:,load_step)
            old_alpha = state(nstate-1,1)
            elastic_dev = x(:,1)-state(1:nmandel,1)
            trace_e = sum(elastic_dev(1:3))
            elastic_dev(1:3) = elastic_dev(1:3)-trace_e/3.0_wp
            stress_scale = max(abs(parameters(1)*trace_e),parameters(3)+parameters(4)*old_alpha, &
                2.0_wp*parameters(2)*sqrt(1.5_wp*dot_product(elastic_dev,elastic_dev)))
            call tracked%eval(x,y,status=result,dsdp=derivative)
            call expect(result)
            call tracked%commit(result)
            call expect(result)
            call tracked%export_committed(state,result,sensitivity_state=state_derivative)
            call expect(result)
            if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(state)) .or. &
                .not. all(ieee_is_finite(derivative)) .or. .not. all(ieee_is_finite(state_derivative))) then
                error stop "Nonfinite ductile-damage history comparison"
            end if
            elastic_dev = x(:,1)-state(1:nmandel,1)
            trace_e = sum(elastic_dev(1:3))
            elastic_dev(1:3) = elastic_dev(1:3)-trace_e/3.0_wp
            available_energy = 0.5_wp*parameters(1)*trace_e**2 + parameters(2)*dot_product(elastic_dev,elastic_dev) &
                + 0.5_wp*parameters(4)*state(nstate-1,1)**2
            ! Keep energy errors separate from stress; include accepted-history roundoff.
            output_scale(1:nmandel) = output_scale(1:nmandel)+stress_scale
            output_scale(nmandel+1:noutput) = output_scale(nmandel+1:noutput)+available_energy &
                + parameters(3)*(state(nstate-1,1)-old_alpha)
            strain_scale = max(maxval(abs(x)),maxval(abs(state(1:nmandel,1))),state(nstate-1,1))
            state_scale(1:nstate-1) = state_scale(1:nstate-1)+strain_scale
            state_scale(nstate) = state_scale(nstate)+state(nstate,1)+parameters(7)*parameters(6)*strain_scale
            do parameter = 1,nparameters
                do refinement = 1,2
                    call plus(parameter,refinement)%eval(x,yp,status=result)
                    call expect(result)
                    call minus(parameter,refinement)%eval(x,ym,status=result)
                    call expect(result)
                    call plus(parameter,refinement)%commit(result)
                    call expect(result)
                    call minus(parameter,refinement)%commit(result)
                    call expect(result)
                    call plus(parameter,refinement)%export_committed(hp,result)
                    call expect(result)
                    call minus(parameter,refinement)%export_committed(hm,result)
                    call expect(result)
                    output_difference(:,refinement) = (yp(:,1)-ym(:,1))/(2.0_wp*steps(parameter,refinement))
                    state_difference(:,refinement) = (hp(:,1)-hm(:,1))/(2.0_wp*steps(parameter,refinement))
                end do
                if (.not. all(ieee_is_finite(output_difference)) .or. .not. all(ieee_is_finite(state_difference))) then
                    error stop "Nonfinite ductile-damage history secants"
                end if
                output_tolerance = abs(output_difference(:,2)-output_difference(:,1))/3.0_wp &
                    + 32.0_wp*epsilon(1.0_wp)*output_scale/steps(parameter,2)
                state_tolerance = abs(state_difference(:,2)-state_difference(:,1))/3.0_wp &
                    + 32.0_wp*epsilon(1.0_wp)*state_scale/steps(parameter,2)
                if (any(abs(output_difference(:,2)-derivative(:,parameter,1)) > output_tolerance)) then
                    print *, "Ductile-damage output history derivative",load_step,parameter
                    error stop "Ductile-damage parameter derivative mismatch"
                end if
                if (any(abs(state_difference(:,2)-state_derivative(:,parameter,1)) > state_tolerance)) then
                    print *, "Ductile-damage internal history derivative",load_step,parameter
                    error stop "Ductile-damage state derivative mismatch"
                end if
            end do
        end do
    end subroutine test_path_sensitivities

    subroutine expect(result)
        integer, intent(in) :: result

        if (result /= cf_ok) then
            print *, "Ductile damage status",result
            error stop "Exact ductile damage evaluation failed"
        end if
    end subroutine expect

    subroutine close(actual, expected)
        real(wp), intent(in) :: actual, expected

        if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) then
            error stop "Nonfinite exact ductile damage comparison"
        end if
        if (abs(actual-expected) > 512.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(expected))) then
            print *, "Expected",expected,"got",actual
            error stop "Exact ductile damage reference mismatch"
        end if
    end subroutine close
end program test_ductile_damage

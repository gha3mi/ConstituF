! SPDX-License-Identifier: Apache-2.0

program chaboche
    !! Chaboche cyclic plasticity with backstress.
    !! SI example: moduli, yield stress and backstress [Pa]. Strains and recovery [-].
    use constituf, only: wp, material_t, cf_ok
    use constituf_chaboche, only: chaboche_t
    implicit none (type, external)

    integer, parameter :: ninput = 6, noutput = 6
    type(material_t) :: material
    real(wp) :: strain(ninput,1), stress(noutput,1)
    real(wp) :: tangent(noutput,ninput,1), dsdp(noutput,1,1)
    real(wp) :: prescribed(ninput)
    real(wp), allocatable :: state(:,:), state_dsdp(:,:,:)
    integer :: status, step

    ! Select one exact analytical parameter derivative; other names can be added.
    call material%init(chaboche_t(100.0e9_wp,40.0e9_wp,180.0e6_wp,1.0e9_wp,9.0e9_wp,12.0_wp),1,status, &
        sensitivity_parameters=[character(len=17) :: "kinematic_modulus"])
    if (status /= cf_ok) error stop "Initialization failed"

    prescribed = [0.015_wp,-0.004_wp,-0.004_wp,0.0_wp,0.0_wp,0.0_wp]
    allocate(state(material%state_size(),1), state_dsdp(material%state_size(),1,1))
    do step = 1,4
        strain(:,1) = real(step,wp)/real(4,wp)*prescribed
        if (step == 4) strain(:,1) = -0.35_wp*prescribed  ! Reverse.
        call material%eval(strain,stress,tangent,status,dsdp=dsdp)
        if (status /= cf_ok) error stop "Evaluation failed"
        print "(a,i0,2(a,es12.4))", "Step ",step," stress(1) [Pa]=",stress(1,1), &
            " d(stress(1))/d(kinematic_modulus) [-]=",dsdp(1,1,1)
        call material%commit(status)
        if (status /= cf_ok) error stop "Commit failed"
    end do

    ! Committed history: plastic strain [-], backstress [Pa] and accumulated plastic strain [-].
    ! Their kinematic-modulus derivatives have units [1/Pa], [-] and [1/Pa], respectively.
    call material%export_committed(state,status,sensitivity_state=state_dsdp)
    if (status /= cf_ok) error stop "State export failed"
    print "(a,*(1x,es12.4))", "Committed plastic strain [-]:", state(1:6,1)
    print "(a,*(1x,es12.4))", "Committed backstress [Pa]:", state(7:12,1)
    print "(a,*(1x,es12.4))", "Accumulated plastic strain [-]:", state(13,1)
    print "(a,*(1x,es12.4))", "d(plastic strain)/d(kinematic_modulus) [1/Pa]:", state_dsdp(1:6,1,1)
    print "(a,*(1x,es12.4))", "d(backstress)/d(kinematic_modulus) [-]:", state_dsdp(7:12,1,1)
    print "(a,*(1x,es12.4))", "d(accumulated strain)/d(kinematic_modulus) [1/Pa]:", state_dsdp(13,1,1)

    print "(a,*(1x,es12.4))", "Cauchy stress [Pa]:", stress(:,1)
    print "(a,*(1x,es12.4))", "Tangent: first output row [Pa]:", tangent(1,:,1)

    ! Rejected trial: leave the committed state unchanged.
    strain(:,1) = 0.9_wp*prescribed
    call material%eval(strain,stress,tangent,status,dsdp=dsdp)
    if (status /= cf_ok) error stop "Trial evaluation failed"
    call material%revert()
end program chaboche

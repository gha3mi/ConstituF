! SPDX-License-Identifier: Apache-2.0

program axisymmetric_j2
    !! Axisymmetric J2 yielding and unloading.
    !! SI example: elastic/hardening moduli and yield stress [Pa]. Strains [-].
    use constituf, only: wp, material_t, cf_ok
    use constituf_axisymmetric_j2, only: axisymmetric_j2_t
    implicit none (type, external)

    integer, parameter :: ninput = 4, noutput = 4
    type(material_t) :: material
    real(wp) :: strain(ninput,1), stress(noutput,1)
    real(wp) :: tangent(noutput,ninput,1), dsdp(noutput,1,1)
    real(wp) :: prescribed(ninput)
    real(wp), allocatable :: state(:,:), state_dsdp(:,:,:)
    integer :: status, step

    ! Select one exact analytical parameter derivative; other names can be added.
    call material%init(axisymmetric_j2_t(150.0e9_wp,75.0e9_wp,350.0e6_wp,2.0e9_wp),1,status, &
        sensitivity_parameters=[character(len=12) :: "yield_stress"])
    if (status /= cf_ok) error stop "Initialization failed"

    prescribed = [0.004_wp,-0.0015_wp,-0.0025_wp,0.00075_wp]
    allocate(state(material%state_size(),1), state_dsdp(material%state_size(),1,1))
    do step = 1,5
        strain(:,1) = real(step,wp)/real(5,wp)*prescribed
        if (step == 5) strain(:,1) = 0.35_wp*prescribed  ! Unload.
        call material%eval(strain,stress,tangent,status,dsdp=dsdp)
        if (status /= cf_ok) error stop "Evaluation failed"
        print "(a,i0,2(a,es12.4))", "Step ",step," stress(1) [Pa]=",stress(1,1), &
            " d(stress(1))/d(yield_stress) [-]=",dsdp(1,1,1)
        call material%commit(status)
        if (status /= cf_ok) error stop "Commit failed"
    end do

    ! Committed history: plastic strain and accumulated plastic strain;
    ! state_dsdp holds the corresponding exact path derivatives.
    call material%export_committed(state,status,sensitivity_state=state_dsdp)
    if (status /= cf_ok) error stop "State export failed"
    print "(a,*(1x,es12.4))", "Committed internal state [-]:", state(:,1)
    print "(a,*(1x,es12.4))", "d(state)/d(yield_stress) [1/Pa]:", state_dsdp(:,1,1)

    print "(a,*(1x,es12.4))", "Cauchy stress [Pa]:", stress(:,1)
    print "(a,*(1x,es12.4))", "Tangent: first output row [Pa]:", tangent(1,:,1)

    ! Rejected trial: leave the committed state unchanged.
    strain(:,1) = 0.9_wp*prescribed
    call material%eval(strain,stress,tangent,status,dsdp=dsdp)
    if (status /= cf_ok) error stop "Trial evaluation failed"
    call material%revert()
end program axisymmetric_j2

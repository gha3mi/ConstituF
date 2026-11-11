! SPDX-License-Identifier: Apache-2.0

program hill_plasticity
    !! Orthotropic Hill plasticity.
    !! SI example: moduli and yield stress [Pa]. Hill coefficients and strains [-].
    use constituf, only: wp, material_t, cf_ok
    use constituf_hill_plasticity, only: hill_plasticity_t
    implicit none (type, external)

    integer, parameter :: ninput = 6, noutput = 6
    type(material_t) :: material
    real(wp) :: strain(ninput,1), stress(noutput,1)
    real(wp) :: tangent(noutput,ninput,1), dsdp(noutput,1,1)
    real(wp) :: prescribed(ninput)
    real(wp), allocatable :: state(:,:), state_dsdp(:,:,:)
    integer :: status, step

    ! Select one exact analytical parameter derivative; other names can be added.
    call material%init(hill_plasticity_t(150.0e9_wp,75.0e9_wp,350.0e6_wp,2.0e9_wp,[0.6_wp,0.4_wp,0.5_wp, &
        1.6_wp,1.3_wp,1.8_wp]),1,status, &
        sensitivity_parameters=[character(len=6) :: "hill_f"])
    if (status /= cf_ok) error stop "Initialization failed"

    prescribed = [0.008_wp,-0.004_wp,-0.004_wp,0.002_wp,0.0_wp,0.0_wp]
    allocate(state(material%state_size(),1), state_dsdp(material%state_size(),1,1))
    do step = 1,4
        strain(:,1) = real(step,wp)/real(4,wp)*prescribed
        if (step == 4) strain(:,1) = 0.35_wp*prescribed  ! Unload.
        call material%eval(strain,stress,tangent,status,dsdp=dsdp)
        if (status /= cf_ok) error stop "Evaluation failed"
        print "(a,i0,2(a,es12.4))", "Step ",step," stress(1) [Pa]=",stress(1,1), &
            " d(stress(1))/d(hill_f) [Pa]=",dsdp(1,1,1)
        call material%commit(status)
        if (status /= cf_ok) error stop "Commit failed"
    end do

    ! Committed history: plastic strain and accumulated plastic strain;
    ! state_dsdp holds the corresponding exact path derivatives.
    call material%export_committed(state,status,sensitivity_state=state_dsdp)
    if (status /= cf_ok) error stop "State export failed"
    print "(a,*(1x,es12.4))", "Committed internal state [-]:", state(:,1)
    print "(a,*(1x,es12.4))", "d(state)/d(hill_f) [-]:", state_dsdp(:,1,1)

    print "(a,*(1x,es12.4))", "Cauchy stress [Pa]:", stress(:,1)
    print "(a,*(1x,es12.4))", "Tangent: first output row [Pa]:", tangent(1,:,1)

    ! Rejected trial: leave the committed state unchanged.
    strain(:,1) = 0.9_wp*prescribed
    call material%eval(strain,stress,tangent,status,dsdp=dsdp)
    if (status /= cf_ok) error stop "Trial evaluation failed"
    call material%revert()
end program hill_plasticity

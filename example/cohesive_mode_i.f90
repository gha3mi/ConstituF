! SPDX-License-Identifier: Apache-2.0

program cohesive_mode_i
    !! Mode-I cohesive softening and unloading.
    !! SI example: stiffness [Pa/m], openings [m] and traction [Pa].
    use constituf, only: wp, material_t, cf_ok
    use constituf_cohesive_mode_i, only: cohesive_mode_i_t
    implicit none (type, external)

    integer, parameter :: ninput = 1, noutput = 1
    type(material_t) :: material
    real(wp) :: opening(ninput,1), traction(noutput,1)
    real(wp) :: tangent(noutput,ninput,1), dsdp(noutput,1,1)
    real(wp) :: prescribed(ninput)
    real(wp), allocatable :: state(:,:), state_dsdp(:,:,:)
    integer :: status, step

    ! Select one exact analytical parameter derivative; other names can be added.
    call material%init(cohesive_mode_i_t(2000.0_wp,0.01_wp,0.05_wp),1,status, &
        sensitivity_parameters=[character(len=13) :: "onset_opening"])
    if (status /= cf_ok) error stop "Initialization failed"

    prescribed = [0.025_wp]
    allocate(state(material%state_size(),1), state_dsdp(material%state_size(),1,1))
    do step = 1,4
        opening(:,1) = real(step,wp)/real(4,wp)*prescribed
        if (step == 4) opening(:,1) = 0.35_wp*prescribed  ! Unload.
        call material%eval(opening,traction,tangent,status,dsdp=dsdp)
        if (status /= cf_ok) error stop "Evaluation failed"
        print "(a,i0,2(a,es12.4))", "Step ",step," traction(1) [Pa]=",traction(1,1), &
            " d(traction(1))/d(onset_opening) [Pa/m]=",dsdp(1,1,1)
        call material%commit(status)
        if (status /= cf_ok) error stop "Commit failed"
    end do

    ! Committed history: maximum normal opening;
    ! state_dsdp holds the corresponding exact path derivatives.
    call material%export_committed(state,status,sensitivity_state=state_dsdp)
    if (status /= cf_ok) error stop "State export failed"
    print "(a,*(1x,es12.4))", "Committed internal state [m]:", state(:,1)
    print "(a,*(1x,es12.4))", "d(state)/d(onset_opening) [-]:", state_dsdp(:,1,1)

    print "(a,*(1x,es12.4))", "Normal interface traction [Pa]:", traction(:,1)
    print "(a,*(1x,es12.4))", "Tangent: first output row [Pa/m]:", tangent(1,:,1)

    ! Rejected trial: leave the committed state unchanged.
    opening(:,1) = 0.9_wp*prescribed
    call material%eval(opening,traction,tangent,status,dsdp=dsdp)
    if (status /= cf_ok) error stop "Trial evaluation failed"
    call material%revert()
end program cohesive_mode_i

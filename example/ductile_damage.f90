! SPDX-License-Identifier: Apache-2.0

program ductile_damage
    !! Ductile plasticity with damage and dissipation.
    !! SI example: moduli/stress [Pa] and energy density [J/m3].
    !! Strains, damage and the damage parameters are dimensionless.
    use constituf, only: wp, material_t, cf_ok
    use constituf_ductile_damage, only: ductile_damage_t
    implicit none (type, external)

    integer, parameter :: ninput = 6, noutput = 8
    type(material_t) :: material
    real(wp) :: strain(ninput,1), response(noutput,1)
    real(wp) :: tangent(noutput,ninput,1), dsdp(noutput,1,1)
    real(wp) :: prescribed(ninput)
    real(wp), allocatable :: state(:,:), state_dsdp(:,:,:)
    integer :: status, step

    ! Select one exact analytical parameter derivative; other names can be added.
    call material%init(ductile_damage_t(150.0e9_wp,75.0e9_wp,350.0e6_wp,2.0e9_wp,0.001_wp,150.0_wp, &
        0.85_wp),1,status, &
        sensitivity_parameters=[character(len=11) :: "damage_rate"])
    if (status /= cf_ok) error stop "Initialization failed"

    prescribed = [0.005_wp,-0.0025_wp,-0.0025_wp,0.001_wp,0.0_wp,0.0_wp]
    allocate(state(material%state_size(),1), state_dsdp(material%state_size(),1,1))
    do step = 1,5
        strain(:,1) = real(step,wp)/real(5,wp)*prescribed
        call material%eval(strain,response,tangent,status,dsdp=dsdp)
        if (status /= cf_ok) error stop "Evaluation failed"
        print "(a,i0,2(a,es12.4))", "Step ",step," response(1) [Pa]=",response(1,1), &
            " d(response(1))/d(damage_rate) [Pa]=",dsdp(1,1,1)
        call material%commit(status)
        if (status /= cf_ok) error stop "Commit failed"
    end do

    ! Committed history: plastic strain, accumulated plastic strain and damage;
    ! state_dsdp holds the corresponding exact path derivatives.
    call material%export_committed(state,status,sensitivity_state=state_dsdp)
    if (status /= cf_ok) error stop "State export failed"
    print "(a,*(1x,es12.4))", "Committed internal state [-]:", state(:,1)
    print "(a,*(1x,es12.4))", "d(state)/d(damage_rate) [-]:", state_dsdp(:,1,1)

    print "(a,*(1x,es12.4))", "Cauchy stress [Pa]:", response(1:6,1)
    print "(a,*(1x,es12.4))", "Helmholtz energy [J/m3]:", response(7,1)
    print "(a,*(1x,es12.4))", "Dissipation increment [J/m3]:", response(8,1)
    print "(a,*(1x,es12.4))", "Tangent: first output row [Pa]:", tangent(1,:,1)

    ! Rejected trial: leave the committed state unchanged.
    strain(:,1) = 0.9_wp*prescribed
    call material%eval(strain,response,tangent,status,dsdp=dsdp)
    if (status /= cf_ok) error stop "Trial evaluation failed"
    call material%revert()
end program ductile_damage

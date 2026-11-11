! SPDX-License-Identifier: Apache-2.0

program kelvin_voigt
    !! Kelvin-Voigt ramp and fixed-strain hold.
    !! SI example: moduli/stress [Pa], viscosities [Pa*s] and dt [s]. Strains [-].
    use constituf, only: wp, material_t, cf_ok
    use constituf_kelvin_voigt, only: kelvin_voigt_t
    implicit none (type, external)

    integer, parameter :: ninput = 6, noutput = 6
    type(material_t) :: material
    real(wp) :: strain(ninput,1), stress(noutput,1)
    real(wp) :: tangent(noutput,ninput,1), dsdp(noutput,1,1)
    real(wp) :: prescribed(ninput)
    real(wp) :: previous(ninput,1), initial(ninput,1)
    real(wp) :: state(0,1), state_dsdp(0,1,1)
    integer :: status, step

    initial = 0.0_wp
    ! Select one exact analytical parameter derivative; other names can be added.
    call material%init(kelvin_voigt_t(100.0e6_wp,40.0e6_wp,15.0e6_wp,8.0e6_wp),1,status,initial_input=initial, &
        sensitivity_parameters=[character(len=15) :: "shear_viscosity"])
    if (status /= cf_ok) error stop "Initialization failed"

    prescribed = [0.01_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp]
    do step = 1,3
        strain(:,1) = prescribed
        if (step == 1) strain(:,1) = 0.5_wp*prescribed  ! Then hold.
        call material%eval(strain,stress,tangent,status,dt=0.5_wp,dsdp=dsdp)
        if (status /= cf_ok) error stop "Evaluation failed"
        print "(a,i0,2(a,es12.4))", "Step ",step," stress(1) [Pa]=",stress(1,1), &
            " d(stress(1))/d(shear_viscosity) [1/s]=",dsdp(1,1,1)
        call material%commit(status)
        if (status /= cf_ok) error stop "Commit failed"
    end do

    ! Kelvin-Voigt retains the preceding input, not internal state variables.
    call material%export_committed(state,status,previous_input=previous,sensitivity_state=state_dsdp)
    if (status /= cf_ok) error stop "Previous-input export failed"
    print "(a,*(1x,es12.4))", "Last accepted strain [-]:", previous(:,1)

    print "(a,*(1x,es12.4))", "Cauchy stress [Pa]:", stress(:,1)
    print "(a,*(1x,es12.4))", "Tangent: first output row [Pa]:", tangent(1,:,1)

    ! Rejected trial: leave the committed state unchanged.
    strain(:,1) = 0.9_wp*prescribed
    call material%eval(strain,stress,tangent,status,dt=0.5_wp,dsdp=dsdp)
    if (status /= cf_ok) error stop "Trial evaluation failed"
    call material%revert()
end program kelvin_voigt

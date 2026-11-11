! SPDX-License-Identifier: Apache-2.0

program darcy_exponential
    !! Pressure-dependent Darcy flow.
    !! SI example: reference permeability [m2], viscosity [Pa*s], pressure coefficient [1/Pa] and reference pressure [Pa].
    !! Inputs are pressure gradient [Pa/m] and absolute pressure [Pa]. Velocity [m/s].
    use constituf, only: wp, material_t, cf_ok
    use constituf_darcy_exponential, only: darcy_exponential_t
    implicit none (type, external)

    integer, parameter :: ninput = 4, noutput = 3
    type(material_t) :: material
    real(wp) :: driving(ninput,1), velocity(noutput,1)
    real(wp) :: tangent(noutput,ninput,1), dsdp(noutput,1,1)
    integer :: status

    ! Select one exact analytical parameter derivative; other names can be added.
    call material%init(darcy_exponential_t(1.0e-12_wp,1.0e-3_wp,1.0e-6_wp,1.0e6_wp),1,status, &
        sensitivity_parameters=[character(len=20) :: "pressure_coefficient"])
    if (status /= cf_ok) error stop "Initialization failed"

    ! Input: pressure gradient [Pa/m] and absolute pressure [Pa].
    driving(:,1) = [1.0e6_wp,0.0_wp,-0.5e6_wp,1.2e6_wp]
    call material%eval(driving,velocity,tangent,status,dsdp=dsdp)
    if (status /= cf_ok) error stop "Evaluation failed"

    print "(a,*(1x,es12.4))", "Darcy velocity [m/s]:", velocity(:,1)
    print "(a,*(1x,es12.4))", "Tangent: row 1, columns 1:3 [m2/(Pa*s)]:", tangent(1,1:3,1)
    print "(a,*(1x,es12.4))", "Tangent: row 1, remaining columns [m/(Pa*s)]:", tangent(1,4:,1)
    print "(a,*(1x,es12.4))", "d(response)/d(pressure_coefficient) [Pa*m/s]:", dsdp(:,1,1)
end program darcy_exponential

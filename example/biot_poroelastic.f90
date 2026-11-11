! SPDX-License-Identifier: Apache-2.0

program biot_poroelastic
    !! Biot strain-pressure coupling.
    !! SI example: elastic/Biot moduli and pore pressure [Pa], Biot coefficient [-].
    !! Strains and fluid content are dimensionless. Total stress [Pa].
    use constituf, only: wp, material_t, cf_ok
    use constituf_biot_poroelastic, only: biot_poroelastic_t
    implicit none (type, external)

    integer, parameter :: ninput = 7, noutput = 7
    type(material_t) :: material
    real(wp) :: driving(ninput,1), response(noutput,1)
    real(wp) :: tangent(noutput,ninput,1), dsdp(noutput,1,1)
    integer :: status

    ! Select one exact analytical parameter derivative; other names can be added.
    call material%init(biot_poroelastic_t(5.0e9_wp,2.0e9_wp,0.85_wp,1.0e10_wp),1,status, &
        sensitivity_parameters=[character(len=16) :: "biot_coefficient"])
    if (status /= cf_ok) error stop "Initialization failed"

    ! Input: six Mandel strains followed by pore pressure [Pa].
    driving(:,1) = [0.001_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,3.0e6_wp]
    call material%eval(driving,response,tangent,status,dsdp=dsdp)
    if (status /= cf_ok) error stop "Evaluation failed"

    print "(a,*(1x,es12.4))", "Total Cauchy stress [Pa]:", response(1:6,1)
    print "(a,*(1x,es12.4))", "Fluid content [-]:", response(7,1)
    print "(a,*(1x,es12.4))", "Tangent: row 1, columns 1:6 [Pa]:", tangent(1,1:6,1)
    print "(a,*(1x,es12.4))", "Tangent: row 1, remaining columns [-]:", tangent(1,7:,1)
    print "(a,*(1x,es12.4))", "d(stress)/d(biot_coefficient) [Pa]:", dsdp(1:6,1,1)
    print "(a,*(1x,es12.4))", "d(fluid content)/d(biot_coefficient) [-]:", dsdp(7,1,1)
end program biot_poroelastic

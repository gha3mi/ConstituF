! SPDX-License-Identifier: Apache-2.0

program isotropic_thermoelastic
    !! Coupled strain-temperature elasticity.
    !! SI example: moduli/stress [Pa], expansion coefficient [1/K] and absolute temperatures [K].
    !! Strains are dimensionless.
    use constituf, only: wp, material_t, cf_ok
    use constituf_isotropic_thermoelastic, only: isotropic_thermoelastic_t
    implicit none (type, external)

    integer, parameter :: ninput = 7, noutput = 6
    type(material_t) :: material
    real(wp) :: driving(ninput,1), stress(noutput,1)
    real(wp) :: tangent(noutput,ninput,1), dsdp(noutput,1,1)
    integer :: status

    ! Select one exact analytical parameter derivative; other names can be added.
    call material%init(isotropic_thermoelastic_t(100.0e9_wp,40.0e9_wp,12.0e-6_wp,293.15_wp),1,status, &
        sensitivity_parameters=[character(len=17) :: "thermal_expansion"])
    if (status /= cf_ok) error stop "Initialization failed"

    driving(:,1) = [0.001_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,303.15_wp]
    call material%eval(driving,stress,tangent,status,dsdp=dsdp)
    if (status /= cf_ok) error stop "Evaluation failed"

    print "(a,*(1x,es12.4))", "Cauchy stress [Pa]:", stress(:,1)
    print "(a,*(1x,es12.4))", "Tangent: row 1, columns 1:6 [Pa]:", tangent(1,1:6,1)
    print "(a,*(1x,es12.4))", "Tangent: row 1, remaining columns [Pa/K]:", tangent(1,7:,1)
    print "(a,*(1x,es12.4))", "d(response)/d(thermal_expansion) [Pa*K]:", dsdp(:,1,1)
end program isotropic_thermoelastic

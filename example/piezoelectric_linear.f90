! SPDX-License-Identifier: Apache-2.0

program piezoelectric_linear
    !! Piezoelectric stress-charge coupling.
    !! SI example: moduli [Pa], piezoelectric coefficients [C/m2] and permittivities [F/m].
    !! Strains [-], electric field [V/m], stress [Pa], displacement [C/m2] and enthalpy [J/m3].
    use constituf, only: wp, material_t, cf_ok
    use constituf_piezoelectric_linear, only: piezoelectric_linear_t
    implicit none (type, external)

    integer, parameter :: ninput = 9, noutput = 10
    type(material_t) :: material
    real(wp) :: driving(ninput,1), response(noutput,1)
    real(wp) :: tangent(noutput,ninput,1), dsdp(noutput,1,1)
    integer :: status

    ! Select one exact analytical parameter derivative; other names can be added.
    call material%init(piezoelectric_linear_t(100.0e9_wp,40.0e9_wp,-2.0_wp,4.0_wp,3.0_wp,1.0e-8_wp, &
        1.5e-8_wp),1,status, &
        sensitivity_parameters=[character(len=3) :: "e33"])
    if (status /= cf_ok) error stop "Initialization failed"

    ! Input: six Mandel strains followed by electric field [V/m].
    driving(:,1) = [0.001_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0e3_wp,0.0_wp,2.0e3_wp]
    call material%eval(driving,response,tangent,status,dsdp=dsdp)
    if (status /= cf_ok) error stop "Evaluation failed"

    print "(a,*(1x,es12.4))", "Cauchy stress [Pa]:", response(1:6,1)
    print "(a,*(1x,es12.4))", "Electric displacement [C/m2]:", response(7:9,1)
    print "(a,*(1x,es12.4))", "Electric enthalpy [J/m3]:", response(10,1)
    print "(a,*(1x,es12.4))", "Tangent: row 1, columns 1:6 [Pa]:", tangent(1,1:6,1)
    print "(a,*(1x,es12.4))", "Tangent: row 1, remaining columns [C/m2]:", tangent(1,7:,1)
    print "(a,*(1x,es12.4))", "d(stress)/d(e33) [V/m]:", dsdp(1:6,1,1)
    print "(a,*(1x,es12.4))", "d(electric displacement)/d(e33) [-]:", dsdp(7:9,1,1)
    print "(a,*(1x,es12.4))", "d(enthalpy)/d(e33) [V/m]:", dsdp(10,1,1)
end program piezoelectric_linear

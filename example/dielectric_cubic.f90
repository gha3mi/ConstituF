! SPDX-License-Identifier: Apache-2.0

program dielectric_cubic
    !! Nonlinear cubic dielectric displacement.
    !! SI example: permittivity [F/m] and cubic coefficient [F*m/V2].
    !! Electric field [V/m] and electric displacement [C/m2].
    use constituf, only: wp, material_t, cf_ok
    use constituf_dielectric_cubic, only: dielectric_cubic_t
    implicit none (type, external)

    integer, parameter :: ninput = 3, noutput = 3
    type(material_t) :: material
    real(wp) :: electric_field(ninput,1), displacement(noutput,1)
    real(wp) :: tangent(noutput,ninput,1), dsdp(noutput,1,1)
    integer :: status

    ! Select one exact analytical parameter derivative; other names can be added.
    call material%init(dielectric_cubic_t(2.0e-11_wp,1.0e-22_wp),1,status, &
        sensitivity_parameters=[character(len=17) :: "cubic_coefficient"])
    if (status /= cf_ok) error stop "Initialization failed"

    electric_field(:,1) = [1.0e5_wp,0.0_wp,-0.5e5_wp]
    call material%eval(electric_field,displacement,tangent,status,dsdp=dsdp)
    if (status /= cf_ok) error stop "Evaluation failed"

    print "(a,*(1x,es12.4))", "Electric displacement [C/m2]:", displacement(:,1)
    print "(a,*(1x,es12.4))", "Tangent: first output row [F/m]:", tangent(1,:,1)
    print "(a,*(1x,es12.4))", "d(response)/d(cubic_coefficient) [V3/m3]:", dsdp(:,1,1)
end program dielectric_cubic

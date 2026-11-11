! SPDX-License-Identifier: Apache-2.0

program fourier_exponential
    !! Temperature-dependent Fourier conduction.
    !! SI example: reference conductivity [W/(m*K)], temperature coefficient [1/K] and reference temperature [K].
    !! Inputs are temperature gradient [K/m] and absolute temperature [K]. Heat flux [W/m2].
    use constituf, only: wp, material_t, cf_ok
    use constituf_fourier_exponential, only: fourier_exponential_t
    implicit none (type, external)

    integer, parameter :: ninput = 4, noutput = 3
    type(material_t) :: material
    real(wp) :: driving(ninput,1), heat_flux(noutput,1)
    real(wp) :: tangent(noutput,ninput,1), dsdp(noutput,1,1)
    integer :: status

    ! Select one exact analytical parameter derivative; other names can be added.
    call material%init(fourier_exponential_t(10.0_wp,0.01_wp,300.0_wp),1,status, &
        sensitivity_parameters=[character(len=23) :: "temperature_coefficient"])
    if (status /= cf_ok) error stop "Initialization failed"

    ! Input: temperature gradient [K/m] and absolute temperature [K].
    driving(:,1) = [1.5_wp,0.0_wp,-0.5_wp,320.0_wp]
    call material%eval(driving,heat_flux,tangent,status,dsdp=dsdp)
    if (status /= cf_ok) error stop "Evaluation failed"

    print "(a,*(1x,es12.4))", "Heat flux [W/m2]:", heat_flux(:,1)
    print "(a,*(1x,es12.4))", "Tangent: row 1, columns 1:3 [W/(m*K)]:", tangent(1,1:3,1)
    print "(a,*(1x,es12.4))", "Tangent: row 1, remaining columns [W/(m2*K)]:", tangent(1,4:,1)
    print "(a,*(1x,es12.4))", "d(response)/d(temperature_coefficient) [W*K/m2]:", dsdp(:,1,1)
end program fourier_exponential

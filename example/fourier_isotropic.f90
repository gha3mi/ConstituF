! SPDX-License-Identifier: Apache-2.0

program fourier_isotropic
    !! Isotropic Fourier conduction.
    !! SI example: conductivity [W/(m*K)], temperature gradient [K/m] and heat flux [W/m2].
    use constituf, only: wp, material_t, cf_ok
    use constituf_fourier_isotropic, only: fourier_isotropic_t
    implicit none (type, external)

    integer, parameter :: ninput = 3, noutput = 3
    type(material_t) :: material
    real(wp) :: temperature_gradient(ninput,1), heat_flux(noutput,1)
    real(wp) :: tangent(noutput,ninput,1), dsdp(noutput,1,1)
    integer :: status

    ! Select one exact analytical parameter derivative; other names can be added.
    call material%init(fourier_isotropic_t(25.0_wp),1,status, &
        sensitivity_parameters=[character(len=12) :: "conductivity"])
    if (status /= cf_ok) error stop "Initialization failed"

    temperature_gradient(:,1) = [2.0_wp,0.0_wp,-1.0_wp]
    call material%eval(temperature_gradient,heat_flux,tangent,status,dsdp=dsdp)
    if (status /= cf_ok) error stop "Evaluation failed"

    print "(a,*(1x,es12.4))", "Heat flux [W/m2]:", heat_flux(:,1)
    print "(a,*(1x,es12.4))", "Tangent: first output row [W/(m*K)]:", tangent(1,:,1)
    print "(a,*(1x,es12.4))", "d(response)/d(conductivity) [K/m]:", dsdp(:,1,1)
end program fourier_isotropic

! SPDX-License-Identifier: Apache-2.0

program fick_anisotropic
    !! Anisotropic Fick diffusion.
    !! SI example: diffusivity [m2/s], molar concentration gradient [mol/m4] and molar flux [mol/(m2*s)].
    use constituf, only: wp, material_t, cf_ok
    use constituf_fick_anisotropic, only: fick_anisotropic_t
    implicit none (type, external)

    integer, parameter :: ninput = 3, noutput = 3
    type(material_t) :: material
    real(wp) :: concentration_gradient(ninput,1), molar_flux(noutput,1)
    real(wp) :: tangent(noutput,ninput,1), dsdp(noutput,1,1)
    integer :: status

    ! Select one exact analytical parameter derivative; other names can be added.
    call material%init(fick_anisotropic_t(3.0e-9_wp,2.0e-9_wp,1.0e-9_wp,0.2e-9_wp,-0.1e-9_wp,0.15e-9_wp), &
        1,status, &
        sensitivity_parameters=[character(len=14) :: "diffusivity_12"])
    if (status /= cf_ok) error stop "Initialization failed"

    concentration_gradient(:,1) = [1.0e2_wp,-2.0e2_wp,0.5e2_wp]
    call material%eval(concentration_gradient,molar_flux,tangent,status,dsdp=dsdp)
    if (status /= cf_ok) error stop "Evaluation failed"

    print "(a,*(1x,es12.4))", "Molar flux [mol/(m2*s)]:", molar_flux(:,1)
    print "(a,*(1x,es12.4))", "Tangent: first output row [m2/s]:", tangent(1,:,1)
    print "(a,*(1x,es12.4))", "d(response)/d(diffusivity_12) [mol/m4]:", dsdp(:,1,1)
end program fick_anisotropic

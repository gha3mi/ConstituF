! SPDX-License-Identifier: Apache-2.0

program ohmic_anisotropic
    !! Anisotropic Ohmic conduction.
    !! SI example: conductivity [S/m], electric field [V/m] and current density [A/m2].
    use constituf, only: wp, material_t, cf_ok
    use constituf_ohmic_anisotropic, only: ohmic_anisotropic_t
    implicit none (type, external)

    integer, parameter :: ninput = 3, noutput = 3
    type(material_t) :: material
    real(wp) :: electric_field(ninput,1), current_density(noutput,1)
    real(wp) :: tangent(noutput,ninput,1), dsdp(noutput,1,1)
    integer :: status

    ! Select one exact analytical parameter derivative; other names can be added.
    call material%init(ohmic_anisotropic_t(4.0_wp,3.0_wp,2.0_wp,0.5_wp,0.0_wp,0.0_wp),1,status, &
        sensitivity_parameters=[character(len=15) :: "conductivity_12"])
    if (status /= cf_ok) error stop "Initialization failed"

    electric_field(:,1) = [1.0_wp,-2.0_wp,0.0_wp]
    call material%eval(electric_field,current_density,tangent,status,dsdp=dsdp)
    if (status /= cf_ok) error stop "Evaluation failed"

    print "(a,*(1x,es12.4))", "Current density [A/m2]:", current_density(:,1)
    print "(a,*(1x,es12.4))", "Tangent: first output row [S/m]:", tangent(1,:,1)
    print "(a,*(1x,es12.4))", "d(response)/d(conductivity_12) [V/m]:", dsdp(:,1,1)
end program ohmic_anisotropic

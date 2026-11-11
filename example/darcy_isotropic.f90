! SPDX-License-Identifier: Apache-2.0

program darcy_isotropic
    !! Isotropic Darcy flow.
    !! SI example: permeability [m2], viscosity [Pa*s], pressure gradient [Pa/m] and velocity [m/s].
    use constituf, only: wp, material_t, cf_ok
    use constituf_darcy_isotropic, only: darcy_isotropic_t
    implicit none (type, external)

    integer, parameter :: ninput = 3, noutput = 3
    type(material_t) :: material
    real(wp) :: pressure_gradient(ninput,1), velocity(noutput,1)
    real(wp) :: tangent(noutput,ninput,1), dsdp(noutput,1,1)
    integer :: status

    ! Select one exact analytical parameter derivative; other names can be added.
    call material%init(darcy_isotropic_t(1.0e-12_wp,1.0e-3_wp),1,status, &
        sensitivity_parameters=[character(len=12) :: "permeability"])
    if (status /= cf_ok) error stop "Initialization failed"

    pressure_gradient(:,1) = [1.0e6_wp,0.0_wp,-0.5e6_wp]
    call material%eval(pressure_gradient,velocity,tangent,status,dsdp=dsdp)
    if (status /= cf_ok) error stop "Evaluation failed"

    print "(a,*(1x,es12.4))", "Darcy velocity [m/s]:", velocity(:,1)
    print "(a,*(1x,es12.4))", "Tangent: first output row [m2/(Pa*s)]:", tangent(1,:,1)
    print "(a,*(1x,es12.4))", "d(response)/d(permeability) [1/(m*s)]:", dsdp(:,1,1)
end program darcy_isotropic

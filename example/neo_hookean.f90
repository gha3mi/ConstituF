! SPDX-License-Identifier: Apache-2.0

program neo_hookean
    !! Neo-Hookean finite elasticity.
    !! SI example: elastic coefficients and first Piola stress [Pa]. Deformation gradient [-].
    use constituf, only: wp, material_t, cf_ok
    use constituf_neo_hookean, only: neo_hookean_t
    implicit none (type, external)

    integer, parameter :: ninput = 9, noutput = 9
    type(material_t) :: material
    real(wp) :: deformation(ninput,1), piola(noutput,1)
    real(wp) :: tangent(noutput,ninput,1), dsdp(noutput,1,1)
    integer :: status

    ! Select one exact analytical parameter derivative; other names can be added.
    call material%init(neo_hookean_t(80.0e6_wp,40.0e6_wp),1,status, &
        sensitivity_parameters=[character(len=13) :: "shear_modulus"])
    if (status /= cf_ok) error stop "Initialization failed"

    ! Deformation gradient F in column-major order: 11,21,31,12,22,32,13,23,33.
    deformation(:,1) = [1.1_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp]
    call material%eval(deformation,piola,tangent,status,dsdp=dsdp)
    if (status /= cf_ok) error stop "Evaluation failed"

    print "(a,*(1x,es12.4))", "First Piola stress [Pa]:", piola(:,1)
    print "(a,*(1x,es12.4))", "Tangent: first output row [Pa]:", tangent(1,:,1)
    print "(a,*(1x,es12.4))", "d(response)/d(shear_modulus) [-]:", dsdp(:,1,1)
end program neo_hookean

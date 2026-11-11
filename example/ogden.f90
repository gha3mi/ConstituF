! SPDX-License-Identifier: Apache-2.0

program ogden
    !! Ogden finite elasticity.
    !! SI example: moduli and first Piola stress [Pa]. Exponent and deformation gradient [-].
    use constituf, only: wp, material_t, cf_ok
    use constituf_ogden, only: ogden_t
    implicit none (type, external)

    integer, parameter :: ninput = 9, noutput = 9
    type(material_t) :: material
    real(wp) :: deformation(ninput,1), piola(noutput,1)
    real(wp) :: tangent(noutput,ninput,1), dsdp(noutput,1,1)
    integer :: status

    ! Select one exact analytical parameter derivative; other names can be added.
    call material%init(ogden_t(1.5e6_wp,1.4_wp,1.0e8_wp),1,status, &
        sensitivity_parameters=[character(len=8) :: "exponent"])
    if (status /= cf_ok) error stop "Initialization failed"

    ! Deformation gradient F in column-major order: 11,21,31,12,22,32,13,23,33.
    deformation(:,1) = [1.10_wp,0.0_wp,0.0_wp,0.15_wp,0.97_wp,0.0_wp,0.0_wp,0.0_wp,1.02_wp]
    call material%eval(deformation,piola,tangent,status,dsdp=dsdp)
    if (status /= cf_ok) error stop "Evaluation failed"

    print "(a,*(1x,es12.4))", "First Piola stress [Pa]:", piola(:,1)
    print "(a,*(1x,es12.4))", "Tangent: first output row [Pa]:", tangent(1,:,1)
    print "(a,*(1x,es12.4))", "d(response)/d(exponent) [Pa]:", dsdp(:,1,1)
end program ogden

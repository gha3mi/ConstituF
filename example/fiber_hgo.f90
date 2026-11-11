! SPDX-License-Identifier: Apache-2.0

program fiber_hgo
    !! Fiber-reinforced finite elasticity.
    !! SI example: shear modulus, Lame lambda and fiber k1 [Pa], fiber k2 [-].
    !! Deformation gradient [-], first Piola stress [Pa] and energy density [J/m3].
    use constituf, only: wp, material_t, cf_ok
    use constituf_fiber_hgo, only: fiber_hgo_t
    implicit none (type, external)

    integer, parameter :: ninput = 9, noutput = 10
    type(material_t) :: material
    real(wp) :: deformation(ninput,1), response(noutput,1)
    real(wp) :: tangent(noutput,ninput,1), dsdp(noutput,1,1)
    integer :: status

    ! Select one exact analytical parameter derivative; other names can be added.
    call material%init(fiber_hgo_t(6.0e6_wp,8.0e6_wp,4.0e6_wp,3.0_wp),1,status, &
        sensitivity_parameters=[character(len=8) :: "fiber_k1"])
    if (status /= cf_ok) error stop "Initialization failed"

    ! Deformation gradient F in column-major order: 11,21,31,12,22,32,13,23,33.
    deformation(:,1) = [1.2_wp,0.0_wp,0.0_wp,0.0_wp,0.95_wp,0.0_wp,0.0_wp,0.0_wp,1.02_wp]
    call material%eval(deformation,response,tangent,status,dsdp=dsdp)
    if (status /= cf_ok) error stop "Evaluation failed"

    print "(a,*(1x,es12.4))", "First Piola stress [Pa]:", response(1:9,1)
    print "(a,*(1x,es12.4))", "Helmholtz energy [J/m3]:", response(10,1)
    print "(a,*(1x,es12.4))", "Tangent: first output row [Pa]:", tangent(1,:,1)
    print "(a,*(1x,es12.4))", "d(response)/d(fiber_k1) [-]:", dsdp(:,1,1)
end program fiber_hgo

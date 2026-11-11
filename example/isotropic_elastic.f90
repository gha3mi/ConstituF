! SPDX-License-Identifier: Apache-2.0

program isotropic_elastic
    !! Isotropic small-strain elasticity.
    !! SI example: bulk/shear moduli and stresses [Pa]. Strains [-].
    use constituf, only: wp, material_t, cf_ok
    use constituf_isotropic_elastic, only: isotropic_elastic_t
    implicit none (type, external)

    integer, parameter :: ninput = 6, noutput = 6
    type(material_t) :: material
    real(wp) :: strain(ninput,1), stress(noutput,1)
    real(wp) :: tangent(noutput,ninput,1), dsdp(noutput,1,1)
    integer :: status

    ! Select one exact analytical parameter derivative; other names can be added.
    call material%init(isotropic_elastic_t(160.0e9_wp,80.0e9_wp),1,status, &
        sensitivity_parameters=[character(len=13) :: "shear_modulus"])
    if (status /= cf_ok) error stop "Initialization failed"

    strain(:,1) = [0.001_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp]
    call material%eval(strain,stress,tangent,status,dsdp=dsdp)
    if (status /= cf_ok) error stop "Evaluation failed"

    print "(a,*(1x,es12.4))", "Cauchy stress [Pa]:", stress(:,1)
    print "(a,*(1x,es12.4))", "Tangent: first output row [Pa]:", tangent(1,:,1)
    print "(a,*(1x,es12.4))", "d(response)/d(shear_modulus) [-]:", dsdp(:,1,1)
end program isotropic_elastic

! SPDX-License-Identifier: Apache-2.0

program orthotropic_elastic
    !! Orthotropic elasticity in principal axes.
    !! SI example: Young/shear moduli and stresses [Pa]. Poisson ratios and strains [-].
    use constituf, only: wp, material_t, cf_ok
    use constituf_orthotropic_elastic, only: orthotropic_elastic_t
    implicit none (type, external)

    integer, parameter :: ninput = 6, noutput = 6
    type(material_t) :: material
    real(wp) :: strain(ninput,1), stress(noutput,1)
    real(wp) :: tangent(noutput,ninput,1), dsdp(noutput,1,1)
    integer :: status

    ! Select one exact analytical parameter derivative; other names can be added.
    call material%init(orthotropic_elastic_t(150.0e9_wp,100.0e9_wp,75.0e9_wp,0.25_wp,0.18_wp,0.15_wp, &
        35.0e9_wp,40.0e9_wp,45.0e9_wp),1,status, &
        sensitivity_parameters=[character(len=7) :: "young_1"])
    if (status /= cf_ok) error stop "Initialization failed"

    strain(:,1) = [0.001_wp,0.0_wp,0.0_wp,0.0002_wp,0.0_wp,0.0_wp]
    call material%eval(strain,stress,tangent,status,dsdp=dsdp)
    if (status /= cf_ok) error stop "Evaluation failed"

    print "(a,*(1x,es12.4))", "Cauchy stress [Pa]:", stress(:,1)
    print "(a,*(1x,es12.4))", "Tangent: first output row [Pa]:", tangent(1,:,1)
    print "(a,*(1x,es12.4))", "d(response)/d(young_1) [-]:", dsdp(:,1,1)
end program orthotropic_elastic

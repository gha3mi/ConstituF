! SPDX-License-Identifier: Apache-2.0

program test_core_verification
    !! General-purpose verification on elastic and history-bearing laws.
    !! Mechanical moduli and stresses use Pa. Time increments and relaxation times use s.
    use constituf, only: wp, material_t, verify_tangent, cf_ok, cf_invalid_shape, cf_invalid_argument
    use constituf_isotropic_elastic, only: isotropic_elastic_t
    use constituf_maxwell_shear, only: maxwell_shear_t
    implicit none (type, external)

    integer, parameter :: nmandel = 6

    type(material_t) :: material
    real(wp) :: x(nmandel,2), y(nmandel,2), initial(nmandel,2)
    real(wp) :: err
    integer :: status, i

    x = 0.0_wp
    x(1,:) = [0.03_wp,-0.02_wp]
    x(4,:) = [0.006_wp,0.003_wp]
    call material%init(isotropic_elastic_t(80.0_wp,40.0_wp),2,status)
    if (status /= cf_ok) error stop "Verification elastic setup"
    call verify_tangent(material,x,err,status)
    if (status /= cf_ok .or. err > 1.0_wp) error stop "Elastic tangent verification"
    call verify_tangent(material,x,err,status,step=[(0.0_wp,i=1,6)])
    if (status /= cf_invalid_argument) error stop "Invalid perturbations accepted"
    call verify_tangent(material,x,err,status,step=[1.0e-5_wp,1.0e-5_wp])
    if (status /= cf_invalid_shape) error stop "Wrong perturbation shape accepted"

    initial = 0.0_wp
    call material%init(maxwell_shear_t(90.0_wp,20.0_wp,35.0_wp,1.5_wp),2,status,initial)
    if (status /= cf_ok) error stop "Verification Maxwell setup"
    call material%eval(x,y,status=status,dt=0.3_wp)
    if (status /= cf_ok) error stop "Verification Maxwell accepted step"
    call material%commit(status)
    if (status /= cf_ok) error stop "Verification Maxwell commit"
    x(1,:) = 1.2_wp*x(1,:)
    call verify_tangent(material,x,err,status,dt=0.2_wp)
    if (status /= cf_ok .or. err > 1.0_wp) error stop "History-dependent tangent verification"
    call material%eval(x,y,status=status,dt=0.2_wp)
    if (status /= cf_ok) error stop "Verification did not preserve committed history"
end program test_core_verification

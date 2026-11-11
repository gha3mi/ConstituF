! SPDX-License-Identifier: Apache-2.0

program test_core_potential
    !! Verify Helmholtz and electric enthalpy conjugacy with correct signs.
    !! Test data use SI units, including J/m3 for potentials and V/m for electric fields.
    use constituf, only: wp, material_t, verify_potential, cf_ok, cf_invalid_argument, cf_invalid_shape
    use constituf_fiber_hgo, only: fiber_hgo_t
    use constituf_piezoelectric_linear, only: piezoelectric_linear_t
    implicit none (type, external)

    integer, parameter :: nfull = 9, npiezo = nfull, nstrain = 6, nresponse = nfull+1
    real(wp), parameter :: difference_scale = epsilon(1.0_wp)**(1.0_wp/3.0_wp)
    real(wp), parameter :: strain_step = max(1.0e-6_wp,0.1_wp*difference_scale)
    real(wp), parameter :: electric_step = max(0.01_wp,1000.0_wp*difference_scale)
    type(material_t) :: material
    real(wp) :: x(nfull,1), y(nresponse,1)
    real(wp) :: err
    integer :: status, i

    x(:,1) = [1.12_wp,0.0_wp,0.0_wp,0.0_wp,0.94_wp,0.0_wp,0.0_wp,0.0_wp,1.04_wp]
    call material%init(fiber_hgo_t(6.0_wp,8.0_wp,4.0_wp,3.0_wp),1,status)
    if (status /= cf_ok) error stop "Potential HGO initialization"
    ! Centered differences balance quadratic truncation against roundoff/h.
    ! These perturbations remain clear of both fiber activation switches.
    call verify_potential(material,x,10,[(i,i=1,nfull)],[(1.0_wp,i=1,nfull)],err,status, &
        relative_tolerance=max(5.0e-5_wp,32.0_wp*difference_scale**2), &
        step=[(max(1.0e-6_wp,difference_scale),i=1,nfull)])
    if (status /= cf_ok .or. err > 1.0_wp) error stop "HGO potential work conjugacy"
    call verify_potential(material,x,11,[(i,i=1,nfull)],[(1.0_wp,i=1,nfull)],err,status)
    if (status /= cf_invalid_argument) error stop "Unknown potential index accepted"
    call verify_potential(material,x,10,[(i,i=1,nstrain)],[(1.0_wp,i=1,nfull)],err,status)
    if (status /= cf_invalid_shape) error stop "Wrong energy conjugate size accepted"
    call verify_potential(material,x,10,[(i,i=1,nfull)],[(0.0_wp,i=1,nfull)],err,status)
    if (status /= cf_invalid_argument) error stop "Zero conjugate signs accepted"

    x(:,1) = [0.004_wp,-0.001_wp,0.0005_wp,0.0008_wp,0.0002_wp,-0.0005_wp,200.0_wp,100.0_wp,-500.0_wp]
    call material%init(piezoelectric_linear_t(150000.0_wp,75000.0_wp,-2.3_wp,4.2_wp,3.5_wp, &
        1.2e-8_wp,2.0e-8_wp),1,status)
    if (status /= cf_ok) error stop "Potential piezoelectric initialization"
    ! Quadratic enthalpy permits finite steps without truncation error.
    ! Strain (1) and electric-field (V/m) steps retain distinct physical scales.
    call verify_potential(material,x,10,[(i,i=1,npiezo)], &
        [(1.0_wp,i=1,nstrain),(-1.0_wp,i=nstrain+1,npiezo)],err,status, &
        relative_tolerance=1.0e-4_wp,step=[(strain_step,i=1,nstrain),electric_step,electric_step,electric_step])
    if (status /= cf_ok .or. err > 1.0_wp) error stop "Piezoelectric signed enthalpy conjugacy"

    call material%eval(x,y,status=status)
    if (status /= cf_ok) error stop "Potential verification damaged material usage"
end program test_core_potential

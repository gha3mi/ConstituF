! SPDX-License-Identifier: Apache-2.0

program test_core_calibration
    !! Recover elastic moduli from multiaxial stress observations using
    !! full native analytical Jacobians, bounds and scaled measurement errors.
    use constituf, only: wp, material_t, calibrate_parameters, cf_ok, cf_invalid_argument, cf_invalid_shape, cf_evaluation_failed
    use constituf_isotropic_elastic, only: isotropic_elastic_t
    implicit none (type, external)

    integer, parameter :: nmandel = 6, npoints = 4, nparameters = 2
    integer, parameter :: nobs = nmandel*npoints
    real(wp), parameter :: truth(nparameters) = [150000.0_wp,75000.0_wp]
    type(material_t) :: reference
    real(wp) :: prescribed(nmandel,npoints), stress(nmandel,npoints)
    real(wp) :: target(nobs), noise(nobs), fitted(nparameters)
    real(wp) :: cost
    real(wp) :: starting(nparameters), bounds_lo(nparameters), bounds_hi(nparameters)
    integer :: status, steps

    prescribed(:,1) = [0.012_wp,0.002_wp,-0.001_wp,0.001_wp,0.0003_wp,0.0_wp]
    prescribed(:,2) = [-0.002_wp,0.009_wp,-0.001_wp,0.0_wp,0.0008_wp,0.001_wp]
    prescribed(:,3) = [0.003_wp,0.002_wp,0.006_wp,-0.0004_wp,0.0002_wp,-0.0005_wp]
    prescribed(:,4) = [0.001_wp,-0.005_wp,0.007_wp,0.0005_wp,-0.0007_wp,0.0003_wp]
    call reference%init(isotropic_elastic_t(truth(1),truth(2)),npoints,status)
    if (status /= cf_ok) error stop "Elastic calibration truth setup"
    call reference%eval(prescribed,stress,status=status)
    if (status /= cf_ok) error stop "Elastic calibration truth evaluation"
    target = reshape(stress,[nobs])
    noise = 1.0_wp
    starting = [95000.0_wp,45000.0_wp]
    bounds_lo = [10000.0_wp,10000.0_wp]
    bounds_hi = [300000.0_wp,200000.0_wp]

    call calibrate_parameters(starting,target,noise,bounds_lo,bounds_hi,evaluate_elastic,fitted,cost,status, &
        parameter_scale=[100000.0_wp,60000.0_wp],max_iterations=80,iterations=steps)
    if (status /= cf_ok) error stop "Exact-Jacobian elastic calibration did not converge"
    if (maxval(abs(fitted-truth)) > 1.0e-4_wp) error stop "Calibrated elastic moduli differ from target"
    if (cost > 1.0e-10_wp) error stop "Calibration residual was not reduced"
    if (steps < 1) error stop "Calibration failed to accept any iteration"

    noise(1) = 0.0_wp
    call calibrate_parameters(starting,target,noise,bounds_lo,bounds_hi,evaluate_elastic,fitted,cost,status)
    if (status /= cf_invalid_argument) error stop "Nonpositive observation uncertainty accepted"
    noise(1) = 1.0_wp
    call calibrate_parameters(starting,target(1:nparameters-1),noise,bounds_lo,bounds_hi, &
        evaluate_elastic,fitted,cost,status)
    if (status /= cf_invalid_shape) error stop "Underdetermined observation count accepted"
    call calibrate_parameters(starting,target,noise,bounds_hi,bounds_lo,evaluate_elastic,fitted,cost,status)
    if (status /= cf_invalid_argument) error stop "Inverted bounds accepted"
    call calibrate_parameters(starting,target,noise,bounds_lo,bounds_hi,evaluate_elastic,fitted,cost,status, &
        parameter_scale=[1.0_wp,0.0_wp])
    if (status /= cf_invalid_argument) error stop "Zero parameter scaling accepted"
    call calibrate_parameters(starting,target,noise,bounds_lo,bounds_hi,evaluate_elastic,fitted,cost,status, &
        max_iterations=0)
    if (status /= cf_invalid_argument) error stop "Zero maximum iteration count accepted"
    call calibrate_parameters(starting,target,noise,bounds_lo,bounds_hi,evaluate_bad,fitted,cost,status)
    if (status /= cf_evaluation_failed) error stop "Invalid forward model was not rejected"

contains

    subroutine evaluate_elastic(parameters,predicted,jacobian,status)
        !! Model callback, with a new material_t for each parameter trial.
        !! All data are in Pa, all strains fixed and dimensionless.
        real(wp), intent(in) :: parameters(:)
        real(wp), intent(out) :: predicted(:),jacobian(:,:)
        integer, intent(out) :: status
        type(material_t) :: handle
        real(wp) :: local_stress(nmandel,npoints), dsdp(nmandel,nparameters,npoints)
        integer :: j

        call handle%init(isotropic_elastic_t(parameters(1),parameters(2)),npoints,status, &
            sensitivity_parameters=[character(len=13) :: "bulk_modulus", "shear_modulus"])
        if (status /= cf_ok) return
        call handle%eval(prescribed,local_stress,status=status,dsdp=dsdp)
        if (status /= cf_ok) return
        predicted = reshape(local_stress,[nobs])
        do j = 1,nparameters
            jacobian(:,j) = reshape(dsdp(:,j,:),[nobs])
        end do
    end subroutine evaluate_elastic

    subroutine evaluate_bad(parameters,predicted,jacobian,status)
        real(wp), intent(in) :: parameters(:)
        real(wp), intent(out) :: predicted(:),jacobian(:,:)
        integer, intent(out) :: status

        predicted = sum(parameters)
        jacobian = 0.0_wp
        status = cf_invalid_argument
    end subroutine evaluate_bad

end program test_core_calibration

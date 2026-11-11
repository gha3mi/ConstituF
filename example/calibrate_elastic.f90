! SPDX-License-Identifier: Apache-2.0

program calibrate_elastic
    !! Offline nonlinear least-squares recovery of isotropic elastic moduli.
    !! The callback evaluates the native Fortran material and its exact
    !! analytical parameter Jacobian; no finite differences enter the fit.
    !! SI example: moduli, stresses, observation noise and parameter bounds/scales [Pa]. Strains [-].
    use constituf, only: wp, material_t, calibrate_parameters, cf_ok
    use constituf_isotropic_elastic, only: isotropic_elastic_t
    implicit none (type, external)

    integer, parameter :: nmandel = 6, npoints = 3, nparameters = 2
    integer, parameter :: nobs = nmandel*npoints
    type(material_t) :: truth_material
    real(wp) :: strain(nmandel,npoints), stress(nmandel,npoints)
    real(wp) :: target(nobs), noise(nobs), fitted(nparameters)
    real(wp) :: cost
    integer :: status

    strain(:,1) = [0.010_wp,0.002_wp,0.0_wp,0.001_wp,0.0_wp,0.0_wp]
    strain(:,2) = [0.0_wp,0.004_wp,-0.006_wp,0.0_wp,0.0005_wp,0.0_wp]
    strain(:,3) = [0.004_wp,0.003_wp,0.002_wp,0.0_wp,0.0_wp,0.001_wp]
    call truth_material%init(isotropic_elastic_t(150000.0_wp,75000.0_wp),npoints,status)
    if (status /= cf_ok) error stop "Calibration synthetic data setup"
    call truth_material%eval(strain,stress,status=status)
    if (status /= cf_ok) error stop "Calibration synthetic data"
    target = reshape(stress,[nobs])
    noise = 1.0_wp

    call calibrate_parameters([95000.0_wp,45000.0_wp],target,noise, &
        [10000.0_wp,10000.0_wp],[300000.0_wp,200000.0_wp],evaluate,fitted,cost,status, &
        parameter_scale=[100000.0_wp,60000.0_wp])
    if (status /= cf_ok) error stop "Elastic calibration did not converge"
    print "(a,2(1x,es12.4))", "Recovered K,G [Pa]: ",fitted
    print "(a,es12.4)", "Weighted least-squares objective: ",cost

contains

    subroutine evaluate(parameters,predicted,jacobian,status)
        real(wp), intent(in) :: parameters(:)
        real(wp), intent(out) :: predicted(:),jacobian(:,:)
        integer, intent(out) :: status
        type(material_t) :: material
        real(wp) :: y(nmandel,npoints), dsdp(nmandel,nparameters,npoints)
        integer :: j

        call material%init(isotropic_elastic_t(parameters(1),parameters(2)),npoints,status, &
            sensitivity_parameters=[character(len=13) :: "bulk_modulus", "shear_modulus"])
        if (status /= cf_ok) return
        call material%eval(strain,y,status=status,dsdp=dsdp)
        if (status /= cf_ok) return
        predicted = reshape(y,[nobs])
        do j = 1,nparameters
            jacobian(:,j) = reshape(dsdp(:,j,:),[nobs])
        end do
    end subroutine evaluate

end program calibrate_elastic

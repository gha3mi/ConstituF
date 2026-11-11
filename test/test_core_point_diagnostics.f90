! SPDX-License-Identifier: Apache-2.0

program test_core_point_diagnostics
    !! Per-point kernel statuses are available only after a point evaluation.
    !! Elastic moduli and stresses in this synthetic test use Pa.
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use constituf, only: wp, material_t, cf_ok, cf_uninitialized, cf_no_trial, cf_invalid_argument, &
                        cf_evaluation_failed, cf_invalid_shape
    use constituf_isotropic_elastic, only: isotropic_elastic_t
    implicit none (type, external)

    integer, parameter :: nmandel = 6, npoints = 3
    type(material_t) :: material, empty
    real(wp) :: x(nmandel,npoints), y(nmandel,npoints)
    integer :: statuses(npoints), wrong(npoints-1)
    integer :: status

    call empty%point_statuses(statuses,status)
    if (status /= cf_uninitialized) error stop "Empty material diagnostics accepted"
    call material%init(isotropic_elastic_t(100.0_wp,50.0_wp),npoints,status)
    if (status /= cf_ok) error stop "Point diagnostics setup"
    call material%point_statuses(statuses,status)
    if (status /= cf_no_trial) error stop "Diagnostics present before initial trial"
    x = 0.0_wp
    x(1,:) = [0.001_wp,0.002_wp,0.003_wp]
    x(1,2) = ieee_value(0.0_wp,ieee_quiet_nan)
    call material%eval(x,y,status=status)
    if (status /= cf_evaluation_failed) error stop "One point must reject batch"
    call material%point_statuses(statuses,status)
    if (status /= cf_ok) error stop "Point diagnostics unavailable after failure"
    if (any(statuses /= [cf_ok,cf_invalid_argument,cf_ok])) error stop "Per-point status mismatch"
    if (material%last_failure_point() /= 2) error stop "First failing point mismatch"
    call material%commit(status)
    if (status /= cf_no_trial) error stop "Failed point trial became committed"
    call material%point_statuses(wrong,status)
    if (status /= cf_invalid_shape) error stop "Wrong diagnostic buffer size accepted"

    x(1,2) = 0.002_wp
    call material%eval(x,y,status=status)
    if (status /= cf_ok) error stop "Healthy batch evaluation"
    call material%point_statuses(statuses,status)
    if (status /= cf_ok .or. any(statuses /= cf_ok)) error stop "Stale kernel diagnostics after success"

    call material%eval(x(:,1:2),y,status=status)
    if (status /= cf_invalid_shape) error stop "Wrong input point count accepted"
    call material%point_statuses(statuses,status)
    if (status /= cf_no_trial) error stop "Stale kernel diagnostics after setup error"
end program test_core_point_diagnostics

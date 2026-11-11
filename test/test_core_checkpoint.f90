! SPDX-License-Identifier: Apache-2.0

program test_core_checkpoint
    !! Committed-only snapshots, typed history and sensitivity restart.
    !! Mechanical data use Pa and s. The synthetic typed model is dimensionless.
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use constituf, only: wp, ik, material_t, cf_ok, cf_invalid_shape, cf_invalid_argument
    use constituf_maxwell_shear, only: maxwell_shear_t
    use constituf_test_typed, only: typed_test_model_t
    implicit none (type, external)

    integer, parameter :: nmandel = 6

    type(material_t) :: original, restored, typed, typed_restored
    real(wp) :: initial(nmandel,2), strain(nmandel,2), next_strain(nmandel,2)
    real(wp) :: y(nmandel,2), y_restored(nmandel,2), derivative(nmandel,1,2), derivative_restored(nmandel,1,2)
    real(wp) :: state(nmandel,2), old_input(nmandel,2), sensitivity(nmandel,1,2), exported(nmandel,2), invalid(nmandel,2)
    real(wp) :: typed_input(1,2), typed_output(1,2), typed_reference(1,2)
    real(wp) :: typed_state(1,2)
    integer(ik) :: counts(1,2)
    logical :: seen(1,2)
    integer :: status

    initial = 0.0_wp
    strain = 0.0_wp
    strain(1,:) = [0.004_wp,-0.006_wp]
    strain(2,:) = [-0.004_wp,0.002_wp]
    strain(4,:) = [0.001_wp,0.003_wp]
    call original%init(maxwell_shear_t(80.0_wp,15.0_wp,30.0_wp,2.0_wp),2,status,initial, &
        sensitivity_parameters=[character(len=15) :: "relaxation_time"])
    if (status /= cf_ok) error stop "Checkpoint setup"
    call original%eval(strain,y,status=status,dt=0.3_wp,dsdp=derivative)
    if (status /= cf_ok) error stop "Checkpoint initial path"
    call original%commit(status)
    if (status /= cf_ok) error stop "Checkpoint first commit"

    call original%export_committed(state,status,previous_input=old_input,sensitivity_state=sensitivity)
    if (status /= cf_ok) error stop "Checkpoint export"
    call original%export_committed(exported,status)
    if (status /= cf_invalid_shape) error stop "Missing checkpoint buffers accepted"

    next_strain = 1.2_wp*strain
    call original%eval(next_strain,y,status=status,dt=0.2_wp,dsdp=derivative)
    if (status /= cf_ok) error stop "Pending uncommitted trial"
    call original%export_committed(exported,status,previous_input=old_input,sensitivity_state=sensitivity)
    if (status /= cf_ok .or. maxval(abs(exported-state)) > 0.0_wp) then
        error stop "Export included unaccepted trial state"
    end if

    call restored%init(maxwell_shear_t(80.0_wp,15.0_wp,30.0_wp,2.0_wp),2,status,initial, &
        sensitivity_parameters=[character(len=15) :: "relaxation_time"])
    if (status /= cf_ok) error stop "Checkpoint restored model setup"
    call restored%import_committed(state,status,previous_input=old_input,sensitivity_state=sensitivity)
    if (status /= cf_ok) error stop "Checkpoint import"
    call restored%eval(next_strain,y_restored,status=status,dt=0.2_wp,dsdp=derivative_restored)
    if (status /= cf_ok) error stop "Restored subsequent evaluation"
    if (maxval(abs(y-y_restored)) > 512.0_wp*epsilon(1.0_wp)) error stop "Restarted stress differs"
    if (maxval(abs(derivative-derivative_restored)) > 512.0_wp*epsilon(1.0_wp)) then
        error stop "Restarted parameter sensitivity differs"
    end if

    invalid = state
    invalid(1,1) = ieee_value(0.0_wp,ieee_quiet_nan)
    call restored%import_committed(invalid,status,previous_input=old_input,sensitivity_state=sensitivity)
    if (status /= cf_invalid_argument) error stop "Nonfinite checkpoint accepted"
    call restored%export_committed(exported,status,previous_input=old_input,sensitivity_state=sensitivity)
    if (status /= cf_ok .or. maxval(abs(exported-state)) > 0.0_wp) error stop "Invalid checkpoint changed state"

    call typed%init(typed_test_model_t(),2,status)
    if (status /= cf_ok) error stop "Typed checkpoint model setup"
    typed_input(1,:) = [2.0_wp,-0.5_wp]
    call typed%eval(typed_input,typed_output,status=status)
    if (status /= cf_ok) error stop "Typed checkpoint trial"
    call typed%commit(status)
    if (status /= cf_ok) error stop "Typed checkpoint commit"
    call typed%export_committed(typed_state,status,integer_state=counts,logical_state=seen)
    if (status /= cf_ok) error stop "Typed checkpoint export"
    if (any(counts(1,:) /= [2_ik,1_ik])) error stop "Typed checkpoint integer history"
    if (.not. seen(1,1) .or. seen(1,2)) error stop "Typed checkpoint logical history"

    call typed_restored%init(typed_test_model_t(),2,status)
    if (status /= cf_ok) error stop "Typed checkpoint restore setup"
    call typed_restored%import_committed(typed_state,status,integer_state=counts,logical_state=seen)
    if (status /= cf_ok) error stop "Typed checkpoint restore"
    typed_input(1,:) = [-1.0_wp,1.0_wp]
    call typed%eval(typed_input,typed_output,status=status)
    call typed_restored%eval(typed_input,typed_reference,status=status)
    if (status /= cf_ok) error stop "Typed checkpoint continuation"
    if (maxval(abs(typed_output-typed_reference)) > 0.0_wp) error stop "Typed checkpoint changed response"
end program test_core_checkpoint

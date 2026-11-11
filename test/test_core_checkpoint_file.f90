! SPDX-License-Identifier: Apache-2.0

program test_core_checkpoint_file
    !! Versioned portable checkpoint regression for real, discrete and
    !! analytical-sensitivity history, plus failure-atomic loading.
    !! Mechanical data use Pa and s. The synthetic typed model is dimensionless.
    use constituf, only: wp, ik, material_t, save_checkpoint, load_checkpoint, cf_ok, cf_invalid_argument, &
                        cf_uninitialized, cf_checkpoint_mismatch, cf_checkpoint_corrupt
    use constituf_isotropic_elastic, only: isotropic_elastic_t
    use constituf_maxwell_shear, only: maxwell_shear_t
    use constituf_test_typed, only: typed_test_model_t
    implicit none (type, external)

    integer, parameter :: nmandel = 6, npoints = 2, nselected = 2
    character(*), parameter :: path = "constituf-test-checkpoint-v1.tmp"
    character(*), parameter :: wrong_path = "constituf-test-checkpoint-corrupt.tmp"
    character(*), parameter :: config = "Maxwell K=80 G_inf=15 G1=30 tau=2; points=left,right"
    type(material_t) :: original, restarted, incompatible, uninitialized, typed, typed_restart, elastic, elastic_restart
    real(wp) :: zero(nmandel,npoints), x(nmandel,npoints), next_x(nmandel,npoints)
    real(wp) :: y(nmandel,npoints), y_restarted(nmandel,npoints), dsdp(nmandel,nselected,npoints)
    real(wp) :: d_restarted(nmandel,nselected,npoints)
    real(wp) :: saved(nmandel,npoints), actual(nmandel,npoints), previous(nmandel,npoints)
    real(wp) :: sensitivity(nmandel,nselected,npoints)
    real(wp) :: typed_x(1,npoints), typed_y(1,npoints), typed_compare(1,npoints), typed_history(1,npoints)
    integer(ik) :: integers(1,npoints)
    logical :: flags(1,npoints)
    integer :: status, unit, io_status, file_size, damage_offset
    character(:), allocatable :: file_text
    character(256) :: io_message

    call save_checkpoint(uninitialized,path,config,status)
    if (status /= cf_uninitialized) error stop "Uninitialized file checkpoint accepted"
    zero = 0.0_wp
    call original%init(maxwell_shear_t(80.0_wp,15.0_wp,30.0_wp,2.0_wp),npoints,status,zero, &
        sensitivity_parameters=[character(len=26) :: "relaxation_time", "branch_shear_modulus"])
    if (status /= cf_ok) error stop "Maxwell checkpoint setup"
    if (original%point_count() /= npoints .or. original%previous_input_size() /= nmandel) then
        error stop "Checkpoint layout inquiries"
    end if
    if (original%sensitivity_parameter_id(1) /= 4 .or. original%sensitivity_parameter_id(2) /= 3) then
        error stop "Selected parameter IDs do not preserve caller order"
    end if
    if (original%sensitivity_parameter_id(3) /= 0) error stop "Out-of-range selected parameter ID"
    x = 0.0_wp
    x(1,:) = [0.010_wp,-0.005_wp]
    x(2,:) = [-0.006_wp,0.004_wp]
    x(4,:) = [0.002_wp,-0.001_wp]
    call original%eval(x,y,status=status,dt=0.25_wp,dsdp=dsdp)
    if (status /= cf_ok) error stop "Maxwell initial trial"
    call original%commit(status)
    if (status /= cf_ok) error stop "Maxwell accepted initial step"
    call original%export_committed(saved,status,previous_input=previous,sensitivity_state=sensitivity)
    if (status /= cf_ok) error stop "In-memory reference snapshot"

    call save_checkpoint(original,path,config,status)
    if (status /= cf_ok) error stop "Portable checkpoint write"

    next_x = 1.2_wp*x
    call original%eval(next_x,y,status=status,dt=0.17_wp,dsdp=dsdp)
    if (status /= cf_ok) error stop "Uncommitted Maxwell trial"

    call load_checkpoint(original,path,"WRONG KEY",status)
    if (status /= cf_checkpoint_mismatch) error stop "Wrong model configuration was not rejected"
    call original%export_committed(actual,status,previous_input=previous,sensitivity_state=sensitivity)
    if (status /= cf_ok) error stop "State export after rejected load"
    if (maxval(abs(actual-saved)) > 0.0_wp) error stop "Wrong-key load corrupted committed history"
    call original%commit(status)
    if (status /= cf_ok) error stop "Rejected checkpoint load destroyed pending trial"

    call restarted%init(maxwell_shear_t(80.0_wp,15.0_wp,30.0_wp,2.0_wp),npoints,status,zero, &
        sensitivity_parameters=[character(len=26) :: "relaxation_time", "branch_shear_modulus"])
    if (status /= cf_ok) error stop "Restarted Maxwell setup"
    call load_checkpoint(restarted,path,config,status)
    if (status /= cf_ok) error stop "Portable Maxwell checkpoint restore"
    call restarted%export_committed(actual,status,previous_input=previous,sensitivity_state=sensitivity)
    if (status /= cf_ok) error stop "Checkpoint restored state export"
    if (maxval(abs(actual-saved)) > 0.0_wp) error stop "Portable state roundtrip differs"
    call restarted%eval(next_x,y_restarted,status=status,dt=0.17_wp,dsdp=d_restarted)
    if (status /= cf_ok) error stop "Maxwell continuation after file restart"
    if (maxval(abs(y_restarted-y)) > 1024.0_wp*epsilon(1.0_wp)) then
        error stop "File restart changed continued Maxwell response"
    end if
    if (maxval(abs(d_restarted-dsdp)) > 1024.0_wp*epsilon(1.0_wp)) then
        error stop "File restart changed analytical sensitivity history"
    end if

    call incompatible%init(maxwell_shear_t(80.0_wp,15.0_wp,30.0_wp,2.0_wp),npoints,status,zero, &
        sensitivity_parameters=[character(len=26) :: "branch_shear_modulus", "relaxation_time"])
    if (status /= cf_ok) error stop "Incompatible parameter-order model setup"
    call load_checkpoint(incompatible,path,config,status)
    if (status /= cf_checkpoint_mismatch) error stop "Selected sensitivity ordering mismatch not rejected"
    call incompatible%init(maxwell_shear_t(80.0_wp,15.0_wp,30.0_wp,2.0_wp),1,status,zero(:,1:1), &
        sensitivity_parameters=[character(len=26) :: "relaxation_time", "branch_shear_modulus"])
    if (status /= cf_ok) error stop "Mismatched point-count setup"
    call load_checkpoint(incompatible,path,config,status)
    if (status /= cf_checkpoint_mismatch) error stop "Point-count incompatibility not rejected"

    ! Change one significant decimal digit without updating the digest.
    ! A syntactically valid but modified payload must never be imported.
    inquire(file=path,size=file_size,iostat=io_status,iomsg=io_message)
    if (io_status /= 0 .or. file_size <= 0) error stop "Checkpoint test file size"
    allocate(character(len=file_size) :: file_text)
    open(newunit=unit,file=path,status="old",action="read",access="stream",form="unformatted", &
         iostat=io_status,iomsg=io_message)
    if (io_status /= 0) error stop "Read checkpoint as bytes"
    read(unit,iostat=io_status,iomsg=io_message) file_text
    if (io_status /= 0) error stop "Read checkpoint bytes"
    close(unit,iostat=io_status,iomsg=io_message)
    if (io_status /= 0) error stop "Close checkpoint byte input"
    damage_offset = index(file_text,"E+")
    if (damage_offset == 0) damage_offset = index(file_text,"E-")
    if (damage_offset < 2) error stop "No scientific value to corrupt in checkpoint"
    if (file_text(damage_offset-1:damage_offset-1) == "5") then
        file_text(damage_offset-1:damage_offset-1) = "6"
    else
        file_text(damage_offset-1:damage_offset-1) = "5"
    end if
    open(newunit=unit,file=wrong_path,status="replace",action="write",access="stream",form="unformatted", &
         iostat=io_status,iomsg=io_message)
    if (io_status /= 0) error stop "Write altered checkpoint"
    write(unit,iostat=io_status,iomsg=io_message) file_text
    if (io_status /= 0) error stop "Write altered checkpoint bytes"
    close(unit,iostat=io_status,iomsg=io_message)
    if (io_status /= 0) error stop "Close altered checkpoint"
    call load_checkpoint(restarted,wrong_path,config,status)
    if (status /= cf_checkpoint_corrupt) error stop "Modified numerical payload bypassed digest"
    call restarted%export_committed(actual,status,previous_input=previous,sensitivity_state=sensitivity)
    if (status /= cf_ok .or. maxval(abs(actual-saved)) > 0.0_wp) then
        error stop "Tampered checkpoint changed the live material"
    end if

    ! A trailing line is corruption even if all physical values are valid.
    open(newunit=unit,file=path,status="old",position="append",action="write",iostat=io_status,iomsg=io_message)
    if (io_status /= 0) error stop "Test checkpoint append"
    write(unit,"(a)") "TRAILING DATA"
    close(unit)
    call load_checkpoint(restarted,path,config,status)
    if (status /= cf_checkpoint_corrupt) error stop "Trailing unverified bytes accepted"
    call restarted%export_committed(actual,status,previous_input=previous,sensitivity_state=sensitivity)
    if (status /= cf_ok .or. maxval(abs(actual-saved)) > 0.0_wp) then
        error stop "Corrupted checkpoint import mutated accepted state"
    end if

    ! A short file, a wrong version, and a nonfinite on-disk value fail.
    open(newunit=unit,file=wrong_path,status="replace",action="write")
    write(unit,"(a)") "ConstituF-checkpoint"
    close(unit)
    call load_checkpoint(restarted,wrong_path,config,status)
    if (status /= cf_checkpoint_corrupt) error stop "Truncated checkpoint accepted"
    open(newunit=unit,file=wrong_path,status="replace",action="write")
    write(unit,"(a)") "ConstituF-checkpoint", "999"
    close(unit)
    call load_checkpoint(restarted,wrong_path,config,status)
    if (status /= cf_checkpoint_mismatch) error stop "Unsupported checkpoint version accepted"
    call save_checkpoint(restarted,wrong_path,"",status)
    if (status /= cf_invalid_argument) error stop "Empty configuration key accepted"
    call save_checkpoint(restarted,wrong_path,"invalid"//achar(10)//"key",status)
    if (status /= cf_invalid_argument) error stop "Multiline configuration key accepted"

    ! Verify typed integer/logical states survive file transfer.
    typed_x(1,:) = [0.5_wp,-1.0_wp]
    call typed%init(typed_test_model_t(),npoints,status)
    if (status /= cf_ok) error stop "Typed checkpoint model setup"
    call typed%eval(typed_x,typed_y,status=status)
    if (status /= cf_ok) error stop "Typed checkpoint trial"
    call typed%commit(status)
    if (status /= cf_ok) error stop "Typed checkpoint commit"
    call typed%export_committed(typed_history,status,integer_state=integers,logical_state=flags)
    if (status /= cf_ok) error stop "Typed checkpoint reference"
    call save_checkpoint(typed,wrong_path,"typed-r1",status)
    if (status /= cf_ok) error stop "Typed checkpoint write"
    call typed_restart%init(typed_test_model_t(),npoints,status)
    if (status /= cf_ok) error stop "Typed checkpoint target setup"
    call load_checkpoint(typed_restart,wrong_path,"typed-r1",status)
    if (status /= cf_ok) error stop "Typed checkpoint restore"
    typed_x(1,:) = [-0.7_wp,1.0_wp]
    call typed%eval(typed_x,typed_y,status=status)
    if (status /= cf_ok) error stop "Typed checkpoint original continuation"
    call typed_restart%eval(typed_x,typed_compare,status=status)
    if (status /= cf_ok) error stop "Typed checkpoint restored continuation"
    if (maxval(abs(typed_y-typed_compare)) > 0.0_wp) error stop "Typed checkpoint states do not match"

    ! A stateless material must still be safely serializable.
    call elastic%init(isotropic_elastic_t(80.0_wp,40.0_wp),1,status)
    if (status /= cf_ok) error stop "Stateless checkpoint setup"
    call save_checkpoint(elastic,wrong_path,"elastic-K80-G40",status)
    if (status /= cf_ok) error stop "Stateless checkpoint write"
    call elastic_restart%init(isotropic_elastic_t(80.0_wp,40.0_wp),1,status)
    if (status /= cf_ok) error stop "Stateless checkpoint target setup"
    call load_checkpoint(elastic_restart,wrong_path,"elastic-K80-G40",status)
    if (status /= cf_ok) error stop "Stateless checkpoint reload"
    if (elastic_restart%previous_input_size() /= 0) error stop "Stateless previous input size"

    call delete_test_file(path)
    call delete_test_file(wrong_path)

contains

    subroutine delete_test_file(filename)
        character(*), intent(in) :: filename
        integer :: local_unit, local_status
        character(256) :: local_message

        open(newunit=local_unit,file=filename,status="old",action="readwrite",iostat=local_status,iomsg=local_message)
        if (local_status /= 0) error stop "Cannot remove checkpoint test file"
        close(local_unit,status="delete",iostat=local_status,iomsg=local_message)
        if (local_status /= 0) error stop "Cannot close/delete checkpoint test file"
    end subroutine delete_test_file

end program test_core_checkpoint_file

! SPDX-License-Identifier: Apache-2.0

module constituf_checkpoint
    !! Portable, versioned checkpoint files for committed material-point state.
    !!
    !! This module is an optional setup/restart utility, not a constitutive
    !! kernel. It performs formatted file I/O, allocates temporary buffers
    !! and delegates state transfer to material_t%export_committed and
    !! %import_committed. It never saves or accepts a pending trial.
    !!
    !! Version 1 is a line-oriented ASCII format. It records the magic,
    !! format version, caller-supplied material configuration identity,
    !! precision (radix digits, decimal range, storage bits), input/output
    !! and history shapes, selected analytical parameter IDs in exact order
    !! and all real/integer/logical/previous-input/sensitivity histories.
    !! Real values use 39 significant decimal digits to round-trip all
    !! supported wp kinds. A rotate/XOR checksum detects accidental changes
    !! to the serialized text. It is NOT a cryptographic authentication code.
    !!
    !! IMPORTANT: Fortran does not provide portable reflection of arbitrary
    !! derived-type model parameters. The mandatory configuration_id is a
    !! non-empty, printable ASCII application-owned identifier, *uniquely*
    !! associated with the model type, parameter values, units, point
    !! ordering and any material conventions. The host must preserve and
    !! check that identity. Reusing it for incompatible models is unsafe.
    !! This module additionally checks precision, shape and selected IDs.
    !! It cannot by itself authenticate a misidentified parameter set.
    !!
    !! Loading is failure-atomic for the live material: every line, the
    !! checksum and the EOF marker are checked before import_committed.
    !! A rejected load leaves committed history and pending trial untouched.
    !! Saving uses STATUS='replace' and is NOT an atomic disk replacement.
    !! Applications requiring crash-atomic files should write to a distinct
    !! path and manage the final filesystem rename themselves.
    !!
    !! A single handle must not be concurrently evaluated or mutated while
    !! these routines operate. The file format is not a security boundary.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: iostat_end
    use constituf_kinds, only: wp, ik, lk
    use constituf_material, only: material_t
    use constituf_status, only: cf_ok, cf_invalid_argument, cf_uninitialized, cf_allocation_failed, cf_checkpoint_corrupt, cf_checkpoint_mismatch, cf_io_error
    implicit none (type, external)
    private

    integer, parameter :: file_version = 1
        !! Current checkpoint serialization format version.
    integer, parameter :: max_identity_length = 256
        !! Maximum number of characters in a checkpoint configuration identifier.
    character(*), parameter :: magic = "ConstituF-checkpoint"
        !! Header identifying a ConstituF checkpoint file.
    character(*), parameter :: end_marker = "END-CONSTITUF-CHECKPOINT"
        !! Marker terminating a complete serialized checkpoint.

    type :: snapshot_t
        !! Owned buffers for all committed history fields in a checkpoint.
        real(wp), allocatable :: state(:,:)
            !! Accepted real material state stored in this checkpoint.
        real(wp), allocatable :: previous(:,:)
            !! Previously committed packed input stored in this checkpoint.
        integer(ik), allocatable :: integers(:,:)
            !! Accepted integer-valued history stored in this checkpoint.
        logical, allocatable :: flags(:,:)
            !! Accepted logical-valued history stored in this checkpoint.
        real(wp), allocatable :: sensitivities(:,:,:)
            !! Accepted parameter derivatives of the real material state.
    end type snapshot_t

    public :: save_checkpoint, load_checkpoint

contains

    subroutine save_checkpoint(material, path, configuration_id, status)
        !! Save the accepted state of an initialized material handle.
        !!
        !! path is the destination file. STATUS='replace' overwrites it.
        !! configuration_id must be an application-stable printable ASCII
        !! token (1..256 characters), unique for the model, its parameter
        !! values, point order and conventions. No identity is inferred.
        !! Returns cf_ok, cf_uninitialized, cf_invalid_argument,
        !! cf_allocation_failed or cf_io_error. The caller must not use an
        !! incomplete checkpoint after a failed write or close operation.
        type(material_t), intent(in) :: material
            !! Material handle whose configuration and history are queried or updated.
        character(*), intent(in) :: path
            !! Filesystem path of the checkpoint file.
        character(*), intent(in) :: configuration_id
            !! Stable application-defined identity for model parameters, units and point ordering.
        integer, intent(out) :: status
            !! Operation status code; `cf_ok` indicates success.
        type(snapshot_t) :: snapshot
        integer :: unit
        integer :: io_status
        integer :: close_status
        integer(lk) :: digest
        character(256) :: line
        character(256) :: io_message

        status = cf_invalid_argument
        if (len_trim(path) == 0 .or. .not. valid_identity(configuration_id)) return
        status = cf_uninitialized
        if (material%point_count() == 0) return

        call prepare_snapshot(material,snapshot,status)
        if (status /= cf_ok) return
        call material%export_committed(snapshot%state,status,previous_input=snapshot%previous, integer_state=snapshot%integers,logical_state=snapshot%flags,sensitivity_state=snapshot%sensitivities)
        if (status /= cf_ok) return

        status = cf_io_error
        open(newunit=unit,file=path,status="replace",action="write",form="formatted",iostat=io_status,iomsg=io_message)
        if (io_status /= 0) return
        digest = 0_lk
        file_write: do
            call output_line(unit,magic,digest,io_status)
            if (io_status /= 0) exit file_write
            write(line,"(i0)") file_version
            call output_line(unit,trim(line),digest,io_status)
            if (io_status /= 0) exit file_write
            call output_line(unit,trim(configuration_id),digest,io_status)
            if (io_status /= 0) exit file_write
            write(line,"(i0,1x,i0,1x,i0)") digits(0.0_wp),range(0.0_wp),storage_size(0.0_wp)
            call output_line(unit,trim(line),digest,io_status)
            if (io_status /= 0) exit file_write
            write(line,"(*(i0,1x))") material%point_count(),material%input_size(),material%output_size(), material%state_size(),material%integer_state_size(),material%logical_state_size(), material%sensitivity_size(),material%previous_input_size()
            call output_line(unit,trim(line),digest,io_status)
            if (io_status /= 0) exit file_write
            call write_selected_parameters(unit,material,digest,io_status)
            if (io_status /= 0) exit file_write
            call write_snapshot(unit,snapshot,digest,io_status)
            if (io_status /= 0) exit file_write
            call output_line(unit,end_marker,digest,io_status)
            if (io_status /= 0) exit file_write
            write(unit,"(i0)",iostat=io_status,iomsg=io_message) digest
            exit file_write
        end do file_write

        if (io_status == 0) status = cf_ok
        close(unit,iostat=close_status,iomsg=io_message)
        if (close_status /= 0) status = cf_io_error
    end subroutine save_checkpoint

    subroutine load_checkpoint(material, path, configuration_id, status)
        !! Validate a version-1 checkpoint and restore accepted history.
        !!
        !! Reinitialize an identically configured material_t (including
        !! selected sensitivity parameter names in the same order) before
        !! loading. Configuration identity, precision and dimensions must
        !! match the saved file exactly. The point ordering is assumed to
        !! match the required application-owned configuration identity.
        !!
        !! A successful restore resets pending trials. An unsuccessful load
        !! cannot modify any part of the initialized material handle.
        !! cf_checkpoint_mismatch indicates wrong version/identity/layout.
        !! Cf_checkpoint_corrupt indicates damaged data, nonfinite history,
        !! checksum failure or trailing content. cf_io_error indicates an
        !! open/close failure. Other statuses indicate misuse or allocation.
        type(material_t), intent(inout) :: material
            !! Material handle whose configuration and history are queried or updated.
        character(*), intent(in) :: path
            !! Filesystem path of the checkpoint file.
        character(*), intent(in) :: configuration_id
            !! Stable application-defined identity for model parameters, units and point ordering.
        integer, intent(out) :: status
            !! Operation status code; `cf_ok` indicates success.
        type(snapshot_t) :: snapshot
        integer :: unit
        integer :: io_status
        integer :: close_status
        integer :: version
        integer :: digits_file
        integer :: range_file
        integer :: storage_file
        integer :: npoints
        integer :: nin
        integer :: nout
        integer :: nreal
        integer :: ninteger
        integer :: nlogical
        integer :: nselected
        integer :: nprevious
        integer :: parameter_id
        integer :: i
        integer(lk) :: digest
        integer(lk) :: expected_digest
        character(512) :: line
        character(256) :: io_message
        logical :: valid

        status = cf_invalid_argument
        if (len_trim(path) == 0 .or. .not. valid_identity(configuration_id)) return
        status = cf_uninitialized
        if (material%point_count() == 0) return

        call prepare_snapshot(material,snapshot,status)
        if (status /= cf_ok) return
        status = cf_io_error
        open(newunit=unit,file=path,status="old",action="read",form="formatted",iostat=io_status,iomsg=io_message)
        if (io_status /= 0) return
        status = cf_checkpoint_corrupt
        digest = 0_lk
        valid = .false.
        file_read: do
            call input_line(unit,line,digest,io_status)
            if (io_status /= 0) exit file_read
            if (trim(line) /= magic) exit file_read
            call input_line(unit,line,digest,io_status)
            if (io_status /= 0) exit file_read
            read(line,*,iostat=io_status,iomsg=io_message) version
            if (io_status /= 0) exit file_read
            if (version /= file_version) then
                status = cf_checkpoint_mismatch
                exit file_read
            end if
            call input_line(unit,line,digest,io_status)
            if (io_status /= 0) exit file_read
            if (trim(line) /= trim(configuration_id)) then
                status = cf_checkpoint_mismatch
                exit file_read
            end if
            call input_line(unit,line,digest,io_status)
            if (io_status /= 0) exit file_read
            read(line,*,iostat=io_status,iomsg=io_message) digits_file,range_file,storage_file
            if (io_status /= 0) exit file_read
            if (digits_file /= digits(0.0_wp) .or. range_file /= range(0.0_wp) .or. storage_file /= storage_size(0.0_wp)) then
                status = cf_checkpoint_mismatch
                exit file_read
            end if
            call input_line(unit,line,digest,io_status)
            if (io_status /= 0) exit file_read
            read(line,*,iostat=io_status,iomsg=io_message) npoints,nin,nout,nreal,ninteger,nlogical,nselected,nprevious
            if (io_status /= 0) exit file_read
            if (npoints /= material%point_count() .or. nin /= material%input_size() .or. nout /= material%output_size() .or. nreal /= material%state_size() .or. ninteger /= material%integer_state_size() .or. nlogical /= material%logical_state_size() .or. nselected /= material%sensitivity_size() .or. nprevious /= material%previous_input_size()) then
                status = cf_checkpoint_mismatch
                exit file_read
            end if
            selected_loop: do i = 1, material%sensitivity_size()
                call input_line(unit,line,digest,io_status)
                if (io_status /= 0) exit selected_loop
                read(line,*,iostat=io_status,iomsg=io_message) parameter_id
                if (io_status /= 0) exit selected_loop
                if (parameter_id /= material%sensitivity_parameter_id(i)) then
                    status = cf_checkpoint_mismatch
                    exit selected_loop
                end if
            end do selected_loop
            if (io_status /= 0 .or. status == cf_checkpoint_mismatch) exit file_read

            call read_snapshot(unit,snapshot,digest,io_status)
            if (io_status /= 0) exit file_read
            call input_line(unit,line,digest,io_status)
            if (io_status /= 0) exit file_read
            if (trim(line) /= end_marker) exit file_read
            read(unit,"(a)",iostat=io_status,iomsg=io_message) line
            if (io_status /= 0) exit file_read
            read(line,*,iostat=io_status,iomsg=io_message) expected_digest
            if (io_status /= 0) exit file_read
            if (digest /= expected_digest) exit file_read
            read(unit,"(a)",iostat=io_status,iomsg=io_message) line
            if (io_status /= iostat_end) exit file_read
            valid = .true.
            exit file_read
        end do file_read

        close(unit,iostat=close_status,iomsg=io_message)
        if (close_status /= 0) then
            status = cf_io_error
            return
        end if
        if (.not. valid) return

        call material%import_committed(snapshot%state,status,previous_input=snapshot%previous, integer_state=snapshot%integers,logical_state=snapshot%flags,sensitivity_state=snapshot%sensitivities)
        if (status /= cf_ok) status = cf_checkpoint_corrupt
    end subroutine load_checkpoint

    subroutine prepare_snapshot(material, snapshot, status)
        !! Allocate only the dimensions of the *initialized* destination.
        !! Never use untrusted on-disk dimension values for an allocation.
        type(material_t), intent(in) :: material
            !! Material handle whose configuration and history are queried or updated.
        type(snapshot_t), intent(out) :: snapshot
            !! Complete copy of the committed material history being serialized or loaded.
        integer, intent(out) :: status
            !! Operation status code; `cf_ok` indicates success.
        integer :: allocation_status
        integer :: npoints
        integer :: previous_points
        character(256) :: allocation_message

        npoints = material%point_count()
        previous_points = 0
        if (material%previous_input_size() > 0) previous_points = npoints
        status = cf_allocation_failed

        allocate(snapshot%state(material%state_size(),npoints),stat=allocation_status,errmsg=allocation_message)
        if (allocation_status /= 0) return
        allocate(snapshot%previous(material%previous_input_size(),previous_points),stat=allocation_status,errmsg=allocation_message)
        if (allocation_status /= 0) return
        allocate(snapshot%integers(material%integer_state_size(),npoints),stat=allocation_status,errmsg=allocation_message)
        if (allocation_status /= 0) return
        allocate(snapshot%flags(material%logical_state_size(),npoints),stat=allocation_status,errmsg=allocation_message)
        if (allocation_status /= 0) return
        allocate(snapshot%sensitivities(material%state_size(),material%sensitivity_size(),npoints), stat=allocation_status,errmsg=allocation_message)
        if (allocation_status /= 0) return
        status = cf_ok
    end subroutine prepare_snapshot

    pure logical function valid_identity(identity) result(valid)
        !! Accept printable one-line ASCII identities, including spaces.
        character(*), intent(in) :: identity
            !! Application-defined configuration identifier to validate.
        integer :: i
        integer :: value

        valid = .false.
        if (len_trim(identity) < 1 .or. len_trim(identity) > max_identity_length) return
        do i = 1, len_trim(identity)
            value = iachar(identity(i:i))
            if (value < 32 .or. value > 126) return
        end do
        valid = .true.
    end function valid_identity

    pure integer(lk) function update_digest(previous, line) result(value)
        !! Rotation/XOR is defined bitwise, avoiding signed integer overflow.
        integer(lk), intent(in) :: previous
            !! Previously accumulated checksum or committed input, as appropriate.
        character(*), intent(in) :: line
            !! One formatted checkpoint record without its line terminator.
        integer :: i

        value = previous
        do i = 1, len_trim(line)
            value = ieor(ishftc(value,5),int(iachar(line(i:i)),lk))
        end do
        value = ieor(ishftc(value,5),10_lk)
    end function update_digest

    subroutine output_line(unit, line, digest, io_status)
        !! Write one checkpoint line and incorporate it into the integrity checksum.
        integer, intent(in) :: unit
            !! Open formatted Fortran file unit used for checkpoint input or output.
        character(*), intent(in) :: line
            !! One formatted checkpoint record without its line terminator.
        integer(lk), intent(inout) :: digest
            !! Running rotate/XOR checkpoint integrity checksum.
        integer, intent(out) :: io_status
            !! Fortran I/O status; zero means success.
        character(256) :: io_message

        write(unit,"(a)",iostat=io_status,iomsg=io_message) trim(line)
        if (io_status == 0) digest = update_digest(digest,trim(line))
    end subroutine output_line

    subroutine input_line(unit, line, digest, io_status)
        !! Read one checkpoint line and update the integrity checksum.
        integer, intent(in) :: unit
            !! Open formatted Fortran file unit used for checkpoint input or output.
        character(512), intent(out) :: line
            !! One formatted checkpoint record without its line terminator.
        integer(lk), intent(inout) :: digest
            !! Running rotate/XOR checkpoint integrity checksum.
        integer, intent(out) :: io_status
            !! Fortran I/O status; zero means success.
        character(256) :: io_message

        read(unit,"(a)",iostat=io_status,iomsg=io_message) line
        if (io_status == 0) digest = update_digest(digest,trim(line))
    end subroutine input_line

    subroutine write_selected_parameters(unit, material, digest, io_status)
        !! Serialize the selected scalar parameter IDs in configuration order.
        integer, intent(in) :: unit
            !! Open formatted Fortran file unit used for checkpoint input or output.
        type(material_t), intent(in) :: material
            !! Material handle whose configuration and history are queried or updated.
        integer(lk), intent(inout) :: digest
            !! Running rotate/XOR checkpoint integrity checksum.
        integer, intent(out) :: io_status
            !! Fortran I/O status; zero means success.
        character(64) :: line
        integer :: i

        io_status = 0
        do i = 1, material%sensitivity_size()
            write(line,"(i0)") material%sensitivity_parameter_id(i)
            call output_line(unit,trim(line),digest,io_status)
            if (io_status /= 0) return
        end do
    end subroutine write_selected_parameters

    subroutine write_real(unit, number, digest, io_status)
        !! Write a real-valued checkpoint record with round-trip precision.
        integer, intent(in) :: unit
            !! Open formatted Fortran file unit used for checkpoint input or output.
        real(wp), intent(in) :: number
            !! Scalar value to serialize or parse in the checkpoint.
        integer(lk), intent(inout) :: digest
            !! Running rotate/XOR checkpoint integrity checksum.
        integer, intent(out) :: io_status
            !! Fortran I/O status; zero means success.
        character(64) :: line

        write(line,"(es50.38e4)") number
        call output_line(unit,trim(adjustl(line)),digest,io_status)
    end subroutine write_real

    subroutine write_integer(unit, number, digest, io_status)
        !! Write an integer-valued checkpoint record.
        integer, intent(in) :: unit
            !! Open formatted Fortran file unit used for checkpoint input or output.
        integer(ik), intent(in) :: number
            !! Scalar value to serialize or parse in the checkpoint.
        integer(lk), intent(inout) :: digest
            !! Running rotate/XOR checkpoint integrity checksum.
        integer, intent(out) :: io_status
            !! Fortran I/O status; zero means success.
        character(64) :: line

        write(line,"(i0)") number
        call output_line(unit,trim(line),digest,io_status)
    end subroutine write_integer

    subroutine write_snapshot(unit, snapshot, digest, io_status)
        !! Write in stable component-fastest, point-major order.
        integer, intent(in) :: unit
            !! Open formatted Fortran file unit used for checkpoint input or output.
        type(snapshot_t), intent(in) :: snapshot
            !! Complete copy of the committed material history being serialized or loaded.
        integer(lk), intent(inout) :: digest
            !! Running rotate/XOR checkpoint integrity checksum.
        integer, intent(out) :: io_status
            !! Fortran I/O status; zero means success.
        integer :: i
        integer :: j
        integer :: k

        io_status = 0
        do j = 1, size(snapshot%state,2)
            do i = 1, size(snapshot%state,1)
                call write_real(unit,snapshot%state(i,j),digest,io_status)
                if (io_status /= 0) return
            end do
        end do
        do j = 1, size(snapshot%previous,2)
            do i = 1, size(snapshot%previous,1)
                call write_real(unit,snapshot%previous(i,j),digest,io_status)
                if (io_status /= 0) return
            end do
        end do
        do j = 1, size(snapshot%integers,2)
            do i = 1, size(snapshot%integers,1)
                call write_integer(unit,snapshot%integers(i,j),digest,io_status)
                if (io_status /= 0) return
            end do
        end do
        do j = 1, size(snapshot%flags,2)
            do i = 1, size(snapshot%flags,1)
                call write_integer(unit,merge(1_ik,0_ik,snapshot%flags(i,j)),digest,io_status)
                if (io_status /= 0) return
            end do
        end do
        do k = 1, size(snapshot%sensitivities,3)
            do j = 1, size(snapshot%sensitivities,2)
                do i = 1, size(snapshot%sensitivities,1)
                    call write_real(unit,snapshot%sensitivities(i,j,k),digest,io_status)
                    if (io_status /= 0) return
                end do
            end do
        end do
    end subroutine write_snapshot

    subroutine read_real(unit, number, digest, io_status)
        !! Parse a real-valued checkpoint record and reject invalid values.
        integer, intent(in) :: unit
            !! Open formatted Fortran file unit used for checkpoint input or output.
        real(wp), intent(out) :: number
            !! Scalar value to serialize or parse in the checkpoint.
        integer(lk), intent(inout) :: digest
            !! Running rotate/XOR checkpoint integrity checksum.
        integer, intent(out) :: io_status
            !! Fortran I/O status; zero means success.
        character(512) :: line
        character(256) :: io_message

        call input_line(unit,line,digest,io_status)
        if (io_status /= 0) return
        read(line,*,iostat=io_status,iomsg=io_message) number
        if (io_status /= 0) return
        if (.not. ieee_is_finite(number)) io_status = 1
    end subroutine read_real

    subroutine read_integer(unit, number, digest, io_status)
        !! Parse an integer-valued checkpoint record.
        integer, intent(in) :: unit
            !! Open formatted Fortran file unit used for checkpoint input or output.
        integer(ik), intent(out) :: number
            !! Scalar value to serialize or parse in the checkpoint.
        integer(lk), intent(inout) :: digest
            !! Running rotate/XOR checkpoint integrity checksum.
        integer, intent(out) :: io_status
            !! Fortran I/O status; zero means success.
        character(512) :: line
        character(256) :: io_message

        call input_line(unit,line,digest,io_status)
        if (io_status /= 0) return
        read(line,*,iostat=io_status,iomsg=io_message) number
        if (io_status /= 0) return
    end subroutine read_integer

    subroutine read_snapshot(unit, snapshot, digest, io_status)
        !! Read the full checkpoint state and validate its saved data.
        integer, intent(in) :: unit
            !! Open formatted Fortran file unit used for checkpoint input or output.
        type(snapshot_t), intent(inout) :: snapshot
            !! Complete copy of the committed material history being serialized or loaded.
        integer(lk), intent(inout) :: digest
            !! Running rotate/XOR checkpoint integrity checksum.
        integer, intent(out) :: io_status
            !! Fortran I/O status; zero means success.
        integer(ik) :: marker
        integer :: i
        integer :: j
        integer :: k

        io_status = 0
        do j = 1, size(snapshot%state,2)
            do i = 1, size(snapshot%state,1)
                call read_real(unit,snapshot%state(i,j),digest,io_status)
                if (io_status /= 0) return
            end do
        end do
        do j = 1, size(snapshot%previous,2)
            do i = 1, size(snapshot%previous,1)
                call read_real(unit,snapshot%previous(i,j),digest,io_status)
                if (io_status /= 0) return
            end do
        end do
        do j = 1, size(snapshot%integers,2)
            do i = 1, size(snapshot%integers,1)
                call read_integer(unit,snapshot%integers(i,j),digest,io_status)
                if (io_status /= 0) return
            end do
        end do
        do j = 1, size(snapshot%flags,2)
            do i = 1, size(snapshot%flags,1)
                call read_integer(unit,marker,digest,io_status)
                if (io_status /= 0) return
                if (marker /= 0_ik .and. marker /= 1_ik) then
                    io_status = 1
                    return
                end if
                snapshot%flags(i,j) = marker == 1_ik
            end do
        end do
        do k = 1, size(snapshot%sensitivities,3)
            do j = 1, size(snapshot%sensitivities,2)
                do i = 1, size(snapshot%sensitivities,1)
                    call read_real(unit,snapshot%sensitivities(i,j,k),digest,io_status)
                    if (io_status /= 0) return
                end do
            end do
        end do
    end subroutine read_snapshot

end module constituf_checkpoint

! SPDX-License-Identifier: Apache-2.0

module constituf_status
    !! Status values used by ConstituF public operations and material kernels.
    implicit none (type, external)
    private

    integer, parameter, public :: cf_ok = 0
        !! Operation completed successfully.
    integer, parameter, public :: cf_invalid_argument = 1
        !! Invalid model parameters, field specification, or point count.
    integer, parameter, public :: cf_invalid_shape = 2
        !! A supplied array has the wrong extent or rank.
    integer, parameter, public :: cf_uninitialized = 3
        !! A material handle has not been initialized.
    integer, parameter, public :: cf_no_trial = 4
        !! There is no successful trial to commit.
    integer, parameter, public :: cf_invalid_time = 5
        !! A time-dependent model received an invalid time increment.
    integer, parameter, public :: cf_unsupported = 6
        !! The selected material does not support the requested operation.
    integer, parameter, public :: cf_evaluation_failed = 7
        !! At least one integration point reported failure.
    integer, parameter, public :: cf_allocation_failed = 8
        !! Material setup could not allocate the requested storage.
    integer, parameter, public :: cf_nondifferentiable = 9
        !! A requested classical derivative is undefined at a branch switch.
    integer, parameter, public :: cf_checkpoint_mismatch = 10
        !! Checkpoint precision, model signature, shape, or sensitivity order differs.
    integer, parameter, public :: cf_checkpoint_corrupt = 11
        !! Checkpoint data are malformed, truncated, nonfinite, or fail integrity.
    integer, parameter, public :: cf_io_error = 12
        !! Unable to open, write, read or close the requested checkpoint file.
end module constituf_status

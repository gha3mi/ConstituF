! SPDX-License-Identifier: Apache-2.0

program test_core_kinds
    !! Check the public kind contract for every supported build configuration.
    use, intrinsic :: iso_fortran_env, only: real32, real64, real128, int32, int64
    use constituf, only: wp, ik, lk
    implicit none (type, external)

#ifdef CONSTITUF_REAL32
    if (wp /= real32) error stop "Unexpected real32 kind"
#elif defined(CONSTITUF_REAL128)
    if (wp /= real128) error stop "Unexpected real128 kind"
#else
    if (wp /= real64) error stop "Unexpected default real64 kind"
#endif
    if (ik /= int32) error stop "Unexpected integer kind"
    if (lk /= int64) error stop "Unexpected large integer kind"
end program test_core_kinds

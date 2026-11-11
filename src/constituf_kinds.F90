! SPDX-License-Identifier: Apache-2.0

module constituf_kinds
    !! Central numerical kinds, selected once at build time.
    use, intrinsic :: iso_fortran_env, only: real32, real64, real128, int32, int64
    implicit none (type, external)
    private

    !! Working real kind. The default is real64. Select a different kind by
    !! enabling the real32 or real128 fpm feature. Rebuild all consumers when
    !! changing precision. Compiled modules of different kinds are incompatible.
#ifdef CONSTITUF_REAL32
#ifdef CONSTITUF_REAL128
#error "CONSTITUF_REAL32 and CONSTITUF_REAL128 are mutually exclusive"
#endif
    integer, parameter, public :: wp = real32
#elif defined(CONSTITUF_REAL128)
    integer, parameter, public :: wp = real128
#else
    integer, parameter, public :: wp = real64
#endif
    integer, parameter, public :: ik = int32
        !! Kind for integer-valued material history and explicit integer data.
    integer, parameter, public :: lk = int64
        !! Kind for large counts and offsets.
end module constituf_kinds

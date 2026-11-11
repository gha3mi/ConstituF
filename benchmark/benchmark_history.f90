! SPDX-License-Identifier: Apache-2.0

module benchmark_history
    !! Accepted dimensionless loading amplitudes shared by both geometries.
    use constituf, only: wp
    implicit none (type, external)
    private

    integer, parameter, public :: npath=5
    real(wp), parameter, public :: path(npath)=[0.25_wp,0.65_wp,1.0_wp,0.40_wp,0.85_wp]
end module benchmark_history

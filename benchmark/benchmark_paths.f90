! SPDX-License-Identifier: Apache-2.0

module benchmark_paths
    !! Single location and naming convention for all benchmark exports.
    !! Kept within FORD media so the documentation uses the same artifacts.
    implicit none (type, external)
    private

    character(*), parameter, public :: benchmark_directory = "docs/assets/benchmarks/"
    character(*), parameter, public :: annulus_directory = benchmark_directory//"annulus/"
    character(*), parameter, public :: cavity3d_directory = benchmark_directory//"cavity3d/"
end module benchmark_paths

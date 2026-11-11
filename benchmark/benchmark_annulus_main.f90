! SPDX-License-Identifier: Apache-2.0

program benchmark_annulus_main
    !! Prescribed two-dimensional annulus benchmark. Optional mode and slug:
    !! benchmark_annulus [all|response|sensitivity] [material_slug]
    use benchmark_annulus_runner, only: run_annulus
    implicit none (type, external)

    character(32) :: mode
    character(64) :: material
    integer :: n

    n=command_argument_count()
    if (n>2) error stop "Usage: benchmark_annulus [all|response|sensitivity] [material_slug]"
    mode="all"
    if (n>=1) call get_command_argument(1,mode)
    if (n==2) then
        call get_command_argument(2,material)
        if (len_trim(material)==0) error stop "Empty annulus material name"
        call run_annulus(trim(mode),trim(material))
    else
        call run_annulus(trim(mode))
    end if
end program benchmark_annulus_main

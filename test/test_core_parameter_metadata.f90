! SPDX-License-Identifier: Apache-2.0

program test_core_parameter_metadata
    !! Ordered scalar metadata is part of the public setup API, independent
    !! of whether analytical sensitivities were selected at initialization.
    !! Test data use conductivity in W/(m*K), elastic moduli in Pa and time in s.
    use constituf, only: wp, material_t, cf_ok
    use constituf_fourier_isotropic, only: fourier_isotropic_t
    use constituf_maxwell_multi, only: maxwell_multi_t
    implicit none (type, external)

    integer, parameter :: nmandel=6, nparameters=6

    type(material_t) :: material
    character(:), allocatable :: name
    real(wp) :: initial(nmandel,1)
    integer :: status,j
    logical :: found
    character(26), parameter :: expected(nparameters)=[character(len=26) :: &
        "bulk_modulus","equilibrium_shear_modulus", &
        "branch_shear_modulus_1","relaxation_time_1", &
        "branch_shear_modulus_2","relaxation_time_2"]

    if (material%parameter_count()/=0) error stop "Uninitialized material has parameters"
    call material%parameter_info(1,name,found)
    if (found.or.len(name)>0) error stop "Uninitialized parameter lookup"

    call material%init(fourier_isotropic_t(25.0_wp),1,status)
    if (status/=cf_ok) error stop "Fourier metadata setup"
    if (material%parameter_count()/=1) error stop "Fourier parameter count"
    call material%parameter_info(1,name,found)
    if (.not.found.or.name/="conductivity") error stop "Fourier parameter info"
    call material%parameter_info(2,name,found)
    if (found.or.len(name)>0) error stop "Invalid parameter info"
    call material%parameter_info(0,name,found)
    if (found.or.len(name)>0) error stop "Zero parameter index"

    initial=0.0_wp
    call material%init(maxwell_multi_t(80.0e6_wp,15.0e6_wp, &
        [30.0e6_wp,10.0e6_wp],[0.4_wp,2.0_wp]),1,status,initial_input=initial)
    if (status/=cf_ok) error stop "Maxwell metadata setup"
    if (material%parameter_count()/=size(expected)) error stop "Maxwell parameter count"
    do j=1,size(expected)
        call material%parameter_info(j,name,found)
        if (.not.found.or.name/=trim(expected(j))) error stop "Maxwell parameter order"
    end do
end program test_core_parameter_metadata

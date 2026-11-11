! SPDX-License-Identifier: Apache-2.0

module benchmark_models
    !! Single source for benchmark model descriptors, reference parameters,
    !! SI unit labels and native ConstituF constructors. Unit choices belong
    !! to these benchmark data. No material law or sensitivity is
    !! reimplemented here.
    use constituf, only: wp, material_t, cf_ok
    use constituf_axisymmetric_elastic, only: axisymmetric_elastic_t
    use constituf_axisymmetric_j2, only: axisymmetric_j2_t
    use constituf_biot_poroelastic, only: biot_poroelastic_t
    use constituf_chaboche, only: chaboche_t
    use constituf_cohesive_mixed_mode, only: cohesive_mixed_mode_t
    use constituf_cohesive_mode_i, only: cohesive_mode_i_t
    use constituf_darcy_anisotropic, only: darcy_anisotropic_t
    use constituf_darcy_exponential, only: darcy_exponential_t
    use constituf_darcy_forchheimer, only: darcy_forchheimer_t
    use constituf_darcy_isotropic, only: darcy_isotropic_t
    use constituf_dielectric_anisotropic, only: dielectric_anisotropic_t
    use constituf_dielectric_cubic, only: dielectric_cubic_t
    use constituf_dielectric_isotropic, only: dielectric_isotropic_t
    use constituf_drucker_prager, only: drucker_prager_t
    use constituf_ductile_damage, only: ductile_damage_t
    use constituf_fiber_hgo, only: fiber_hgo_t
    use constituf_fick_anisotropic, only: fick_anisotropic_t
    use constituf_fick_isotropic, only: fick_isotropic_t
    use constituf_fourier_anisotropic, only: fourier_anisotropic_t
    use constituf_fourier_exponential, only: fourier_exponential_t
    use constituf_fourier_isotropic, only: fourier_isotropic_t
    use constituf_gent, only: gent_t
    use constituf_hill_plasticity, only: hill_plasticity_t
    use constituf_isotropic_elastic, only: isotropic_elastic_t
    use constituf_isotropic_thermoelastic, only: isotropic_thermoelastic_t
    use constituf_j2_plane, only: j2_plane_t
    use constituf_j2_plasticity, only: j2_plasticity_t
    use constituf_kelvin_voigt, only: kelvin_voigt_t
    use constituf_maxwell_multi, only: maxwell_multi_t
    use constituf_maxwell_shear, only: maxwell_shear_t
    use constituf_mooney_rivlin, only: mooney_rivlin_t
    use constituf_neo_hookean, only: neo_hookean_t
    use constituf_ogden, only: ogden_t
    use constituf_ohmic_anisotropic, only: ohmic_anisotropic_t
    use constituf_ohmic_isotropic, only: ohmic_isotropic_t
    use constituf_orthotropic_elastic, only: orthotropic_elastic_t
    use constituf_perzyna, only: perzyna_t
    use constituf_piezoelectric_linear, only: piezoelectric_linear_t
    use constituf_plane_elastic, only: plane_elastic_t
    use constituf_yeoh, only: yeoh_t
    implicit none (type, external)
    private

    integer, parameter, public :: max_parameters=10
    integer, parameter, public :: load_mandel=1, load_plane=2, load_axisym=3, load_F=4
    integer, parameter, public :: load_cohesive_i=5, load_cohesive_mixed=6, load_thermoelastic=7
    integer, parameter, public :: load_heat=8, load_heat_temperature=9, load_concentration=10
    integer, parameter, public :: load_electric=11, load_piezo=12, load_pressure=13
    integer, parameter, public :: load_pressure_pressure=14, load_biot=15
    integer, parameter :: palette_vik=1
    integer, parameter :: palette_magma=2
    integer, parameter :: palette_viridis=3
    integer, parameter :: palette_batlow=4

    type, public :: benchmark_case_spec_t
        character(32) :: name
        integer :: response_component
        character(32) :: response_label
        character(24) :: response_units
        logical :: has_history
        integer :: input_kind
        integer :: colormap_kind
        integer :: parameter_count
        character(24) :: parameter_units(max_parameters)
        real(wp) :: reference_values(max_parameters)
    end type benchmark_case_spec_t

    ! Register a case here and add its native constructor below.
    ! SI parameter labels and reference values share the model's declared order.
    type(benchmark_case_spec_t), parameter, public :: cases(*) = [ &
        benchmark_case_spec_t( &
            name="isotropic_elastic", response_component=1, response_label="sigma_rr", &
            response_units="Pa", has_history=.false., &
            input_kind=load_mandel, &
            colormap_kind=palette_vik, parameter_count=2, &
            parameter_units=[character(len=24) :: "Pa", "Pa", "", "", "", "", "", "", "", ""], &
            reference_values=[ &
                150.0e9_wp, 75.0e9_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="plane_elastic", response_component=2, response_label="sigma_rr", &
            response_units="Pa", has_history=.false., &
            input_kind=load_plane, &
            colormap_kind=palette_vik, parameter_count=2, &
            parameter_units=[character(len=24) :: "Pa", "Pa", "", "", "", "", "", "", "", ""], &
            reference_values=[ &
                150.0e9_wp, 75.0e9_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="axisymmetric_elastic", response_component=3, response_label="sigma_rr", &
            response_units="Pa", has_history=.false., &
            input_kind=load_axisym, &
            colormap_kind=palette_vik, parameter_count=2, &
            parameter_units=[character(len=24) :: "Pa", "Pa", "", "", "", "", "", "", "", ""], &
            reference_values=[ &
                150.0e9_wp, 75.0e9_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="orthotropic_elastic", response_component=1, response_label="sigma_rr", &
            response_units="Pa", has_history=.false., &
            input_kind=load_mandel, &
            colormap_kind=palette_vik, parameter_count=9, &
            parameter_units=[character(len=24) :: "Pa", "Pa", "Pa", "1", "1", "1", "Pa", "Pa", "Pa", ""], &
            reference_values=[ &
                150.0e9_wp, 100.0e9_wp, 75.0e9_wp, 0.25_wp, 0.18_wp, 0.15_wp, 35.0e9_wp, 40.0e9_wp, 45.0e9_wp, &
                0.0_wp]), &
        benchmark_case_spec_t( &
            name="neo_hookean", response_component=4, response_label="P_rR", &
            response_units="Pa", has_history=.false., &
            input_kind=load_F, &
            colormap_kind=palette_vik, parameter_count=2, &
            parameter_units=[character(len=24) :: "Pa", "Pa", "", "", "", "", "", "", "", ""], &
            reference_values=[ &
                8.0e6_wp, 4.0e6_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="mooney_rivlin", response_component=4, response_label="P_rR", &
            response_units="Pa", has_history=.false., &
            input_kind=load_F, &
            colormap_kind=palette_vik, parameter_count=3, &
            parameter_units=[character(len=24) :: "Pa", "Pa", "Pa", "", "", "", "", "", "", ""], &
            reference_values=[ &
                1.5e6_wp, 0.5e6_wp, 1.0e8_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="gent", response_component=4, response_label="P_rR", &
            response_units="Pa", has_history=.false., &
            input_kind=load_F, &
            colormap_kind=palette_vik, parameter_count=3, &
            parameter_units=[character(len=24) :: "Pa", "1", "Pa", "", "", "", "", "", "", ""], &
            reference_values=[ &
                1.5e6_wp, 8.0_wp, 1.0e8_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="fiber_hgo", response_component=4, response_label="P_rR", &
            response_units="Pa", has_history=.false., &
            input_kind=load_F, &
            colormap_kind=palette_vik, parameter_count=4, &
            parameter_units=[character(len=24) :: "Pa", "Pa", "Pa", "1", "", "", "", "", "", ""], &
            reference_values=[ &
                6.0e6_wp, 8.0e6_wp, 4.0e6_wp, 3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="yeoh", response_component=4, response_label="P_rR", &
            response_units="Pa", has_history=.false., &
            input_kind=load_F, &
            colormap_kind=palette_vik, parameter_count=4, &
            parameter_units=[character(len=24) :: "Pa", "Pa", "Pa", "Pa", "", "", "", "", "", ""], &
            reference_values=[ &
                1.5e6_wp, 0.3e6_wp, 0.05e6_wp, 1.0e8_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="ogden", response_component=4, response_label="P_rR", &
            response_units="Pa", has_history=.false., &
            input_kind=load_F, &
            colormap_kind=palette_vik, parameter_count=3, &
            parameter_units=[character(len=24) :: "Pa", "1", "Pa", "", "", "", "", "", "", ""], &
            reference_values=[ &
                1.5e6_wp, 1.4_wp, 1.0e8_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="j2_plasticity", response_component=1, response_label="sigma_rr", &
            response_units="Pa", has_history=.true., &
            input_kind=load_mandel, &
            colormap_kind=palette_vik, parameter_count=4, &
            parameter_units=[character(len=24) :: "Pa", "Pa", "Pa", "Pa", "", "", "", "", "", ""], &
            reference_values=[ &
                150.0e9_wp, 75.0e9_wp, 350.0e6_wp, 2.0e9_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="ductile_damage", response_component=1, response_label="sigma_rr", &
            response_units="Pa", has_history=.true., &
            input_kind=load_mandel, &
            colormap_kind=palette_vik, parameter_count=7, &
            parameter_units=[character(len=24) :: "Pa", "Pa", "Pa", "Pa", "1", "1", "1", "", "", ""], &
            reference_values=[ &
                150.0e9_wp, 75.0e9_wp, 350.0e6_wp, 2.0e9_wp, 0.001_wp, 150.0_wp, 0.85_wp, 0.0_wp, 0.0_wp, &
                0.0_wp]), &
        benchmark_case_spec_t( &
            name="perzyna", response_component=1, response_label="sigma_rr", &
            response_units="Pa", has_history=.true., &
            input_kind=load_mandel, &
            colormap_kind=palette_vik, parameter_count=5, &
            parameter_units=[character(len=24) :: "Pa", "Pa", "Pa", "Pa", "Pa.s", "", "", "", "", ""], &
            reference_values=[ &
                160.0e9_wp, 80.0e9_wp, 350.0e6_wp, 2.0e9_wp, 50.0e6_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="j2_plane", response_component=2, response_label="sigma_rr", &
            response_units="Pa", has_history=.true., &
            input_kind=load_plane, &
            colormap_kind=palette_vik, parameter_count=4, &
            parameter_units=[character(len=24) :: "Pa", "Pa", "Pa", "Pa", "", "", "", "", "", ""], &
            reference_values=[ &
                150.0e9_wp, 75.0e9_wp, 350.0e6_wp, 2.0e9_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="axisymmetric_j2", response_component=3, response_label="sigma_rr", &
            response_units="Pa", has_history=.true., &
            input_kind=load_axisym, &
            colormap_kind=palette_vik, parameter_count=4, &
            parameter_units=[character(len=24) :: "Pa", "Pa", "Pa", "Pa", "", "", "", "", "", ""], &
            reference_values=[ &
                150.0e9_wp, 75.0e9_wp, 350.0e6_wp, 2.0e9_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="hill_plasticity", response_component=1, response_label="sigma_rr", &
            response_units="Pa", has_history=.true., &
            input_kind=load_mandel, &
            colormap_kind=palette_vik, parameter_count=10, &
            parameter_units=[character(len=24) :: "Pa", "Pa", "Pa", "Pa", "1", "1", "1", "1", "1", "1"], &
            reference_values=[ &
                150.0e9_wp, 75.0e9_wp, 350.0e6_wp, 2.0e9_wp, 0.6_wp, 0.4_wp, 0.5_wp, 1.6_wp, 1.3_wp, 1.8_wp]), &
        benchmark_case_spec_t( &
            name="chaboche", response_component=1, response_label="sigma_rr", &
            response_units="Pa", has_history=.true., &
            input_kind=load_mandel, &
            colormap_kind=palette_vik, parameter_count=6, &
            parameter_units=[character(len=24) :: "Pa", "Pa", "Pa", "Pa", "Pa", "1", "", "", "", ""], &
            reference_values=[ &
                100.0e9_wp, 40.0e9_wp, 180.0e6_wp, 1.0e9_wp, 9.0e9_wp, 12.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="drucker_prager", response_component=1, response_label="sigma_rr", &
            response_units="Pa", has_history=.true., &
            input_kind=load_mandel, &
            colormap_kind=palette_vik, parameter_count=5, &
            parameter_units=[character(len=24) :: "Pa", "Pa", "Pa", "1", "Pa", "", "", "", "", ""], &
            reference_values=[ &
                3.0e9_wp, 1.0e9_wp, 10.0e6_wp, 0.15_wp, 0.1e9_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="kelvin_voigt", response_component=1, response_label="sigma_rr", &
            response_units="Pa", has_history=.true., &
            input_kind=load_mandel, &
            colormap_kind=palette_vik, parameter_count=4, &
            parameter_units=[character(len=24) :: "Pa", "Pa", "Pa*s", "Pa*s", "", "", "", "", "", ""], &
            reference_values=[ &
                100.0e6_wp, 40.0e6_wp, 15.0e6_wp, 8.0e6_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="maxwell_multi", response_component=1, response_label="sigma_rr", &
            response_units="Pa", has_history=.true., &
            input_kind=load_mandel, &
            colormap_kind=palette_vik, parameter_count=6, &
            parameter_units=[character(len=24) :: "Pa", "Pa", "Pa", "s", "Pa", "s", "", "", "", ""], &
            reference_values=[ &
                80.0e6_wp, 15.0e6_wp, 30.0e6_wp, 0.4_wp, 10.0e6_wp, 2.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="maxwell_shear", response_component=1, response_label="sigma_rr", &
            response_units="Pa", has_history=.true., &
            input_kind=load_mandel, &
            colormap_kind=palette_vik, parameter_count=4, &
            parameter_units=[character(len=24) :: "Pa", "Pa", "Pa", "s", "", "", "", "", "", ""], &
            reference_values=[ &
                160.0e6_wp, 40.0e6_wp, 30.0e6_wp, 2.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="cohesive_mode_i", response_component=6, response_label="traction_normal", &
            response_units="Pa", has_history=.true., &
            input_kind=load_cohesive_i, &
            colormap_kind=palette_magma, parameter_count=3, &
            parameter_units=[character(len=24) :: "Pa/m", "m", "m", "", "", "", "", "", "", ""], &
            reference_values=[ &
                2.0e12_wp, 1.0e-5_wp, 5.0e-5_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="cohesive_mixed_mode", response_component=6, response_label="traction_normal", &
            response_units="Pa", has_history=.true., &
            input_kind=load_cohesive_mixed, &
            colormap_kind=palette_magma, parameter_count=4, &
            parameter_units=[character(len=24) :: "Pa/m", "m", "m", "1", "", "", "", "", "", ""], &
            reference_values=[ &
                2.0e12_wp, 1.0e-5_wp, 5.0e-5_wp, 1.3_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="isotropic_thermoelastic", response_component=1, response_label="sigma_rr", &
            response_units="Pa", has_history=.false., &
            input_kind=load_thermoelastic, &
            colormap_kind=palette_vik, parameter_count=4, &
            parameter_units=[character(len=24) :: "Pa", "Pa", "1/K", "K", "", "", "", "", "", ""], &
            reference_values=[ &
                100.0e9_wp, 40.0e9_wp, 12.0e-6_wp, 293.15_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="fourier_isotropic", response_component=5, response_label="heat_flux_radial", &
            response_units="W/m^2", has_history=.false., &
            input_kind=load_heat, &
            colormap_kind=palette_magma, parameter_count=1, &
            parameter_units=[character(len=24) :: "W/(m*K)", "", "", "", "", "", "", "", "", ""], &
            reference_values=[ &
                25.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="fourier_anisotropic", response_component=5, response_label="heat_flux_radial", &
            response_units="W/m^2", has_history=.false., &
            input_kind=load_heat, &
            colormap_kind=palette_magma, parameter_count=6, &
            parameter_units=[character(len=24) :: "W/(m*K)", "W/(m*K)", "W/(m*K)", "W/(m*K)", "W/(m*K)", "W/(m*K)", "", "", "", ""], &
            reference_values=[ &
                10.0_wp, 8.0_wp, 6.0_wp, 2.0_wp, -1.0_wp, 1.5_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="fourier_exponential", response_component=5, response_label="heat_flux_radial", &
            response_units="W/m^2", has_history=.false., &
            input_kind=load_heat_temperature, &
            colormap_kind=palette_magma, parameter_count=3, &
            parameter_units=[character(len=24) :: "W/(m*K)", "1/K", "K", "", "", "", "", "", "", ""], &
            reference_values=[ &
                10.0_wp, 0.01_wp, 300.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="fick_isotropic", response_component=5, response_label="molar_flux_radial", &
            response_units="mol/(m^2 s)", has_history=.false., &
            input_kind=load_concentration, &
            colormap_kind=palette_viridis, parameter_count=1, &
            parameter_units=[character(len=24) :: "m2/s", "", "", "", "", "", "", "", "", ""], &
            reference_values=[ &
                1.5e-9_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="fick_anisotropic", response_component=5, response_label="molar_flux_radial", &
            response_units="mol/(m^2 s)", has_history=.false., &
            input_kind=load_concentration, &
            colormap_kind=palette_viridis, parameter_count=6, &
            parameter_units=[character(len=24) :: "m2/s", "m2/s", "m2/s", "m2/s", "m2/s", "m2/s", "", "", "", ""], &
            reference_values=[ &
                3.0e-9_wp, 2.0e-9_wp, 1.0e-9_wp, 0.2e-9_wp, -0.1e-9_wp, 0.15e-9_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                0.0_wp]), &
        benchmark_case_spec_t( &
            name="ohmic_isotropic", response_component=5, response_label="current_density_radial", &
            response_units="A/m^2", has_history=.false., &
            input_kind=load_electric, &
            colormap_kind=palette_batlow, parameter_count=1, &
            parameter_units=[character(len=24) :: "S/m", "", "", "", "", "", "", "", "", ""], &
            reference_values=[ &
                3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="ohmic_anisotropic", response_component=5, response_label="current_density_radial", &
            response_units="A/m^2", has_history=.false., &
            input_kind=load_electric, &
            colormap_kind=palette_batlow, parameter_count=6, &
            parameter_units=[character(len=24) :: "S/m", "S/m", "S/m", "S/m", "S/m", "S/m", "", "", "", ""], &
            reference_values=[ &
                4.0_wp, 3.0_wp, 2.0_wp, 0.5_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="dielectric_isotropic", response_component=5, response_label="electric_displacement_radial", &
            response_units="C/m^2", has_history=.false., &
            input_kind=load_electric, &
            colormap_kind=palette_batlow, parameter_count=1, &
            parameter_units=[character(len=24) :: "F/m", "", "", "", "", "", "", "", "", ""], &
            reference_values=[ &
                3.0e-11_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="dielectric_anisotropic", response_component=5, response_label="electric_displacement_radial", &
            response_units="C/m^2", has_history=.false., &
            input_kind=load_electric, &
            colormap_kind=palette_batlow, parameter_count=6, &
            parameter_units=[character(len=24) :: "F/m", "F/m", "F/m", "F/m", "F/m", "F/m", "", "", "", ""], &
            reference_values=[ &
                4.0e-11_wp, 3.0e-11_wp, 2.0e-11_wp, 0.5e-11_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="dielectric_cubic", response_component=5, response_label="electric_displacement_radial", &
            response_units="C/m^2", has_history=.false., &
            input_kind=load_electric, &
            colormap_kind=palette_batlow, parameter_count=2, &
            parameter_units=[character(len=24) :: "F/m", "F*m/V2", "", "", "", "", "", "", "", ""], &
            reference_values=[ &
                2.0e-11_wp, 1.0e-22_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="piezoelectric_linear", response_component=1, response_label="sigma_rr", &
            response_units="Pa", has_history=.false., &
            input_kind=load_piezo, &
            colormap_kind=palette_vik, parameter_count=7, &
            parameter_units=[character(len=24) :: "Pa", "Pa", "C/m2", "C/m2", "C/m2", "F/m", "F/m", "", "", ""], &
            reference_values=[ &
                100.0e9_wp, 40.0e9_wp, -2.0_wp, 4.0_wp, 3.0_wp, 1.0e-8_wp, 1.5e-8_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="darcy_isotropic", response_component=5, response_label="darcy_velocity_radial", &
            response_units="m/s", has_history=.false., &
            input_kind=load_pressure, &
            colormap_kind=palette_viridis, parameter_count=2, &
            parameter_units=[character(len=24) :: "m2", "Pa*s", "", "", "", "", "", "", "", ""], &
            reference_values=[ &
                1.0e-12_wp, 1.0e-3_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="darcy_anisotropic", response_component=5, response_label="darcy_velocity_radial", &
            response_units="m/s", has_history=.false., &
            input_kind=load_pressure, &
            colormap_kind=palette_viridis, parameter_count=7, &
            parameter_units=[character(len=24) :: "m2", "m2", "m2", "m2", "m2", "m2", "Pa*s", "", "", ""], &
            reference_values=[ &
                4.0e-12_wp, 3.0e-12_wp, 2.0e-12_wp, 0.5e-12_wp, 0.0_wp, 0.0_wp, 1.0e-3_wp, 0.0_wp, 0.0_wp, &
                0.0_wp]), &
        benchmark_case_spec_t( &
            name="darcy_forchheimer", response_component=5, response_label="darcy_velocity_radial", &
            response_units="m/s", has_history=.false., &
            input_kind=load_pressure, &
            colormap_kind=palette_viridis, parameter_count=2, &
            parameter_units=[character(len=24) :: "Pa*s/m2", "Pa*s2/m3", "", "", "", "", "", "", "", ""], &
            reference_values=[ &
                2.0e8_wp, 4.0e9_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="darcy_exponential", response_component=5, response_label="darcy_velocity_radial", &
            response_units="m/s", has_history=.false., &
            input_kind=load_pressure_pressure, &
            colormap_kind=palette_viridis, parameter_count=4, &
            parameter_units=[character(len=24) :: "m2", "Pa*s", "1/Pa", "Pa", "", "", "", "", "", ""], &
            reference_values=[ &
                1.0e-12_wp, 1.0e-3_wp, 1.0e-6_wp, 1.0e6_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]), &
        benchmark_case_spec_t( &
            name="biot_poroelastic", response_component=1, response_label="sigma_rr", &
            response_units="Pa", has_history=.false., &
            input_kind=load_biot, &
            colormap_kind=palette_vik, parameter_count=4, &
            parameter_units=[character(len=24) :: "Pa", "Pa", "1", "Pa", "", "", "", "", "", ""], &
            reference_values=[ &
                5.0e9_wp, 2.0e9_wp, 0.85_wp, 1.0e10_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]) &
    ]
    integer, parameter, public :: ncases=size(cases)
    character(32), parameter, public :: names(ncases)=cases%name

    public :: initialize_model, reference_parameters, response_kind
    public :: interface_case, history_case
    public :: latex_response, latex_units, latex_number

contains

    subroutine initialize_model(id,material,npoint,selected_parameters)
        integer, intent(in) :: id,npoint
        type(material_t), intent(inout) :: material
        character(*), intent(in), optional :: selected_parameters(:)
        real(wp) :: p(max_parameters)
        real(wp), allocatable :: initial(:,:)
        integer :: status

        if (id<1.or.id>ncases) error stop "Unknown annulus model"
        p=cases(id)%reference_values
        select case(trim(cases(id)%name))
        case("isotropic_elastic")
            call material%init(isotropic_elastic_t(p(1), p(2)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("plane_elastic")
            call material%init(plane_elastic_t(p(1), p(2), .true.), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("axisymmetric_elastic")
            call material%init(axisymmetric_elastic_t(p(1), p(2)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("orthotropic_elastic")
            call material%init(orthotropic_elastic_t(p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8), p(9)), npoint, &
                              status, sensitivity_parameters=selected_parameters)
        case("neo_hookean")
            call material%init(neo_hookean_t(p(1), p(2)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("mooney_rivlin")
            call material%init(mooney_rivlin_t(p(1), p(2), p(3)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("gent")
            call material%init(gent_t(p(1), p(2), p(3)), npoint, status, sensitivity_parameters=selected_parameters)
        case("fiber_hgo")
            call material%init(fiber_hgo_t(p(1), p(2), p(3), p(4)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("yeoh")
            call material%init(yeoh_t(p(1), p(2), p(3), p(4)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("ogden")
            call material%init(ogden_t(p(1), p(2), p(3)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("j2_plasticity")
            call material%init(j2_plasticity_t(p(1), p(2), p(3), p(4)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("ductile_damage")
            call material%init(ductile_damage_t(p(1), p(2), p(3), p(4), p(5), p(6), p(7)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("perzyna")
            allocate(initial(6,npoint))
            initial=0.0_wp
            call material%init(perzyna_t(p(1), p(2), p(3), p(4), p(5), 0.0004_wp), npoint, status, &
                              initial_input=initial(1:6, :), sensitivity_parameters=selected_parameters)
        case("j2_plane")
            call material%init(j2_plane_t(p(1), p(2), p(3), p(4), .true.), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("axisymmetric_j2")
            call material%init(axisymmetric_j2_t(p(1), p(2), p(3), p(4)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("hill_plasticity")
            call material%init(hill_plasticity_t(p(1), p(2), p(3), p(4), [p(5), p(6), p(7), p(8), p(9), p(10)]), &
                              npoint, status, sensitivity_parameters=selected_parameters)
        case("chaboche")
            call material%init(chaboche_t(p(1), p(2), p(3), p(4), p(5), p(6)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("drucker_prager")
            call material%init(drucker_prager_t(p(1), p(2), p(3), p(4), p(5)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("kelvin_voigt")
            allocate(initial(6,npoint))
            initial=0.0_wp
            call material%init(kelvin_voigt_t(p(1), p(2), p(3), p(4)), npoint, status, initial_input=initial(1:6, &
                              :), sensitivity_parameters=selected_parameters)
        case("maxwell_multi")
            allocate(initial(6,npoint))
            initial=0.0_wp
            call material%init(maxwell_multi_t(p(1), p(2), [p(3), p(5)], [p(4), p(6)]), npoint, status, &
                              initial_input=initial(1:6, :), sensitivity_parameters=selected_parameters)
        case("maxwell_shear")
            allocate(initial(6,npoint))
            initial=0.0_wp
            call material%init(maxwell_shear_t(p(1), p(2), p(3), p(4)), npoint, status, initial_input=initial(1:6, &
                              :), sensitivity_parameters=selected_parameters)
        case("cohesive_mode_i")
            call material%init(cohesive_mode_i_t(p(1), p(2), p(3)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("cohesive_mixed_mode")
            call material%init(cohesive_mixed_mode_t(p(1), p(2), p(3), p(4)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("isotropic_thermoelastic")
            call material%init(isotropic_thermoelastic_t(p(1), p(2), p(3), p(4)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("fourier_isotropic")
            call material%init(fourier_isotropic_t(p(1)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("fourier_anisotropic")
            call material%init(fourier_anisotropic_t(p(1), p(2), p(3), p(4), p(5), p(6)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("fourier_exponential")
            call material%init(fourier_exponential_t(p(1), p(2), p(3)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("fick_isotropic")
            call material%init(fick_isotropic_t(p(1)), npoint, status, sensitivity_parameters=selected_parameters)
        case("fick_anisotropic")
            call material%init(fick_anisotropic_t(p(1), p(2), p(3), p(4), p(5), p(6)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("ohmic_isotropic")
            call material%init(ohmic_isotropic_t(p(1)), npoint, status, sensitivity_parameters=selected_parameters)
        case("ohmic_anisotropic")
            call material%init(ohmic_anisotropic_t(p(1), p(2), p(3), p(4), p(5), p(6)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("dielectric_isotropic")
            call material%init(dielectric_isotropic_t(p(1)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("dielectric_anisotropic")
            call material%init(dielectric_anisotropic_t(p(1), p(2), p(3), p(4), p(5), p(6)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("dielectric_cubic")
            call material%init(dielectric_cubic_t(p(1), p(2)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("piezoelectric_linear")
            call material%init(piezoelectric_linear_t(p(1), p(2), p(3), p(4), p(5), p(6), p(7)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("darcy_isotropic")
            call material%init(darcy_isotropic_t(p(1), p(2)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("darcy_anisotropic")
            call material%init(darcy_anisotropic_t(p(1), p(2), p(3), p(4), p(5), p(6), p(7)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("darcy_forchheimer")
            call material%init(darcy_forchheimer_t(p(1), p(2)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("darcy_exponential")
            call material%init(darcy_exponential_t(p(1), p(2), p(3), p(4)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case("biot_poroelastic")
            call material%init(biot_poroelastic_t(p(1), p(2), p(3), p(4)), npoint, status, &
                              sensitivity_parameters=selected_parameters)
        case default
            error stop "Annulus model has no ConstituF constructor"
        end select
        if (status/=cf_ok) then
            print "(a,1x,a,1x,i0)", "Annulus model initialization:",trim(names(id)),status
            error stop 1
        end if
    end subroutine initialize_model

    subroutine reference_parameters(id,values,n)
        integer, intent(in) :: id
        real(wp), intent(out) :: values(:)
        integer, intent(out) :: n

        if (id<1.or.id>ncases) error stop "Unknown annulus reference parameters"
        n=cases(id)%parameter_count
        if (size(values)<n) error stop "Reference parameter array too small"
        values=0.0_wp
        values(:n)=cases(id)%reference_values(:n)
    end subroutine reference_parameters

    pure logical function interface_case(id) result(value)
        integer, intent(in) :: id
        value=.false.
        if (id>=1.and.id<=ncases) then
            value=cases(id)%input_kind==load_cohesive_i.or.cases(id)%input_kind==load_cohesive_mixed
        end if
    end function interface_case

    pure logical function history_case(id) result(value)
        integer, intent(in) :: id
        value=.false.
        if (id>=1.and.id<=ncases) value=cases(id)%has_history
    end function history_case

    subroutine response_kind(id,component,label,units)
        integer, intent(in) :: id
        integer, intent(out) :: component
        character(:), allocatable, intent(out) :: label,units

        if (id<1.or.id>ncases) error stop "Unknown annulus response kind"
        component=cases(id)%response_component
        label=trim(cases(id)%response_label)
        units=trim(cases(id)%response_units)
    end subroutine response_kind

    function latex_response(label) result(text)
        !! Describe a response measure in generated Markdown with its mathematical symbol.
        character(*), intent(in) :: label
        character(:), allocatable :: text

        select case(trim(label))
        case("sigma_rr")
            text="radial Cauchy stress \(\sigma_{rr}\)"
        case("P_rR")
            text="radial first Piola-Kirchhoff stress \(P_{rR}\)"
        case("heat_flux_radial")
            text="radial heat flux \(q_r\)"
        case("molar_flux_radial")
            text="radial molar flux \(j_r\)"
        case("current_density_radial")
            text="radial electric current density \(J_r\)"
        case("electric_displacement_radial")
            text="radial electric displacement \(D_r\)"
        case("darcy_velocity_radial")
            text="radial Darcy velocity \(v_r\)"
        case("traction_normal")
            text="normal cohesive traction \(t_n\)"
        case default
            text="`"//trim(label)//"`"
        end select
    end function latex_response

    function latex_units(units) result(text)
        !! Format the benchmark's SI labels for Markdown without changing stored unit strings.
        character(*), intent(in) :: units
        character(:), allocatable :: text

        select case(trim(units))
        case("1")
            text="1"
        case("1/K")
            text="\mathrm{K}^{-1}"
        case("1/Pa")
            text="\mathrm{Pa}^{-1}"
        case("A/m^2")
            text="\mathrm{A}\,\mathrm{m}^{-2}"
        case("C/m2","C/m^2")
            text="\mathrm{C}\,\mathrm{m}^{-2}"
        case("F*m/V2")
            text="\mathrm{F}\,\mathrm{m}\,\mathrm{V}^{-2}"
        case("F/m")
            text="\mathrm{F}\,\mathrm{m}^{-1}"
        case("K")
            text="\mathrm{K}"
        case("Pa")
            text="\mathrm{Pa}"
        case("Pa*s","Pa.s")
            text="\mathrm{Pa}\,\mathrm{s}"
        case("Pa*s/m2")
            text="\mathrm{Pa}\,\mathrm{s}\,\mathrm{m}^{-2}"
        case("Pa*s2/m3")
            text="\mathrm{Pa}\,\mathrm{s}^{2}\,\mathrm{m}^{-3}"
        case("Pa/m")
            text="\mathrm{Pa}\,\mathrm{m}^{-1}"
        case("S/m")
            text="\mathrm{S}\,\mathrm{m}^{-1}"
        case("W/(m*K)")
            text="\mathrm{W}\,\mathrm{m}^{-1}\,\mathrm{K}^{-1}"
        case("W/m^2")
            text="\mathrm{W}\,\mathrm{m}^{-2}"
        case("m")
            text="\mathrm{m}"
        case("m/s")
            text="\mathrm{m}\,\mathrm{s}^{-1}"
        case("m2")
            text="\mathrm{m}^{2}"
        case("m2/s")
            text="\mathrm{m}^{2}\,\mathrm{s}^{-1}"
        case("mol/(m^2 s)")
            text="\mathrm{mol}\,\mathrm{m}^{-2}\,\mathrm{s}^{-1}"
        case("s")
            text="\mathrm{s}"
        case default
            text="`"//trim(units)//"`"
            return
        end select
        text="\("//text//"\)"
    end function latex_units

    function latex_number(value) result(text)
        !! Scientific notation for Markdown tables, retaining five decimal places in the mantissa.
        real(wp), intent(in) :: value
        character(:), allocatable :: text
        character(24) :: number
        character(8) :: power
        integer :: marker,power_value

        write(number,"(es18.5e4)") value
        number=adjustl(number)
        marker=index(number,"E")
        if (marker==0) then
            text="`"//trim(number)//"`"
            return
        end if
        read(number(marker+1:),*) power_value
        write(power,"(i0)") power_value
        text="\("//number(:marker-1)//"\times10^{"//trim(power)//"}\)"
    end function latex_number
end module benchmark_models

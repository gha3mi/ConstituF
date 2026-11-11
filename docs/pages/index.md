---
title: Materials catalog
---

Choose a model by physical domain, then follow its FORD module link for
the governing equations, parameter domains, physical dimensions, discrete update and
derivative contracts. The table uses the exact declared field keys.
Packed component counts and ordering are defined in
[field conventions](|page|/conventions.html).

To evaluate a model, begin with [getting started](|page|/getting-started.html)
and its linked [material example](|page|/material-examples.html).
For solver coupling, see [FEM integration](|page|/fem-integration.html).
For reproducible response contours and analytical sensitivity heatmaps,
see the [benchmark studies](|page|/benchmarks.html).

## Solid and interface mechanics

| Model and FORD reference | Declared input → output | Distinguishing feature | Example |
| --- | --- | --- | --- |
| Isotropic elasticity — [[constituf_isotropic_elastic(module)]] | `strain.small` → `stress.cauchy` | Linear, isotropic small strain | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/isotropic_elastic.f90) |
| Plane elasticity — [[constituf_plane_elastic(module)]] | `strain.small` → `stress.cauchy` | Isotropic plane stress or plane strain with 2D Mandel fields | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/plane_elastic.f90) |
| Axisymmetric elasticity — [[constituf_axisymmetric_elastic(module)]] | `strain.small` → `stress.cauchy` | Axisymmetric Mandel response with four components | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/axisymmetric_elastic.f90) |
| Orthotropic elasticity — [[constituf_orthotropic_elastic(module)]] | `strain.small` → `stress.cauchy` | Linear, material principal axes | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/orthotropic_elastic.f90) |
| Neo-Hookean — [[constituf_neo_hookean(module)]] | `deformation.f` → `stress.piola` | Compressible finite strain | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/neo_hookean.f90) |
| Mooney–Rivlin — [[constituf_mooney_rivlin(module)]] | `deformation.f` → `stress.piola` | Compressible finite strain | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/mooney_rivlin.f90) |
| Gent — [[constituf_gent(module)]] | `deformation.f` → `stress.piola` | Compressible finite strain, limiting chain stretch | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/gent.f90) |
| Fiber-reinforced HGO — [[constituf_fiber_hgo(module)]] | `deformation.f` → `stress.piola`, `energy.helmholtz` | Two-fiber hyperelasticity with Helmholtz energy output | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/fiber_hgo.f90) |
| Yeoh — [[constituf_yeoh(module)]] | `deformation.f` → `stress.piola` | Three-term modified-invariant hyperelasticity | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/yeoh.f90) |
| Ogden — [[constituf_ogden(module)]] | `deformation.f` → `stress.piola` | One-term spectral hyperelasticity, general exponent | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/ogden.f90) |
| \(J_2\) plasticity — [[constituf_j2_plasticity(module)]] | `strain.small` → `stress.cauchy` | Associative plasticity with hardening | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/j2_plasticity.f90) |
| Scalar ductile damage — [[constituf_ductile_damage(module)]] | `strain.small` → `stress.cauchy`, `energy.helmholtz`, `energy.dissipation_increment` | \(J_2\) plasticity with strain-threshold damage and energy outputs | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/ductile_damage.f90) |
| Perzyna viscoplasticity — [[constituf_perzyna(module)]] | `strain.small` → `stress.cauchy` | Rate-dependent \(J_2\) return with implicit substepping | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/perzyna.f90) |
| Plane \(J_2\) plasticity — [[constituf_j2_plane(module)]] | `strain.small` → `stress.cauchy` | Small-strain \(J_2\) with plane-stress constraint or plane strain | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/j2_plane.f90) |
| Axisymmetric \(J_2\) plasticity — [[constituf_axisymmetric_j2(module)]] | `strain.small` → `stress.cauchy` | Axisymmetric \(J_2\) return with history-dependent sensitivities | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/axisymmetric_j2.f90) |
| Hill 1948 plasticity — [[constituf_hill_plasticity(module)]] | `strain.small` → `stress.cauchy` | Anisotropic yielding with history-dependent parameter sensitivities | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/hill_plasticity.f90) |
| Chaboche — [[constituf_chaboche(module)]] | `strain.small` → `stress.cauchy` | Single-backstress cyclic kinematic hardening | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/chaboche.f90) |
| Drucker–Prager — [[constituf_drucker_prager(module)]] | `strain.small` → `stress.cauchy` | Associated pressure-sensitive cone with apex return | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/drucker_prager.f90) |
| Kelvin–Voigt — [[constituf_kelvin_voigt(module)]] | `strain.small` → `stress.cauchy` | Rate-dependent spring and dashpot | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/kelvin_voigt.f90) |
| Multi-branch Maxwell — [[constituf_maxwell_multi(module)]] | `strain.small` → `stress.cauchy` | Multiple relaxation times, full-history derivatives | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/maxwell_multi.f90) |
| Single-branch Maxwell — [[constituf_maxwell_shear(module)]] | `strain.small` → `stress.cauchy` | Shear relaxation with history | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/maxwell_shear.f90) |
| Mode I cohesive — [[constituf_cohesive_mode_i(module)]] | `opening.normal` → `traction.normal` | Tension-only bilinear interface softening | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/cohesive_mode_i.f90) |
| Mixed-mode cohesive — [[constituf_cohesive_mixed_mode(module)]] | `opening.jump` → `traction.jump` | Effective-opening damage, frictionless contact | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/cohesive_mixed_mode.f90) |
| Isotropic thermoelasticity — [[constituf_isotropic_thermoelastic(module)]] | `strain.small`, `temperature.absolute` → `stress.cauchy` | Small-strain thermal expansion | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/isotropic_thermoelastic.f90) |

## Heat and species transport

| Model and FORD reference | Declared input → output | Distinguishing feature | Example |
| --- | --- | --- | --- |
| Isotropic Fourier — [[constituf_fourier_isotropic(module)]] | `temperature.gradient` → `heat.flux` | Constant isotropic conductivity | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/fourier_isotropic.f90) |
| Anisotropic Fourier — [[constituf_fourier_anisotropic(module)]] | `temperature.gradient` → `heat.flux` | Symmetric conductivity tensor | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/fourier_anisotropic.f90) |
| Temperature-dependent Fourier — [[constituf_fourier_exponential(module)]] | `temperature.gradient`, `temperature.absolute` → `heat.flux` | Exponential conductivity | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/fourier_exponential.f90) |
| Isotropic Fick — [[constituf_fick_isotropic(module)]] | `concentration.gradient` → `species.molar_flux` | Isotropic molar diffusion | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/fick_isotropic.f90) |
| Anisotropic Fick — [[constituf_fick_anisotropic(module)]] | `concentration.gradient` → `species.molar_flux` | Symmetric diffusivity tensor | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/fick_anisotropic.f90) |

## Electrical and electrostatic response

| Model and FORD reference | Declared input → output | Distinguishing feature | Example |
| --- | --- | --- | --- |
| Isotropic Ohmic — [[constituf_ohmic_isotropic(module)]] | `electric.field` → `electric.current_density` | Isotropic conduction | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/ohmic_isotropic.f90) |
| Anisotropic Ohmic — [[constituf_ohmic_anisotropic(module)]] | `electric.field` → `electric.current_density` | Symmetric conductivity tensor | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/ohmic_anisotropic.f90) |
| Isotropic dielectric — [[constituf_dielectric_isotropic(module)]] | `electric.field` → `electric.displacement` | Linear electrostatics | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/dielectric_isotropic.f90) |
| Anisotropic dielectric — [[constituf_dielectric_anisotropic(module)]] | `electric.field` → `electric.displacement` | Symmetric permittivity tensor | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/dielectric_anisotropic.f90) |
| Cubic-field dielectric — [[constituf_dielectric_cubic(module)]] | `electric.field` → `electric.displacement` | Nonlinear, reversible electrostatics | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/dielectric_cubic.f90) |

## Coupled piezoelectric response

| Model and FORD reference | Declared input → output | Distinguishing feature | Example |
| --- | --- | --- | --- |
| Linear piezoelectricity — [[constituf_piezoelectric_linear(module)]] | `strain.small`, `electric.field` → `stress.cauchy`, `electric.displacement`, `energy.electric_enthalpy` | Stress-charge coupling with electric enthalpy output | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/piezoelectric_linear.f90) |

## Porous media flow

| Model and FORD reference | Declared input → output | Distinguishing feature | Example |
| --- | --- | --- | --- |
| Isotropic Darcy — [[constituf_darcy_isotropic(module)]] | `pressure.gradient` → `fluid.darcy_velocity` | Constant scalar permeability | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/darcy_isotropic.f90) |
| Anisotropic Darcy — [[constituf_darcy_anisotropic(module)]] | `pressure.gradient` → `fluid.darcy_velocity` | Symmetric permeability tensor | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/darcy_anisotropic.f90) |
| Darcy–Forchheimer — [[constituf_darcy_forchheimer(module)]] | `pressure.gradient` → `fluid.darcy_velocity` | Nonlinear inertial resistance | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/darcy_forchheimer.f90) |
| Pressure-dependent Darcy — [[constituf_darcy_exponential(module)]] | `pressure.gradient`, `pressure.absolute` → `fluid.darcy_velocity` | Exponential pressure dependence | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/darcy_exponential.f90) |

## Coupled poromechanics

| Model and FORD reference | Declared input → output | Distinguishing feature | Example |
| --- | --- | --- | --- |
| Linear Biot poroelasticity — [[constituf_biot_poroelastic(module)]] | `strain.small`, `pressure.pore` → `stress.total`, `fluid.content` | Coupled total stress and fluid-content response | [Fortran](https://github.com/gha3mi/ConstituF/blob/main/example/biot_poroelastic.f90) |

## Use or extend a model

Construct a [[constituf_material(module):material_t(type)]], evaluate a
trial and commit only accepted solver steps. Time increments, initial
fields and derivative capabilities depend on the chosen model.
See [adding a material](|page|/adding-materials.html) to extend the
catalog and [working precision](|page|/precision.html) for build options.

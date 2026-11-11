---
title: Reading material examples
---

Each of the 40 material models has an independent Fortran program in
[`example/`](https://github.com/gha3mi/ConstituF/tree/main/example).
The programs show how to call ConstituF. The material laws themselves
are implemented in `src/materials/`. The separate
[`calibrate_elastic` example](https://github.com/gha3mi/ConstituF/blob/main/example/calibrate_elastic.f90)
demonstrates parameter fitting.

## Run an example

```console
pixi run -e gnu fpm run --example isotropic_elastic --profile debug
pixi run -e gnu fpm run --example j2_plasticity --profile debug
```

Run all programs with `pixi run -e gnu example`.

## Units chosen by the examples

These programs choose SI units. Their comments identify parameter and input
units and their output labels identify responses, tangent blocks and selected
parameter derivatives. Strains, deformation gradients and dimensionless
material coefficients are marked `[-]`.

ConstituF stores only the numerical values. A caller can use another
consistent system without a unit selector or conversion call. Keep all
parameters, fields, time increments and accepted history consistent.
Tangent units are output units divided by input units. Parameter and
history derivative units are the differentiated quantity's units divided
by the selected parameter's units.

## Common evaluation sequence

Every material example uses the same public `material_t` interface:

1. Construct a model and call `init` for one point. If requesting
   analytical sensitivities, name the scalar parameters in
   `sensitivity_parameters`.
2. Set the physical input, call `eval` and check `status == cf_ok`
   before reading the response, optional `tangent`, or `dsdp`.
3. For history-dependent materials, `commit` each accepted increment.
   Inspect accepted history with `export_committed` and use `revert`
   to discard a rejected trial. Supply `dt` and `initial_input`
   when the model requires them.

All examples use one material point, but the same calls accept multiple
points without a different single-point API.

## Array shapes

| Array | One-point shape | Meaning |
| --- | --- | --- |
| Input | `(ninput,1)` | All declared driving fields in packed order |
| Response | `(noutput,1)` | All declared physical output fields |
| `tangent` | `(noutput,ninput,1)` | Jacobian of the current trial response with respect to the current input |
| `dsdp` | `(noutput,nselected,1)` | Analytical response derivatives for selected parameters |
| Accepted `state` | `(nstate,1)` | Real internal variables, if declared |
| `sensitivity_state` | `(nstate,nselected,1)` | Accepted parameter derivatives of real state |

The last index is the point index. For models with discrete state,
additional integer or logical history arrays may be needed during export.
Some rate models retain a `previous_input` rather than a nonempty real
state vector.

The tangent holds committed history fixed. `dsdp` includes derivatives
propagated through accepted constitutive history at fixed prescribed
input/time-step history. Neither represents the derivative of a
parameter-dependent global FEM solution.

## Track more than one parameter

Example programs request one representative derivative for readability.
To select both isotropic elastic moduli, use:

```fortran
real(wp) :: dsdp(6,2,1)

call material%init(isotropic_elastic_t(160.0e9_wp,80.0e9_wp),1,status, &
    sensitivity_parameters=[character(len=13) :: "bulk_modulus","shear_modulus"])
if (status /= cf_ok) error stop "Initialization failed"

call material%eval(strain,stress,tangent,status,dsdp=dsdp)
if (status /= cf_ok) error stop "Evaluation failed"
```

Here `dsdp(:,1,1)` is the bulk-modulus derivative and
`dsdp(:,2,1)` the shear-modulus derivative. Ask
`material%parameter_info` about available parameters during setup.
Each concrete module reference gives their physical dimensions and domains.
For enabled sensitivity tracking, pass `dsdp` on **every** trial,
including trials later rejected.

## Models with history

| Examples | Behavior |
| --- | --- |
| \(J_2\), Hill, Chaboche, Perzyna and damage | Accepted plastic or damage state, plus derivative history where supported |
| Cohesive laws | Maximum opening and irreversible unloading |
| Maxwell laws | Stress relaxation and retained branch overstress |
| Kelvin–Voigt | Rate calculated from the last accepted strain and a positive `dt`. No real internal-state components |

For coupled laws, the response and tangent contain several physical
fields with different units. Identify their offsets through
[the catalog](|page|/index.html) and
[field conventions](|page|/conventions.html) rather than assuming stress
always has six components.

For restart and imported state, see
[checkpointing](|page|/checkpointing.html). For integration-point
batching and failure handling, see
[FEM integration](|page|/fem-integration.html). Numerical checks live in
[`test/`](https://github.com/gha3mi/ConstituF/tree/main/test).

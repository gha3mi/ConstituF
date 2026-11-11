---
title: Public interface
---

Applications import `material_t`, the numeric kinds `wp`, `ik`
and `lk` and status constants from `constituf`. Concrete material
constructors are imported from their individual
`constituf_<model>` modules. The
[materials catalog](|page|/index.html) lists the available models.
The host application supplies the physical driving fields and handles
global equilibrium or conservation equations.

The API accepts and returns plain numerical values. The caller owns unit
consistency across parameters, fields, time increments and history. Field
metadata describes physical measures and tensor layouts without unit labels.
See [unit conventions](|page|/conventions.html) for the dimensional contract.

## Evaluation and accepted history

| Call | Purpose |
| --- | --- |
| `material%init(model,npoints,status,...)` | Initialize a model and its point history. Optional `initial_input` and `sensitivity_parameters` |
| `material%eval(x,y,tangent,status,dt=...,dsdp=...)` | Evaluate a trial response and optional Jacobians |
| `material%commit(status)` | Accept the latest successful trial, including history derivatives |
| `material%revert()` | Discard the trial without modifying committed history |

For \(N\) material points and \(m\) selected scalar parameters:

| Buffer | Shape |
| --- | --- |
| Prescribed input `x` | `(ninput,N)` |
| Material response `y` | `(noutput,N)` |
| Optional `tangent` | `(noutput,ninput,N)` |
| Optional `dsdp` | `(noutput,m,N)` |

The last index counts **material points**, including the single-point
case \(N=1\). The number and meaning of packed components depend on the
model. The algorithmic tangent differentiates the *discrete* response
with respect to the current input at fixed committed history.
Analytical `dsdp` tracks the material-parameter derivative through
accepted history at fixed prescribed inputs and time increments. These
are distinct derivatives.

Select sensitivity parameter names during `init` and pass `dsdp`
on **every** trial of that handle. Some materials require a positive
`dt` or a prescribed `initial_input`. Consult the material's
module reference. Check `status == cf_ok` before reading outputs or
committing. Failed batches yield no usable response and cannot be
committed.

The [getting-started tutorial](|page|/getting-started.html) provides a
complete one-point program. [Material examples](|page|/material-examples.html)
show stateful loading histories.

## Field and parameter queries

Use the following inquiries during setup rather than in the solver's
integration-point loop:

| Methods on `material_t` | Meaning |
| --- | --- |
| `input_size()`, `output_size()`, `point_count()` | Total packed input/output components and initialized points |
| `input_offset(name)`, `output_offset(name)` | One-based field offset. Zero if the name is absent |
| `input_components(name)`, `output_components(name)` | Number of packed values in a field |
| `field_matches(...)`, `field_info(...)` | Match or inspect the field's representation, measure, frame and physical shape |
| `parameter_count()`, `parameter_info(index,name,found)` | Declared scalar parameter names |
| `sensitivity_size()`, `sensitivity_parameter_id(index)` | Selected parameter count and declared IDs |

Field roles are `input`, `output`, `state`,
`integer_state` and `logical_state`. `field_info` allocates
returned metadata. See [field conventions](|page|/conventions.html)
for Mandel and full-tensor packing. The
[[constituf_material(module):material_t(type)]] reference documents
the complete procedure signatures and status contracts.

## State transfer and diagnostics

| Methods on `material_t` | Purpose |
| --- | --- |
| `state_size()`, `integer_state_size()`, `logical_state_size()` | Real, integer and logical state components per point |
| `previous_input_size()` | Retained prior-input components per point. Zero if none |
| `export_committed(...)`, `import_committed(...)` | Transfer accepted history between compatible handles |
| `point_statuses(codes,status)`, `last_failure_point()` | Diagnose the last batch evaluation |
| `state_storage_bytes()` | Size of retained history buffers, excluding caller arrays and some internal storage |

The handle manages state internally. Normal `eval` calls need no
caller-owned history arrays. Explicit exports are for restart or
interoperability. See [checkpointing](|page|/checkpointing.html).

## Other procedures exported by `constituf`

| Procedure | Task |
| --- | --- |
| `verify_tangent` | Check a discrete tangent against numerical trial perturbations |
| `verify_potential` | Compare energy or enthalpy gradients with conjugate responses |
| `save_checkpoint`, `load_checkpoint` | Save and restore committed state in a versioned file |
| `calibrate_parameters` | Fit parameters with a user-provided prediction/analytical-Jacobian callback |

Calibration is not automatic for every model. The application must
provide observations, a suitable experiment, model construction,
parameter bounds and an analytical response Jacobian.
See [calibration](|page|/calibration.html) and
[verification](|page|/verification.html).

Status constants include `cf_ok`, `cf_invalid_argument`,
`cf_invalid_shape`, `cf_invalid_time`, `cf_unsupported`,
`cf_nondifferentiable`, `cf_evaluation_failed`,
`cf_no_trial` and checkpoint-related errors. See
[[constituf_status(module)]] for the full list.

To implement new materials rather than use existing ones, import
[[constituf_author(module)]] and extend
[[constituf_model(module):material_model_t(type)]].

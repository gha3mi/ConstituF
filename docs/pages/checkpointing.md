---
title: Restarting committed material histories
---

This how-to covers transferring or saving accepted ConstituF
material-point state. A checkpoint preserves **committed history**, not
an unaccepted trial response. Constitutive state is owned by
[[constituf_material(module):material_t(type)]], while the host solver
owns its global solution, mesh and integration-point ordering.

## Choose a restart method

- **Portable file:** [[constituf_checkpoint(module)]] provides the
  public `save_checkpoint` and `load_checkpoint` routines. Use this
  when the solver needs a versioned on-disk record.
- **Caller-owned arrays:** `material_t%export_committed` and
  `%import_committed` copy accepted buffers between handles with no
  file I/O. Use these when the solver has its own checkpoint format.

In either case, initialize the destination handle with the **same**
model type, numerical parameters, precision, point count, field
conventions, initial conditions and selected sensitivity parameters
in the same order *before* importing the state.

## Save and restore a checkpoint file

After accepting a converged global step with `commit`, save a
versioned checkpoint:

```fortran
use constituf, only: cf_ok, save_checkpoint, load_checkpoint

! material and restarted are identically initialized material_t handles.
call save_checkpoint(material, "material.cf", &
    "maxwell_shear.K80.G15.G1_30.tau2.N12", status)
if (status /= cf_ok) error stop "Checkpoint save failed"

call load_checkpoint(restarted, "material.cf", &
    "maxwell_shear.K80.G15.G1_30.tau2.N12", status)
if (status /= cf_ok) error stop "Checkpoint load failed"
```

The configuration ID above is illustrative. Applications must
derive their own stable ID from the **complete** model configuration,
parameter values and units, point ordering and material conventions.
The library cannot reflect arbitrary user-defined model constructors.
Reusing an ID for different configurations defeats the consistency
check.

Format version 1 records working precision, layout and history
extents, selected parameter IDs in order and a text-integrity
checksum. A malformed or truncated file returns
`cf_checkpoint_corrupt`. Incompatible configuration, precision or
shapes return `cf_checkpoint_mismatch`. File I/O errors return
`cf_io_error`. The checksum detects accidental corruption, **not**
malicious modification. A rejected load does not modify existing
committed or pending trial state.

Saving replaces the destination file and is **not crash-atomic**.
Applications needing durable crash recovery should write to a distinct
temporary path and manage an atomic rename and file synchronization
themselves. Check numerical restart behavior when crossing compiler
families, even with the same working precision.

## Transfer the committed arrays directly

For a model with real internal history, allocate
`state(material%state_size(),material%point_count())`.
If `material%previous_input_size()` is nonzero, also allocate
`previous(material%previous_input_size(),material%point_count())`.

```fortran
call material%export_committed(state, status, previous_input=previous)
if (status /= cf_ok) error stop "State export failed"

call restarted%import_committed(state, status, previous_input=previous)
if (status /= cf_ok) error stop "State import failed"
```

This excerpt assumes both arrays have the correct extents and that
previous input is required. Omit `previous_input` when the model does
not retain it. Models with integer or logical state also require
`integer_state` and `logical_state` buffers. Enabled analytical
sensitivities require `sensitivity_state(nstate,nselected,npoints)`.
Missing or mismatched buffers are rejected. Direct array transfers
do **not** verify physical model identity or provide serialization.

After a successful import, pending trials are invalidated. A
rejected import leaves the live material unchanged. See
[FEM integration](|page|/fem-integration.html) for the trial/commit
lifecycle and [verification](|page|/verification.html) for the
checkpoint roundtrip regression.

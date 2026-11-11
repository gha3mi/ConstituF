---
title: Integrating with a finite-element solver
---

This how-to explains how an existing finite-element (FEM) program can call
ConstituF at quadrature points. The FEM solver remains responsible for
kinematics, global residual assembly, boundary conditions, nonlinear
iteration and convergence. ConstituF supplies local responses and, for
models that support them, local algorithmic tangent blocks.

## Connect the library

For an FEM solver using [fpm](https://fpm.fortran-lang.org/), add the
GitHub dependency to the solver's `fpm.toml`:

```toml
[dependencies]
constituf = { git = "https://github.com/gha3mi/ConstituF.git" }
```

To request a different working precision, add `features = ["real32"]`
or `features = ["real128"]` to the dependency.
Other Fortran build systems must compile and link ConstituF with a
compatible compiler, kind and module configuration. A Fortran `.mod`
file is generally compiler-specific.

The user-facing handle is [[constituf_material(module):material_t(type)]]
(exported from `constituf`). Concrete model constructors are imported
from their own modules, for example
[[constituf_isotropic_elastic(module)]]. See the
[public interface](|page|/public-interface.html) for available operations.

## Bind fields once during solver setup

1. Select the material model and construct a handle for a chosen batch of
   quadrature points. Call `material%init(model, npoints, status, ...)`.
2. Query `input_size()` and `output_size()` and allocate packed arrays of
   shape `(ninput,npoints)` and `(noutput,npoints)`.
3. For each declared field, inspect or validate its exact role, key,
   representation, measure, frame and **physical** shape. Resolve named
   offsets before entering the integration-point loop.
4. Supply all parameters, fields and history in the solver's consistent
   unit system. Map tensor conventions explicitly where necessary.
   ConstituF does not store unit labels or convert numerical units or
   tensor conventions.

For example, the isotropic small-strain material declares an input field
`strain.small`, with dimensionless strain, Mandel representation,
`small_strain` measure, `spatial` frame and physical shape `[3,3]`:

```fortran
if (.not. material%field_matches( &
    "input", "strain.small", "mandel", &
    "small_strain", "spatial", [3,3])) then
    error stop "FEM/material field mismatch"
end if
```

`field_matches` requires an initialized handle. Use `field_info` to
retrieve metadata when building generic material adapters. It allocates
returned strings and shape and therefore belongs in setup, not the hot
path. These checks do not validate numerical unit consistency.
Packed offsets are one-based. See
[conventions](|page|/conventions.html) for the Mandel and full-tensor order.

For **coupled laws**, bind every input and output field separately before
packing the entire batch. For example,
[[constituf_biot_poroelastic(module)]] declares six strain components
followed by pore pressure, then returns six total-stress components
followed by fluid content. Its local \(7\times7\) tangent includes
stress/pressure and fluid-content/strain cross-derivatives with different
physical units. The FEM solver—not ConstituF—must map these blocks into
its coupled residual and global tangent.

## Evaluate, accept or reject a step

During a trial FEM load/time step, pack each integration-point input into
`x(ninput,npoints)`, then call:

```fortran
call material%eval(x, response, tangent, status, dt=dt)
```

Here `response` has shape `(noutput,npoints)` and optional `tangent`
has shape `(noutput,ninput,npoints)`. The tangent is the derivative of the
**discrete material response** with respect to the current packed input
at fixed committed history and previous input. It is not the assembled
global FEM Jacobian.

Pass a positive `dt` for a model that requires time integration. Other
models may not require it. Some models also require `initial_input`
during `init`. Do not request an unavailable tangent or an analytical
sensitivity that the model does not implement.

All trial calls start from the same committed material state. A new
Newton iterate can replace the previous trial without committing it.

- If `status == cf_ok`, use the response and optional tangent in the solver.
  Once the **global** load/time step is accepted, call
  `material%commit(status)` and check `status == cf_ok`.
- If the solver rejects the step, call `material%revert()`. The committed
  history remains unchanged.
- If evaluation fails, discard **all** response and tangent entries.
  A failed trial cannot be committed. `last_failure_point()` identifies
  the first failing point for a kernel-level point failure. Setup or shape
  failures may instead return zero.

### Interpret status codes

Status constants are exported by `constituf` and described in
[[constituf_status(module)]]. Important cases for solver integration are:

| Status | Meaning for the caller |
| --- | --- |
| `cf_invalid_argument` | Invalid constructor parameters or setup arguments |
| `cf_invalid_shape` | Incorrect packed input, response, tangent or sensitivity extents |
| `cf_invalid_time` | Missing, nonfinite or nonpositive required time increment |
| `cf_unsupported` | Requested derivative or capability is not offered by this model |
| `cf_nondifferentiable` | A requested classical derivative is undefined at a material branch switch |
| `cf_evaluation_failed` | One or more integration-point kernels failed. Discard all trial outputs |
| `cf_no_trial` | `commit` was requested without a successful pending trial |

A solver must handle each failed step explicitly, not reuse or commit
trial data after an error. Use `last_failure_point()` for point-level
diagnosis, not as a per-point result array.

The core owns each model's real, integer and logical committed/trial
history internally. **No caller-supplied history array is required**, even
for path-dependent materials. Stateless models still use the same `init`
and `eval` interface.

## Point-level evaluation diagnostics

The Fortran handle records a per-point status vector from the latest
material kernel call. If the *batch* returns a nonzero status, discard
**all** response and tangent arrays. Use this query only to determine
which points failed:

```fortran
integer :: point_codes(npoints)
call material%point_statuses(point_codes, status)
if (status == cf_ok) then
    ! point_codes(i) is the status from the latest material kernel
end if
```

`point_statuses` returns `cf_no_trial` if no kernel has run or if
the latest `eval` failed its input-shape/time/capability checks before
calling a material kernel. It returns `cf_invalid_shape` for the wrong
buffer length. A successful checkpoint import invalidates old per-point
diagnostics. Use `last_failure_point()` to inspect only the first failure.

## Integration and performance boundaries

- Pack contiguous batches with component index first and point index second.
  The caller owns the `x`, `response` and `tangent` arrays. They must not
  alias. Match physical units and tensor measures before passing data.
- Keep each mutable `material_t` handle exclusive to one concurrent
  evaluation. Separate worker threads can own separate handles or use
  disjoint solver-managed batches.
- Request `tangent` only when necessary. Optional parameter sensitivities
  are enabled during `init(sensitivity_parameters=...)`. An enabled handle
  must receive `dsdp(noutput,nselected,npoints)` on every `eval`.
  These are derivatives at fixed prescribed input history, not derivatives
  of a parameter-dependent global solution.
- Obtain `state_storage_bytes()` when estimating retained local-history
  memory. It excludes caller-owned arrays and some internal metadata.

Consult the actual material module's FORD reference for fields,
parameter bounds, time-step requirements, differentiability and supported
derivatives. See the [\(J_2\) plasticity example](https://github.com/gha3mi/ConstituF/blob/main/example/j2_plasticity.f90)
for a committed plastic step followed by a reverted unloading trial.

## Related workflows

- [Restarting committed material histories](|page|/checkpointing.html)
  explains direct state transfer, portable files, compatibility and
  failure-atomic restore.
- [Field and numerical conventions](|page|/conventions.html) explains
  plane-stress, axisymmetric, coupled, energy-conjugate and damage-model
  output layouts. Each concrete material module remains the authority
  on its supported fields and derivatives.
- [Calibrating material parameters](|page|/calibration.html) covers the
  offline analytical-sensitivity least-squares utility.
- [Verification and reproducibility](|page|/verification.html) covers
  tangent checks, work-conjugate energy checks and independent tests.
- [Choosing working precision](|page|/precision.html) covers a
  consistent `real(wp)` build of the material and the host solver.


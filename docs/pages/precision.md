---
title: Choosing working precision
---

ConstituF selects its Fortran working real kind **at build time**.
All material modules and any application that uses `real(wp)` must be
built with a compatible kind. Precision is not a runtime setting.

## Select a precision

The default `wp` is `real64`. ConstituF provides two alternative
fpm features, `real32` and `real128`, corresponding to the
`debug-real32`, `release-real32`, `debug-real128` and
`release-real128` build profiles in [`fpm.toml`](https://github.com/gha3mi/ConstituF/blob/main/fpm.toml).
Do not enable both features.

For an application depending on ConstituF through Git, select one
feature in the application's `fpm.toml`:

```toml
[dependencies]
constituf = { git = "https://github.com/gha3mi/ConstituF.git", features = ["real32"] }
```

Use `features = ["real128"]` for extended precision, or omit
`features` for `real64`. Add a published `tag` or exact `rev`
for reproducible dependency resolution.

## Check the selected build

From a ConstituF checkout, run the kind, elastic-response and checkpoint
roundtrip smoke tests for the selected precision:

```console
pixi run -e gnu fpm test --profile debug-real32 --target test_core_kinds test_isotropic_elastic test_core_checkpoint_file
```

For a compiler that supports `real128`, use:

```console
pixi run -e gnu fpm test --profile debug-real128 --target test_core_kinds test_isotropic_elastic test_core_checkpoint_file
```

These targeted tests do **not** establish precision-specific validation
of every material model. Recheck tolerances, valid parameter scales,
convergence and physical results for the selected application. The
`real128` kind depends on compiler support and is not guaranteed by
every platform.

## Rebuild dependent code

When the selected kind changes, cleanly rebuild ConstituF, its model
modules and the consuming solver together. Do not mix Fortran module
files, array kinds or checkpoint snapshots produced with different
precision configurations. For public kind definitions, see
[[constituf_kinds(module)]]. For restart compatibility, see
[checkpointing](|page|/checkpointing.html).

Return to [getting started](|page|/getting-started.html) for the
first material-point evaluation.

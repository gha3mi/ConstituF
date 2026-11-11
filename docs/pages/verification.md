---
title: Verification
---

Verification asks whether a material implementation reproduces its stated
constitutive law and discrete update. Validation requires separate evidence
that the law describes a physical material in its intended regime.
Plausible response figures are not experimental validation.

## Constitutive responses and histories

Compare the response with an independent analytical result, a published
example or a separately derived limiting case. State the field measures,
units, parameter values, loading history and justified numerical tolerances.

For stateful models, verify the accepted internal variables after each
increment. Check that `revert` discards trial history and that failed
evaluations do not change previously committed state.

## Algorithmic tangents

The local algorithmic Jacobian differentiates the **discrete** trial
response \(\boldsymbol y\) with respect to the current prescribed input
\(\boldsymbol x\), at fixed committed state \(\boldsymbol s_n\),
previous input \(\boldsymbol x_n\), parameters \(\boldsymbol p\)
and time increment \(\Delta t\):

\[
J_{ij}^{\mathrm{alg}} =
\left.\frac{\partial y_i}{\partial x_j}\right|_{
\boldsymbol s_n,\boldsymbol x_n,\boldsymbol p,\Delta t}.
\]

The public [[constituf_verify(module)]] routine compares this matrix
with centered perturbations:

\[
J_{ij}^{\mathrm{FD}}\approx
\frac{y_i(\boldsymbol x+h_j\boldsymbol e_j)
-y_i(\boldsymbol x-h_j\boldsymbol e_j)}{2h_j}.
\]

For an initialized `material` and compatible packed input `x`:

```fortran
use constituf, only: wp, cf_ok, verify_tangent

real(wp) :: scaled_error
integer :: status

call verify_tangent(material, x, scaled_error, status)
if (status /= cf_ok) error stop "Tangent check failed"
```

A successful call returns `scaled_error` no greater than \(1\), where the error is
normalized by the specified absolute/relative tolerances. Supply
`dt=dt` for time-dependent materials and a physically appropriate
`step(ninput)` for mixed-unit inputs. Verification allocates workspace,
performs repeated material trials and discards pending trials. It
does not generate production tangents.

Default perturbations and error normalization use a numerical reference
magnitude of one in the caller's representation. Choose explicit steps and
tolerances when changing units or comparing fields with different scales.
The verification helpers do not infer a unit system.

## Analytical parameter derivatives

For a parameter \(p_j\), the `dsdp` result includes derivatives
propagated through all **accepted** constitutive history steps, at fixed
prescribed input and time-step history. To check it independently,
perturb \(p_j\), reinitialize the material and replay the complete
accepted loading path for each parameter perturbation. Changing only
the final step omits derivatives of previous internal variables.

Avoid yield-surface, contact, damage and other branch switches where
a classical derivative may not exist. Finite differences here are
verification comparisons, not the library's sensitivity method.

## Energy and enthalpy conjugacy

[[constituf_potential_verify(module)]] provides `verify_potential`
to compare a declared scalar potential \(\Phi\) with its conjugate
response components \(y_{c_j}\). The expected relation is

\[
\frac{\partial\Phi}{\partial x_j}=s_j\,y_{c_j},
\qquad s_j\in\{-1,+1\}.
\]

For a hyperelastic Helmholtz energy,
\(\partial\Psi/\partial F_{iJ}=P_{iJ}\). For electric enthalpy,
\(\partial\mathcal H/\partial E_i=-D_i\).
The caller supplies the potential output index, conjugate output
indices and signs. Use input-specific perturbation scales where
physical units differ. The routine does not commit trials.

## Tests and reproducibility

```console
pixi run lint
pixi run -e gnu test
pixi run -e gnu example
```

The [test sources](https://github.com/gha3mi/ConstituF/tree/main/test)
contain independent material checks, constitutive limiting cases,
algorithmic tangents, history derivatives, field metadata, checkpoint
roundtrips and failure handling. Each test defines its own numerical
criteria. On supported platforms, CI also runs Intel Fortran and
targeted precision checks. These do not establish experimental
validation or full-model testing at every precision.

For reproducible prescribed-field response and sensitivity figures:

```console
pixi run -e gnu benchmarks
pixi run docs
```

The [annulus and spherical-cavity benchmarks](|page|/benchmarks.html)
demonstrate reproducible material-point evaluation under prescribed
fields. Independent constitutive-law verification is provided by tests,
not by comparing colors or plots across materials. See the
[materials catalog](|page|/index.html) for model-specific equations
and limitations.

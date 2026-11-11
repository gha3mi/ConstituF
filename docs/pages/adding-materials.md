---
title: Adding a material model
---

ConstituF separates the common material-point interface from
individual constitutive laws. New material models extend
[[constituf_model(module):material_model_t(type)]] through the
[[constituf_author(module)]] author API. The core requires no
model-specific dispatch branch.

## 1. Define the law and its fields

Identify an original publication or authoritative textbook. Specify:

- Independent and response fields, stress or flux measures, reference
  frames, physical dimensions, tensor representations and sign conventions.
- Parameters and admissible ranges, together with meaningful limiting
  cases and constraints on time increments or inputs.
- Internal variables, reference conditions and the discrete update,
  including its behavior on rejected steps.
- Consistent current-input tangents, optional analytical parameter
  sensitivities and branches where derivatives are not defined.

For internal state \(\boldsymbol s_n\), current input
\(\boldsymbol x_{n+1}\), previous input \(\boldsymbol x_n\),
parameters \(\boldsymbol p\) and time increment \(\Delta t\),
the model implements a local transition of the form

\[
\begin{aligned}
\boldsymbol s_{n+1}
 &=\boldsymbol G(\boldsymbol x_{n+1},\boldsymbol x_n,
   \boldsymbol s_n,\boldsymbol p,\Delta t),\\
\boldsymbol y_{n+1}
 &=\boldsymbol H(\boldsymbol x_{n+1},\boldsymbol x_n,
   \boldsymbol s_n,\boldsymbol p,\Delta t).
\end{aligned}
\]

Distinguish the published continuum law from its implemented
discretization and any model-specific adaptations.

## 2. Implement the model module

Create `src/materials/<family>/constituf_<model>.f90` with a
matching module name and a public derived type extending
`material_model_t`.

| Required binding | Responsibility |
| --- | --- |
| `describe(self,spec,status)` | Validate constructor parameters and declare all fields, scalar parameters, state types and supported capabilities |
| `eval_batch(...)` | Compute trial outputs and states for every point and set a point status, without modifying committed history |

Use [[constituf_spec(module):material_spec_t(type)]] to declare
field metadata: exact keys, packed representation, physical
measure, frame and shape. A material is not assumed to have
six strain inputs or six stress outputs.

The API carries no unit labels. Document dimensional relationships and
evaluate the caller's numerical values directly in their consistent unit
system. Dimensional constants and tolerances must follow that system.

Override `initialize_state` if the default zero state is not
appropriate. Integer and logical history uses the typed hooks.
For analytical sensitivities, implement the sensitivity evaluation
hook, including derivatives propagated through accepted state history.
Do not report unsupported capabilities as available.

For implementation examples, compare the stateless
[[constituf_fourier_isotropic(module)]] with the history-dependent
[[constituf_maxwell_shear(module)]].

## 3. Document the formulation

Write FORD comments on the public model and constructor, defining
the governing equations in LaTeX, all symbols and physical dimensions, parameter
domains, discrete algorithm, tangent and sensitivity definitions,
nonsmooth cases and known limitations. Include original
bibliographic references and identify any adaptations.

Add the model to the [materials catalog](|page|/index.html)
with its exact input/output field names and module reference.

## 4. Add an example and independent tests

Add `example/<model>.f90` that initializes `material_t`,
checks evaluation status and shows accepted history where applicable.
Choose and label a consistent unit system in the example and tests so
readers can interpret the numerical values and derivatives.
Use [Fourier's example](https://github.com/gha3mi/ConstituF/blob/main/example/fourier_isotropic.f90)
as a minimal pattern.

Add `test/test_<model>.f90` with independent response or limiting-case
checks, invalid parameter cases, tangent comparisons and, when
appropriate, commit/revert and full-history parameter sensitivity
tests. See [verification](|page|/verification.html).

Record code provenance and confirm the right to redistribute any
reused third-party implementation. Citing a paper does not grant
permission to copy source code.

## 5. Register a benchmark case

Register the material in
[`benchmark/benchmark_models.f90`](https://github.com/gha3mi/ConstituF/blob/main/benchmark/benchmark_models.f90)
with its native constructor, nominal parameters in declared order,
input kind, projected response measure, benchmark-owned SI labels for
parameters and responses and colormap.
Both the [annulus](|page|/benchmark-annulus.html) and
[3D spherical-cavity](|page|/benchmark-cavity3d.html) benchmarks share
this catalog and use ConstituF to evaluate the response and its analytical
parameter derivatives. Plane and axisymmetric laws remain on the annulus.
If the law introduces a new input layout or output projection, extend the
prescribed fields in
[`benchmark_annulus_fields.f90`](https://github.com/gha3mi/ConstituF/blob/main/benchmark/benchmark_annulus_fields.f90)
and, when the law is compatible with three dimensions, in
[`benchmark_cavity3d_fields.f90`](https://github.com/gha3mi/ConstituF/blob/main/benchmark/benchmark_cavity3d_fields.f90).
These benchmarks do not reimplement constitutive laws.

## 6. Run the checks

```console
pixi run lint
pixi run -e gnu test
pixi run -e gnu example
pixi run -e gnu benchmarks
pixi run docs
```

Inspect the generated FORD module, rendered equations, catalog and
benchmark links. The [contributor instructions](https://github.com/gha3mi/ConstituF/blob/main/AGENTS.md)
describe remaining implementation and licensing requirements.

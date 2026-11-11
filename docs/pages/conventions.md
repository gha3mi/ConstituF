---
title: Field and numerical conventions
---

This cross-model reference defines the common library-level conventions.
Model-specific parameter restrictions and physical measures remain
authoritative in the corresponding source module's FORD reference.
ConstituF performs **local** constitutive evaluation, not global FEM
equilibrium or transport assembly.

## Packed fields and metadata

All batch inputs and outputs use `(component,point)` order, with the
point index last. Each declared field carries a key, physical measure,
frame, storage convention and physical shape. The number
of packed components can differ from the physical tensor size.
Several fields are concatenated in their declared order.
The mappings below use vector components \(v_i\) and tensor components
\(A_{ij}\).

| Representation | Packed components | Typical mapping |
| --- | --- | --- |
| `scalar` | 1 | One scalar per point |
| `vector` | 3 for the existing 3D models | \(v_1,v_2,v_3\) |
| `mandel` | 3 for 2D symmetric tensors | \(A_{11},A_{22},\sqrt{2}A_{12}\). Physical shape `[2,2]` |
| `mandel` | 6 for 3D symmetric tensors | \(A_{11},A_{22},A_{33},\sqrt{2}A_{23},\sqrt{2}A_{13},\sqrt{2}A_{12}\) |
| `custom.axisymmetric_mandel` | 4 | \(A_{rr},A_{zz},A_{\theta\theta},\sqrt{2}A_{rz}\). Physical shape `[3,3]` |
| `tensor_full` | 9 for the existing 3D full tensors | \(A_{11},A_{21},A_{31},A_{12},A_{22},A_{32},A_{13},A_{23},A_{33}\) |

For a symmetric three-dimensional tensor \(\boldsymbol A\), Mandel
storage is
\[
\operatorname{mandel}(\boldsymbol A)=
[A_{11},A_{22},A_{33},\sqrt{2}A_{23},\sqrt{2}A_{13},\sqrt{2}A_{12}]^{\mathsf T}.
\]
This preserves the tensor inner product:
\(\boldsymbol A:\boldsymbol B
=\operatorname{mandel}(\boldsymbol A)^{\mathsf T}
\operatorname{mandel}(\boldsymbol B)\).
Engineering shear strains instead use a factor of two, so they must
not be passed as Mandel components. Full tensors use Fortran
column-major packing \(i+3(j-1)\) for components \((i,j)\).
Finite-strain deformation gradients use this full-tensor layout.

A field's representation, measure and frame are matched as exact strings
at setup. Units are not part of the metadata. For example, the small-strain
model declares `strain.small` with convention `"mandel"`, measure `"small_strain"`,
frame `"spatial"` and physical shape `[3,3]`.
Neo-Hookean `deformation.f` declares `"tensor_full"`,
`"deformation_gradient"` and `"spatial_reference"` with physical
shape `[3,3]`. These different metadata sets must not be conflated.

Use `input_offset`, `output_offset` and `field_matches` on
[[constituf_material(module):material_t(type)]] when binding solver fields.
The offsets are one-based and zero means the name is not declared.
`field_info` returns the exact metadata when an adapter must discover it.

## Reduced-dimensional and multi-field models

[[constituf_plane_elastic(module)]] and
[[constituf_j2_plane(module)]] accept the three in-plane Mandel
components, but **plane strain** constrains out-of-plane strain
whereas **plane stress** enforces zero out-of-plane traction.
In nonlinear plane-stress plasticity, the out-of-plane strain is
solved locally, and the tangent is a discrete Schur complement,
not a truncated three-dimensional tangent. The two-dimensional
response omits \(\sigma_{33}\).

[[constituf_axisymmetric_elastic(module)]] and
[[constituf_axisymmetric_j2(module)]] use the four-component
`custom.axisymmetric_mandel` representation shown above.
Axisymmetric solvers must provide hoop strain, integration weights
and regularity near the symmetry axis. No FEM geometry is supplied
by the material model. [[constituf_perzyna(module)]] is a separate
rate-dependent model requiring a positive `dt` and initial strain.
Its substep limit and tangent contract are documented at the source.

Some models return **multiple output fields** instead of a stress-only
vector. [[constituf_fiber_hgo(module)]] provides first Piola stress
and Helmholtz energy (local tangent shape \(10\times9\)).
[[constituf_piezoelectric_linear(module)]] couples stress, electric
displacement and electric enthalpy. Its declared potential has
\(\partial\mathcal H/\partial\boldsymbol E=-\boldsymbol D\),
not \(+\boldsymbol D\). [[constituf_ductile_damage(module)]]
returns Cauchy stress, Helmholtz energy and a nonnegative
**dissipation-increment diagnostic**, not an exact finite-increment
energy balance. Its local softening law does not by itself provide
mesh-objective fracture regularization.

For every such model, bind fields by their exact metadata and
differentiate outputs with respect to the corresponding declared
inputs. Refer to the model's FORD module reference for limitations
rather than assuming a universal strain/stress layout.

## Coupled fields and tangent blocks

A material can declare several inputs and outputs of different physical
types. For example, [[constituf_biot_poroelastic(module)]] declares
`strain.small` (six Mandel strain components) followed by
`pressure.pore` (one pressure scalar) and returns `stress.total`
(six Mandel Cauchy stress components in the same pressure unit) followed by
`fluid.content` (one dimensionless scalar). For a batch of \(N\) points,

- input `x` and output `response` each have shape `(7,N)`.
- optional tangent has shape `(7,7,N)`, ordered by the declared fields.
- rows 1–6/column 7 represent stress differentiated with respect to pore
  pressure (dimensionless), while row 7/columns 1–6 represent fluid
  content differentiated with respect to strain (dimensionless).
- entry \((7,7)\) differentiates fluid content with respect to pressure
  and has inverse-pressure units, for example \(\mathrm{Pa}^{-1}\) when the
  caller chooses \(\mathrm{Pa}\) for pore pressure and stress.

The signs of the cross-blocks follow the **total-stress tension-positive**
and **pore-pressure compression-positive** conventions of this model.
This local constitutive Jacobian does not assemble the global displacement/
pressure residual or implement fluid balance. See
[FEM integration](|page|/fem-integration.html).

## Physical units and signs

ConstituF has a unitless numerical API. It stores no unit labels and
performs no unit matching or conversion. The caller must use a consistent
unit system for material parameters, prescribed fields, time increments
and retained history. Physical quantities keep the dimensions specified
by each material's equations and documentation.

For example, elastic moduli supplied in \(\mathrm{MPa}\) give stress and stress/strain
tangents in \(\mathrm{MPa}\). A relaxation time supplied in milliseconds requires
`dt` in milliseconds. Cohesive traction in \(\mathrm{MPa}\) and opening in \(\mathrm{mm}\) require
stiffness in \(\mathrm{MPa}/\mathrm{mm}\). Dimensionless quantities, such as strain and the
deformation gradient, keep their ordinary numerical values.

Each tangent entry has output units divided by input units. Each parameter
sensitivity has output units divided by the selected parameter's units.
Coupled arrays can therefore contain blocks with different dimensions.
Preserve the chosen unit convention when transferring history or restoring
a checkpoint. Metadata checks cannot detect inconsistent numerical units.

Absolute temperature remains a physical requirement where a model requests
`temperature.absolute`. Use a scale measured from absolute zero and express
temperature-dependent parameters consistently. A Celsius temperature value
cannot be substituted for an absolute temperature. Temperature differences
and temperature gradients have the interval units of the chosen scale.

The supplied examples and benchmarks choose SI units and label their data:
\(\mathrm{Pa}\) for stress, \(\mathrm{m}\) for opening, \(\mathrm{s}\) for time,
\(\mathrm{K}\) for absolute temperature, \(\mathrm{W}\,\mathrm{m}^{-2}\) for heat flux,
\(\mathrm{mol}\,\mathrm{m}^{-2}\,\mathrm{s}^{-1}\) for molar flux,
\(\mathrm{V}\,\mathrm{m}^{-1}\) for electric field and
\(\mathrm{m}\,\mathrm{s}^{-1}\) for Darcy's specific discharge. These are choices of the programs,
not requirements imposed by the library.

Example conventions from the material implementations:

- Small-strain elasticity uses tensile-positive Cauchy stress and
  infinitesimal strain in the shared spatial Cartesian frame.
- Finite-strain hyperelasticity uses deformation gradient
  \(\boldsymbol F\) and first Piola--Kirchhoff stress
  \(\boldsymbol P\), work-conjugate with respect to reference volume.
- Fourier heat flux follows \(\boldsymbol q=-\boldsymbol K\nabla T\).
  Fick diffusion and Darcy flow also return down-gradient responses.
  Electrical current density follows its model's declared convention,
  not an assumed Fourier-style sign rule.

Check each module before making a thermodynamic interpretation:
heat flux, dissipative power, energy density, strain and stress are
not interchangeable output measures.

## History, time and derivatives

A successful `eval` creates a **trial** response and trial history from
the committed state. `commit` accepts it and `revert` discards it.
ConstituF owns real, integer and logical history buffers internally.
Stateless models have no meaningful evolving internal state.
Some rate-dependent models require a positive `dt` and
`initial_input`. The concrete model documents those requirements.

For packed response \(\boldsymbol y\) and current input
\(\boldsymbol x_{n+1}\), the optional **general** algorithmic
Jacobian (not necessarily a fourth-order elastic modulus) is

\[
\boldsymbol J_{\mathrm{alg}} =
\left.\frac{\partial\boldsymbol y_{n+1}}
{\partial\boldsymbol x_{n+1}}\right|_{
\boldsymbol s_n,\boldsymbol x_n,\boldsymbol p,\Delta t}.
\]

It is packed as `(noutput,ninput,npoints)` and includes rectangular
blocks for coupled fields. It is a **discrete local derivative**, not
necessarily a continuum modulus or the global FEM Jacobian.

Analytical parameter sensitivities are separately opt-in at initialization
and are returned in `dsdp(noutput,nselected,npoints)`. They propagate
through accepted local material history at fixed prescribed inputs and
time increments. They do **not** account for design-parameter dependence
of global solver fields. Model authors must report nonsmooth derivatives
and unsupported capabilities rather than silently replacing them with
finite differences.

See the underlying
[[constituf_model(module):material_model_t(type)]] for the exact batch
update and derivative contracts and
[the FEM integration guide](|page|/fem-integration.html) for use in a solver.

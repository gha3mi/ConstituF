---
title: Annulus benchmark (2D)
---

The annulus study evaluates constitutive responses and analytical parameter
sensitivities under prescribed mechanical and scalar fields. The annulus is a
sampling geometry, not a solution of global equilibrium or transport equations.
ConstituF evaluates each registered material using its declared field inputs.

See the [benchmark overview](|page|/benchmarks.html) and the
[annulus material results](|page|/benchmark-annulus-results.html).

## Geometry

The annulus has inner and outer radii \(a\) and \(b\):

\[
a=10\,\mathrm{mm},\qquad b=30\,\mathrm{mm},\qquad
N_r=24,\qquad N_\theta=64.
\]

Bulk responses are sampled over the annular area. Cohesive laws are
evaluated on a circular interface of radius \(R=20\,\mathrm{mm}\).

## Mechanical inputs

Prescribe the radial displacement

\[
u_r(r)=\alpha u_0\frac{a}{r},\qquad u_0=0.25\,\mathrm{mm},
\]

where \(\alpha\) is the loading amplitude. The nonzero polar
infinitesimal-strain components are

\[
\varepsilon_{rr}=\frac{\mathrm{d}u_r}{\mathrm{d}r}
=-\alpha u_0\frac{a}{r^2},\qquad
\varepsilon_{\theta\theta}=\frac{u_r}{r}
=\alpha u_0\frac{a}{r^2}.
\]

Small-strain laws receive the declared Mandel, plane or axisymmetric
components. Finite-strain laws receive
\(\boldsymbol F=\boldsymbol I+\nabla\boldsymbol u\), packed as a full tensor.

![Radial Cauchy stress from isotropic elasticity](|media|/benchmarks/annulus/isotropic_elastic_response.svg)

![Radial first Piola–Kirchhoff stress from Gent hyperelasticity](|media|/benchmarks/annulus/gent_response.svg)

## Thermal, electrical and transport inputs

For temperature, molar concentration, electric potential and pressure,
the prescribed scalar field \(\phi\) is

\[
\phi(r)=\phi_a+
\frac{\phi_b-\phi_a}{\ln(b/a)}\ln\left(\frac r a\right),
\qquad
\nabla\phi=
\frac{\phi_b-\phi_a}{r\ln(b/a)}\boldsymbol e_r.
\]

| Field | Inner value \(\phi_a\) | Outer value \(\phi_b\) |
| --- | ---: | ---: |
| Temperature | \(360\,\mathrm K\) | \(300\,\mathrm K\) |
| Molar concentration | \(1000\,\mathrm{mol/m^3}\) | \(500\,\mathrm{mol/m^3}\) |
| Electric potential | \(500\,\mathrm V\) | \(0\,\mathrm V\) |
| Pressure | \(2\,\mathrm{MPa}\) | \(1\,\mathrm{MPa}\) |

The gradient and, where declared, the scalar value are passed to the
material. Electrical laws use \(\boldsymbol E=-\nabla\phi\).
Coupled laws receive their declared mechanical and scalar fields together.

![Radial Fourier heat flux under the prescribed temperature gradient](|media|/benchmarks/annulus/fourier_isotropic_response.svg)

![Radial Darcy velocity under the prescribed pressure gradient](|media|/benchmarks/annulus/darcy_isotropic_response.svg)

## Loading history

History-dependent materials follow five accepted mechanical loading
amplitudes:

\[
\alpha=(0.25,0.65,1.00,0.40,0.85),\qquad
\Delta t=0.1\,\mathrm s.
\]

Each successful trial is committed before the next increment. The
figures show the last accepted response. Temperature, concentration,
electric potential and pressure remain fixed throughout the loading
history. Stateless laws use \(\alpha=1\).

![Radial J2 stress after the accepted loading history](|media|/benchmarks/annulus/j2_plasticity_response.svg)

## Cohesive interfaces

At \(R=20\,\mathrm{mm}\), prescribe the normal opening

\[
\delta_n(\theta)=\alpha(30\,\mu\mathrm m)
\left[1+0.15\cos(2\theta)\right].
\]

The mixed-mode law also receives a tangential opening of
\(0.2\delta_n\). Interface figures describe traction, not bulk stress
or damage.

![Normal traction along the mode-I cohesive interface](|media|/benchmarks/annulus/cohesive_mode_i_response.svg)

## Constitutive response

The plotted response is a selected radial component of Cauchy stress,
first Piola–Kirchhoff stress, a transport flux, electric response or
interface traction. The Cauchy component \(\sigma_{rr}\) and first Piola
component \(P_{rR}\) are different physical measures.

The linear colorbar reports the measure's physical units. Mechanical
stress figures use diverging `vik` with symmetric negative and positive
limits. Other quantities use sequential `magma`, `viridis` or `batlow`.
`response_manifest.csv` records actual data extrema and plotting limits.
Matching colors across models do not imply matching values or units.

These palettes are supplied by [ForColormap](https://github.com/vmagnin/forcolormap).
`vik` and `batlow` belong to the
[Scientific Colour Maps](https://www.fabiocrameri.ch/colourmaps/) collection.
`magma` and `viridis` originate from
[BIDS](https://github.com/BIDS/colormap).

## Analytical parameter sensitivity

For each declared scalar parameter \(p_j\), `material_t%eval` supplies
the analytical derivative of the selected physical response \(y\).
Prescribed fields and time increments remain fixed. For
history-dependent models, derivatives propagate through every accepted
constitutive increment. Define the scaled derivative

\[
q_j(r,\theta)=s_j
\left.\frac{\mathrm d y(r,\theta)}{\mathrm d p_j}
\right|_{\text{prescribed input history}}
\]

and the dimensionless relative RMS measure

\[
\eta_j=\frac{\operatorname{RMS}_{\Omega}(q_j)}
{\operatorname{RMS}_{\Omega}(y)}.
\]

The positive parameter reference scale is \(s_j=|p_j|\) when the
nominal parameter is nonzero. Otherwise, it is the largest nonzero
nominal magnitude among parameters with the same SI units in the benchmark
registry. The unit labels describe the chosen reference data.
At zero, this is a reference-scale derivative, not a logarithmic derivative.

For bulk materials, RMS uses annular area weights
\(\mathrm d\Omega=r\,\mathrm dr\,\mathrm d\theta\). For cohesive laws,
it uses interface arc length \(\mathrm d\ell=R\,\mathrm d\theta\).
The bulk heatmap entry for parameter \(j\) and radial station \(r_i\) is

\[
H_{j,i}=\log_{10}\left(
\frac{\sqrt{\langle q_j(r_i,\theta)^2\rangle_\theta}}
{\operatorname{RMS}_{\Omega}(y)}
\right).
\]

For cohesive interfaces, columns correspond to interface angles
\(\theta_i\) and the numerator is \(|q_j(\theta_i)|\). Heatmaps
encode **magnitude, not sign**. Signed scaled derivatives are exported
separately. Zero sensitivities, roundoff-level values and undefined
ratios are identified explicitly.

These are local constitutive derivatives, not global equilibrium
sensitivities.

## Reproduce and inspect

Generate both studies and their result pages:

```console
pixi run -e gnu benchmarks
pixi run docs
```

To run only the annulus:

```console
pixi run -e gnu annulus
```

The `benchmark_annulus` executable defaults to `all`, which generates
responses and analytical sensitivities for every registered material.
It also accepts `response` or `sensitivity` and an optional exact model slug:

```console
pixi run -e gnu fpm run --target benchmark_annulus --profile debug -- response
pixi run -e gnu fpm run --target benchmark_annulus --profile debug -- sensitivity gent
pixi run -e gnu fpm run --target benchmark_annulus --profile debug -- all gent
```

A model slug selects that material in any mode. Filtered runs preserve
the aggregate `response_manifest.csv` and `sensitivity_metrics.csv` reports
and the material-results gallery.

The outputs under
`docs/assets/benchmarks/annulus/` include
`<model>_response.svg`, `<model>_response.csv`,
`<model>_response.vtk`, `<model>_sensitivity.svg`,
`<model>_sensitivity.csv` and `<model>_sensitivity_fields.csv`.
The response and sensitivity reports are stored in the same directory.
The heatmap CSV contains normalized sensitivity magnitudes. The signed
scaled derivatives are in `<model>_sensitivity_fields.csv`.

Open the [annulus material results](|page|/benchmark-annulus-results.html)
to select a material and inspect its response contour, sensitivity
heatmap, reference parameters and exports. The gallery and material pages
are generated with the benchmarks and remain untracked, along with the
figures and numerical exports. See
[verification](|page|/verification.html) for independent material-law
checks and the [materials catalog](|page|/index.html) for field conventions
and parameter definitions.


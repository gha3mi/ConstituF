---
title: Spherical-cavity benchmark (3D)
---

The spherical-cavity study evaluates constitutive responses and analytical
parameter sensitivities under prescribed three-dimensional mechanical and
scalar fields. The cube and cavity are sampling geometry, not a solution of
global equilibrium or transport equations. Only compatible 3D bulk and
interface models are included. Plane and axisymmetric models remain in the
annulus study.

See the [benchmark overview](|page|/benchmarks.html) and the
[spherical-cavity material results](|page|/benchmark-cavity3d-results.html).

## Geometry

The outer cube occupies \([-L,L]^3\), with a concentric spherical cavity:

\[
L=30\,\mathrm{mm},\qquad a=12\,\mathrm{mm}.
\]

Six mapped cube-face patches, each with \(10\times10\) cells, form
12 conforming layers through the shell. The mesh contains 7,200
hexahedra, 600 spherical-boundary quadrilaterals and 7,826 distinct
nodes. The inner boundary approximates a sphere with planar facets.

For \(\boldsymbol q\) on the unit-cube surface, the mapping is

\[
\boldsymbol x(\boldsymbol q,\lambda)
=(1-\lambda)a\frac{\boldsymbol q}{\|\boldsymbol q\|}
+\lambda L\boldsymbol q,\qquad 0\leq\lambda\leq1.
\]

Bulk laws are sampled at hexahedral centers. Cohesive laws use
inner-boundary quadrilateral centers.

## Mechanical inputs

Prescribe an affine displacement with a cavity-centered radial
perturbation:

\[
\begin{aligned}
\boldsymbol u(\boldsymbol x,\alpha)
&=\alpha\left(\boldsymbol E_0\boldsymbol x
+\frac{u_0a^2}{r^3}\boldsymbol x\right),\\
r&=\|\boldsymbol x\|,\qquad u_0=30\,\mu\mathrm m,\\
\boldsymbol E_0&=
10^{-3}\begin{bmatrix}
1.2 & 0.35 & -0.20\\
0.35 & -0.50 & 0.25\\
-0.20 & 0.25 & 0.20
\end{bmatrix}.
\end{aligned}
\]

The exact gradient is

\[
\nabla\boldsymbol u=
\alpha\left[
\boldsymbol E_0+
\frac{u_0a^2}{r^3}
\left(\boldsymbol I-3\frac{\boldsymbol x\otimes\boldsymbol x}{r^2}\right)
\right].
\]

Small-strain laws receive \(\operatorname{sym}\nabla\boldsymbol u\)
in their declared three-dimensional Mandel representation.
Finite-strain laws receive
\(\boldsymbol F=\boldsymbol I+\nabla\boldsymbol u\).

## Thermal, electrical and transport inputs

Prescribe scalar fields with radial and Cartesian contributions:

\[
\begin{aligned}
T &=300+45\frac{a}{r}+15\frac{x}{L}
&&[\mathrm K],\\
c &=500+500\frac{a}{r}+80\frac{z}{L}
&&[\mathrm{mol\,m^{-3}}],\\
p &=10^6+10^6\frac{a}{r}+0.2\times10^6\frac{y}{L}
&&[\mathrm{Pa}],\\
\phi &=400\frac{a}{r}+70\frac{x}{L}
&&[\mathrm V].
\end{aligned}
\]

The exact spatial gradients provide thermal, concentration and
pressure-gradient inputs. Electrical laws use
\(\boldsymbol E=-\nabla\phi\). Coupled models receive all declared fields.

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

## Cohesive interfaces

On the spherical boundary, prescribe normal opening

\[
\delta_n=\alpha(30\,\mu\mathrm m)
\left[1+0.15\,\frac{3n_z^2-1}{2}\right],
\]

where \(n_z\) is the axial component of the cavity-radial unit vector.
The mixed-mode law additionally receives
\(\delta_{t1}=0.2\delta_n\) and \(\delta_{t2}=0\).
Cohesive results show normal traction, not a bulk stress field.

## Constitutive response

For mechanical laws, the plotted response is the cavity-radial
projection using \(\boldsymbol n=\boldsymbol x/r\):
\(\boldsymbol n\cdot\boldsymbol\sigma\boldsymbol n\) for Cauchy
stress or \(\boldsymbol n\cdot\boldsymbol P\boldsymbol n\) for
first Piola–Kirchhoff stress. These are different physical measures.
Vector fluxes and electric responses are projected along
\(\boldsymbol n\). Each figure identifies the selected measure and units.

Bulk VTK files contain the full uncut hexahedral mesh. For visibility,
bulk SVG previews omit the positive octant to expose the cavity and
internal faces. Interface VTK and SVG files contain the spherical
boundary only. The cutaway changes the preview, not the sampled
geometry or exported response. VTK displacement vectors are
prescribed, not globally solved. Plotting palettes follow the
[annulus conventions](|page|/benchmark-annulus.html).

## Analytical parameter sensitivity

For every declared scalar parameter \(p_j\), `material_t%eval` computes
the analytical derivative of the selected physical response \(y\)
at fixed prescribed input history. Derivatives of history-dependent
laws include all accepted increments. Define

\[
q_j=s_j\left.\frac{\mathrm d y}{\mathrm d p_j}
\right|_{\text{prescribed input history}},\qquad
\eta_j=\frac{\operatorname{RMS}_\Omega(q_j)}
{\operatorname{RMS}_\Omega(y)}.
\]

For a nonzero nominal parameter, \(s_j=|p_j|\). At zero,
\(s_j\) is the largest nonzero nominal magnitude among parameters
with the same SI units in the benchmark registry. These labels describe
the chosen reference data. Bulk RMS uses cell volumes.
Interface RMS uses facet areas.

Heatmap rows correspond to parameters and columns to 12 spatial bins.
For bulk models, bins are the radial mesh layers from cavity to
outer boundary. For interfaces, bins are grouped by
\(n_z=\cos\theta\). Each cell displays the base-10 logarithm of
the bin RMS of \(q_j\) divided by the global response RMS.
Zeros, roundoff-level values and undefined ratios are identified
explicitly. Heatmaps show **magnitude, not sign**. Signed scaled
derivatives are exported separately.

These are local constitutive derivatives, not global equilibrium
sensitivities.

## Reproduce and inspect

Generate both studies and their result pages:

```console
pixi run -e gnu benchmarks
pixi run docs
```

To run only the 3D study:

```console
pixi run -e gnu cavity-3d
```

The `benchmark_cavity3d` executable defaults to `all`, which generates
responses and analytical sensitivities for every compatible 3D material.
It also accepts `response` or `sensitivity` and an optional exact model slug:

```console
pixi run -e gnu fpm run --target benchmark_cavity3d --profile debug -- response
pixi run -e gnu fpm run --target benchmark_cavity3d --profile debug -- sensitivity gent
pixi run -e gnu fpm run --target benchmark_cavity3d --profile debug -- all gent
```

A model slug selects that material in any mode. Filtered runs preserve
the aggregate `response_manifest.csv` and `sensitivity_metrics.csv` reports
and the material-results gallery.

Under `docs/assets/benchmarks/cavity3d/`, compatible materials
produce `<model>_response.svg`, `<model>_response.csv`,
`<model>_response.vtk`, `<model>_sensitivity.svg`,
`<model>_sensitivity.csv` and `<model>_sensitivity_fields.csv`.
`response_manifest.csv` and `sensitivity_metrics.csv` summarize the
results. Open VTK volumes or interfaces with ParaView or VisIt to
inspect the uncut data.

Open the [spherical-cavity material results](|page|/benchmark-cavity3d-results.html)
to select a material and inspect its response preview, sensitivity
heatmap, reference parameters and exports. The gallery and material pages
are generated with the benchmarks and remain untracked, along with the
figures and numerical exports. See
[verification](|page|/verification.html) for independent material-law
checks.


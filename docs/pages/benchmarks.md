---
title: Benchmarks
---

ConstituF provides two prescribed-field studies for constitutive response and
analytical parameter sensitivity. Both evaluate compatible material models
through `material_t%eval` at prescribed sampling points. Neither solves a
global equilibrium or transport boundary-value problem. The plots show local
material responses to reproducible inputs, not experimental validation.

Both studies choose SI units for the values passed to the unitless API.
Their shared model registry owns the parameter and response unit labels.
Lengths shown in \(\mathrm{mm}\) in the formulations are represented in \(\mathrm{m}\) in the
numerical inputs. Figures, parameter tables and exported data retain their
unit labels so the results can be interpreted independently of the code.

## Studies

<div class="cf-benchmark-studies">
  <section class="cf-benchmark-study" aria-labelledby="cf-annulus-heading">
    <h3 id="cf-annulus-heading">Annulus (2D)</h3>
    <img src="../media/readme/annulus_isotropic_elastic.svg" width="1080" height="810" alt="Radial Cauchy stress from isotropic elasticity on the annulus">
    <div class="cf-benchmark-study-body">
      <p>Circular annulus with prescribed radial displacement and scalar fields. Includes plane, axisymmetric and interface models.</p>
      <p class="cf-benchmark-study-links"><a href="benchmark-annulus.html">Formulation</a> · <a href="benchmark-annulus-results.html">Material results</a></p>
    </div>
  </section>
  <section class="cf-benchmark-study" aria-labelledby="cf-cavity-heading">
    <h3 id="cf-cavity-heading">Spherical cavity (3D)</h3>
    <img src="../media/readme/cavity3d_isotropic_elastic.svg" width="1080" height="810" alt="Radial Cauchy stress from isotropic elasticity on the spherical-cavity cutaway">
    <div class="cf-benchmark-study-body">
      <p>Cube containing a spherical cavity, with prescribed three-dimensional fields. Includes compatible 3D bulk and interface models.</p>
      <p class="cf-benchmark-study-links"><a href="benchmark-cavity3d.html">Formulation</a> · <a href="benchmark-cavity3d-results.html">Material results</a></p>
    </div>
  </section>
</div>

Each formulation specifies the geometry, field measures, loading history,
sensitivity normalization and exports. Material result pages show response
figures, sensitivity heatmaps and numerical data.

## Reproduce the results

From the repository root, generate both studies before building FORD:

```console
pixi run -e gnu benchmarks
pixi run docs
```

For an individual study:

| Command | Produces |
| --- | --- |
| `pixi run -e gnu annulus` | Annulus responses and sensitivities |
| `pixi run -e gnu cavity-3d` | Spherical-cavity responses and sensitivities |

Annulus results are exported to `docs/assets/benchmarks/annulus/` and
spherical-cavity results to `docs/assets/benchmarks/cavity3d/`.

Each study has one fpm target: `benchmark_annulus` or `benchmark_cavity3d`.
Both default to `all` and accept `response` or `sensitivity` followed by an
optional exact model slug. See the commands in the
[annulus formulation](|page|/benchmark-annulus.html) and
[3D formulation](|page|/benchmark-cavity3d.html).

Generated material pages, CSV data, VTK fields and full-size SVG figures are
not committed. The two previews above are versioned examples. Build the
benchmark results before FORD to include the material result pages.

## Interpret the results

Response plots show a selected stress, flux, electric response or traction in
the indicated physical units. Heatmaps show the **magnitude** of the normalized
analytical derivative with respect to each declared scalar parameter. Signed
derivative fields are available in CSV files. Comparisons between different
models must account for their response measures, units and reference parameters.

Prescribed-field studies demonstrate repeatable material-point evaluation,
not independent verification of each law. See [verification](|page|/verification.html)
for material tests and [conventions](|page|/conventions.html) for field
representations.


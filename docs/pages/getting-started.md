---
title: Getting started
---

This tutorial evaluates an isotropic elastic material at one point, checks
the result against its analytical solution and requests a tangent and
parameter derivative. ConstituF computes local constitutive responses.
The calling application supplies the fields.

## Run the example

Install [Pixi](https://pixi.sh/) and clone the repository. From its root:

```console
pixi run -e gnu fpm run --example isotropic_elastic --profile debug
```

The source is [`example/isotropic_elastic.f90`](https://github.com/gha3mi/ConstituF/blob/main/example/isotropic_elastic.f90).
For all material examples, run `pixi run -e gnu example`.

## Evaluate a trial

The following self-contained program uses the same public interface.
This program chooses pascals for bulk modulus \(K\), shear modulus \(G\)
and stress. ConstituF accepts the numerical values directly. A caller may
choose another consistent unit system. The model has no internal state.

```fortran
program elastic_point
    use constituf, only: wp, material_t, cf_ok
    use constituf_isotropic_elastic, only: isotropic_elastic_t
    implicit none (type, external)

    type(material_t) :: material
    real(wp) :: strain(6,1), stress(6,1)
    real(wp) :: tangent(6,6,1), dsdp(6,1,1)
    integer :: status

    call material%init(isotropic_elastic_t(160.0e9_wp,80.0e9_wp),1,status, &
        sensitivity_parameters=[character(len=13) :: "shear_modulus"])
    if (status /= cf_ok) error stop "Initialization failed"

    strain = 0.0_wp
    strain(1,1) = 1.0e-3_wp
    call material%eval(strain,stress,tangent,status,dsdp=dsdp)
    if (status /= cf_ok) error stop "Evaluation failed"

    print '(a,es12.4)', "sigma_11 [Pa]: ",stress(1,1)
    print '(a,es12.4)', "d(sigma_11)/d(eps_11) [Pa]: ",tangent(1,1,1)
    print '(a,es12.4)', "d(sigma_11)/dG [-]: ",dsdp(1,1,1)

    call material%commit(status)
    if (status /= cf_ok) error stop "Commit failed"
end program elastic_point
```

The final index in every array identifies the material point, so
`(6,1)` is six components at **one point**, not an additional
tensor direction. A three-dimensional symmetric tensor \(\boldsymbol A\)
uses the Mandel vector
\([A_{11},A_{22},A_{33},\sqrt{2}A_{23},\sqrt{2}A_{13},\sqrt{2}A_{12}]^{\mathsf T}\).
The `tangent` array has shape
`(noutput,ninput,npoints)`. `dsdp` has
`(noutput,nselected,npoints)`.

## Check the result

The imposed strain has
\(\varepsilon_{11}=10^{-3}\) and
\(\varepsilon_{22}=\varepsilon_{33}=0\).
This is *uniaxial strain*, not a traction-free uniaxial-stress test.
For isotropic elasticity,

\[
\boldsymbol{\sigma}
=K\,\operatorname{tr}(\boldsymbol{\varepsilon})\boldsymbol{I}
+2G\,\operatorname{dev}(\boldsymbol{\varepsilon}).
\]

Consequently,

\[
\begin{aligned}
\sigma_{11}
&=\left(K+\frac43G\right)\varepsilon_{11}
 \approx 2.667\times10^8\ \mathrm{Pa},\\
\frac{\partial\sigma_{11}}{\partial\varepsilon_{11}}
&=K+\frac43G\approx2.667\times10^{11}\ \mathrm{Pa},\\
\left.\frac{\partial\sigma_{11}}{\partial G}\right|_{\boldsymbol\varepsilon}
&=\frac43\varepsilon_{11}\approx1.333\times10^{-3}.
\end{aligned}
\]

The last derivative is the response sensitivity to the selected scalar
parameter \(G\) at fixed prescribed strain. It is not the material tangent.
Lateral stresses do not vanish because the lateral strains are constrained.

## Continue with another material

A successful `eval` creates a trial. Call `commit` only after an
accepted solver step. `revert` discards a rejected trial. A stateless
material such as this one has no evolving history, but the same lifecycle
works for plasticity, damage and viscoelasticity.

Use the [materials catalog](|page|/index.html) to find a model and its
source equations. [Reading material examples](|page|/material-examples.html)
explains histories and multiple parameter derivatives.
[field conventions](|page|/conventions.html) defines tensor layouts and
units. For the available methods and utilities, see the
[public interface](|page|/public-interface.html).
For integration points inside an existing solver, continue with
[FEM integration](|page|/fem-integration.html).

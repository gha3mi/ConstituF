---
title: Calibrating material parameters
---

This how-to uses the optional
[[constituf_calibration(module)]] utility to fit a small set of
constitutive parameters to prescribed observations. ConstituF does
not fit a model automatically: the application must specify the
experiment, model constructor, bounds, units and a callback returning
predictions and **analytical parameter derivatives**.

## Run the supplied example

From the repository root:

```console
pixi run -e gnu fpm run --example calibrate_elastic --profile debug
```

The [Fortran example](https://github.com/gha3mi/ConstituF/blob/main/example/calibrate_elastic.f90)
generates synthetic three-point isotropic elastic stress data with
known moduli, then recovers bulk modulus \(K\) and shear modulus
\(G\). It is a repeatable algorithm check, **not** an experimental
material calibration.

## Set up the inverse problem

Provide observations `observed(nobs)` and strictly positive
`noise_scale(nobs)` in the **same physical units** as the corresponding
observations. The objective is dimensionless:

\[
\mathcal L(\boldsymbol{p})=\frac12\sum_i
\left[\frac{y_i(\boldsymbol{p})-y_i^{\mathrm{obs}}}{s_i}\right]^2,
\]

where \(s_i\) is the observation scale and \(p_j\) denotes a
selected scalar material parameter. Set finite, physically
admissible `lower_bound` and `upper_bound` arrays. Optional
`parameter_scale` values must have the same units as their parameters.

The optimizer calls
`evaluate(parameters,predicted,jacobian,status)`, where
`predicted(nobs)` and
`jacobian(nobs,nparameters)` use the observation order and
\(\partial y_i/\partial p_j\) at fixed prescribed experiment inputs.
For parameter scales \(a_j>0\) expressed in the units of \(p_j\),
the scaled optimization Jacobian is dimensionless:

\[
J_{ij}=\frac{a_j}{s_i}\frac{\partial y_i}{\partial p_j}.
\]

The parameter scales improve the conditioning of a fit involving
different physical parameter magnitudes.
For models supporting analytical sensitivities, initialize a fresh
`material_t` with selected `sensitivity_parameters` and pass
`dsdp` to every evaluation. History-dependent models must replay
**the entire accepted loading history for each parameter candidate**.
Do not reuse a committed state from a different candidate.

## Inspect the result

[[constituf_calibration(module)]] supplies a bounded,
damped Levenberg--Marquardt-style local least-squares solve.
Check `cf_ok` and inspect `fitted` and the dimensionless
`objective` before accepting a result. A failure may still return
the best accepted estimate, not a scientifically validated fit.
The routine does not estimate statistical confidence intervals,
guarantee uniqueness or replace experimental validation. Assess
identifiability, units, residuals and data uncertainty independently.

See [verification](|page|/verification.html) for tests of analytical
derivatives and [conventions](|page|/conventions.html) for the
distinction between material parameter derivatives and global
solver sensitivities.

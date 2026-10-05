# API reference

## Main functions

```@docs
ParameterEstimation.estimate
ParameterEstimation.check_identifiability
ParameterEstimation.sample_data
```

## Result types

```@docs
ParameterEstimation.EstimationResult
ParameterEstimation.IdentifiabilityData
ParameterEstimation.Interpolant
```

## Internals

These functions are used by `estimate`. They are not needed for normal use.

```@docs
ParameterEstimation.filter_solutions
ParameterEstimation.solve_ode
ParameterEstimation.solve_ode!
ParameterEstimation.interpolate
ParameterEstimation.rational_interpolation_coefficients
```

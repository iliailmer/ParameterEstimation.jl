# Parameter Estimation

## Introduction

In this tutorial, we provide a general overview of using `ParameterEstimation.jl`.

Assume we have a simple ODE model as below

```
``

```julia
using ParameterEstimation
using ModelingToolkit

# Input:
# -- Differential model
@parameters mu
@variables t x(t) y(t)
D = Differential(t)
@named Sigma = ODESystem([D(x) ~ -mu * x],
                         t, [x], [mu])
outs = [y ~ x^2 + x]

# -- Data
data = Dict(
  "t"     => [0.000, 0.333, 0.666, 1.000],
  x^2 + x => [2.000, 1.563, 1.229, 0.974])

# Run
res = estimate(Sigma, outs, data);
```
# Parameter Estimation

## Introduction

In this tutorial, we provide a general overview of using `ParameterEstimation.jl`.

Assume we have a simple ODE model with output as below

$$\begin{cases}\dot{x} = -\mu x,\\y = x^2+x\end{cases}$$

If we collect the sample at 4 time points between 0 and 1, we obtain a collection:

```
"t"     => [0.000, 0.333, 0.666, 1.000],
x^2 + x => [2.000, 1.563, 1.229, 0.974]
```

This is all that is needed for the program: a symbolic model (ODE and outputs) and a dictionary of data.

## Code
Below is the working code example:

```@example tutorial
using ParameterEstimation
using ModelingToolkit

# Input:
# -- MTK model
@parameters mu
@variables t x(t) y(t)
D = Differential(t)
@named Sigma = ODESystem([D(x) ~ -mu * x],
                         t, [x], [mu])
outs = [y ~ x^2 + x]

# -- Data
data = Dict(
  "t"     => [0.000, 0.333, 0.666, 1.000],
  x^2 + x => [2.000, 1.563, 1.229, 0.974]
)

# Run
res = estimate(Sigma, outs, data)
```
# Parameter Estimation

## Introduction

In this tutorial, we provide a general overview of using `ParameterEstimation.jl`.

Assume we have a simple ODE model with output as below

```math
\begin{cases}\dot{x} = -\mu x,\\y = x^2+x\end{cases}
```

If we collect the sample at 4 time points between 0 and 1, we obtain a collection:

```
"t"     => [0.000, 0.333, 0.666, 1.000],
  x^2 + x => [2.000, 1.563, 1.229, 0.974]
```

This is all that is needed for the program: a symbolic model (ODE and outputs) and a dictionary of data.

## Code

The code below defines the model and the data, and runs the estimation.

```@example tutorial
using ParameterEstimation
using ModelingToolkit

# Input:
# -- Differential model
@parameters mu
@independent_variables t
@variables x(t) y(t)
D = Differential(t)
@named Sigma = System([D(x) ~ -mu * x],
    t, [x], [mu])
outs = [y ~ x^2 + x]

# -- Data
data = Dict(
    "t" => [0.000, 0.333, 0.666, 1.000],
    x^2 + x => [2.000, 1.563, 1.229, 0.974])

# Run
res = estimate(Sigma, outs, data);
nothing # hide
```

The log shows the identifiability analysis: both `mu` and the initial value of `x` are
globally identifiable, so they can be estimated uniquely.

`estimate` returns a vector of [`ParameterEstimation.EstimationResult`](@ref), one for each
estimate that fits the data, best first:

```@example tutorial
for r in res
    println(r)
end
```

The estimates are sorted by `Error`, the mean absolute difference between the data and the
ODE solution with the estimated values. Use the first estimate. The exact values are
`mu = 0.5` and `x(0) = 1`; the data has 3 decimal places, so the best estimate is accurate
to about 3 digits.

The last estimate comes from a second solution of the equations: `x^2 + x = 2` holds for
both `x = 1` and `x = -2`. It fits these four rounded data points 15 times worse than the
best estimate. Always compare the errors when `estimate` returns more than one result.

The estimated values can be read by indexing with the parameter or the state:

```@example tutorial
best = res[1]
best[mu], best[x], best.err
```

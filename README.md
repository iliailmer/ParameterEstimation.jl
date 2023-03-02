# ParameterEstimation.jl


[![Tests](https://github.com/iliailmer/ParameterEstimation.jl/actions/workflows/tests.yml/badge.svg)](https://github.com/iliailmer/ParameterEstimation.jl/actions/workflows/tests.yml) [![Documentation](https://github.com/iliailmer/ParameterEstimation.jl/actions/workflows/Documentation.yml/badge.svg)](https://github.com/iliailmer/ParameterEstimation.jl/actions/workflows/Documentation.yml)
<!-- [![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle) -->

<p><a href="https://GitHub.com/iliailmer/ParameterEstimation.jl/releases/"><img alt="GitHub release" src="https://img.shields.io/github/release/iliailmer/ParameterEstimation.jl.svg"></a> <a href="https://GitHub.com/iliailmer/ParameterEstimation.jl/stargazers/"> <img alt="GitHub stars" src="https://img.shields.io/github/stars/iliailmer/ParameterEstimation.jl.svg?style=social&amp;label=Star&amp;maxAge=2592000"></a> </p>


Symbolic-Numeric package for parameter estimation in ODEs

## Installation

Currently is installable via


```julia

using Pkg
Pkg.add("ParameterEstimation.jl")
```

or

```julia

using Pkg
Pkg.add(url="https://github.com/iliailmer/ParameterEstimation.jl")
```

## Toy Example

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

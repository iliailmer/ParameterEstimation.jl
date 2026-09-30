# ParameterEstimation.jl

[![Tests](https://github.com/iliailmer/ParameterEstimation.jl/actions/workflows/tests.yml/badge.svg)](https://github.com/iliailmer/ParameterEstimation.jl/actions/workflows/tests.yml)
[![Documentation](https://github.com/iliailmer/ParameterEstimation.jl/actions/workflows/Documentation.yml/badge.svg)](https://github.com/iliailmer/ParameterEstimation.jl/actions/workflows/Documentation.yml)
[![GitHub release](https://img.shields.io/github/release/iliailmer/ParameterEstimation.jl.svg)](https://github.com/iliailmer/ParameterEstimation.jl/releases/)
[![Downloads](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv2%2Ftotal_downloads%2FParameterEstimation&query=%24.total_requests&label=downloads&color=blue)](https://juliapkgstats.com/pkg/ParameterEstimation)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)
[![GitHub stars](https://img.shields.io/github/stars/iliailmer/ParameterEstimation.jl.svg?style=social&label=Star&maxAge=2592000)](https://github.com/iliailmer/ParameterEstimation.jl/stargazers/)

Symbolic-Numeric package for parameter estimation in ODEs

## Installation

Install from the General registry:

```julia
using Pkg
Pkg.add("ParameterEstimation")
```

## Toy Example

```julia
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
```

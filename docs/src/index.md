# ParameterEstimation.jl

Software for Parameter Estimation Based on Identifiability Information and Real Data

## Installation

To install ParameterEstimation.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("ParameterEstimation")
```

## Citation

If you use ParameterEstimation.jl in your work, please cite:

O. Bassik, Y. Berman, S. Go, H. Hong, I. Ilmer, A. Ovchinnikov, C. Rackauckas, P. Soto, C. Yap.
Robust parameter estimation for rational ordinary differential equations.
*Applied Mathematics and Computation* 509 (2026), 129638.
[doi:10.1016/j.amc.2025.129638](https://doi.org/10.1016/j.amc.2025.129638)
(preprint: [arXiv:2303.02159](https://arxiv.org/abs/2303.02159))

```bibtex
@article{bassik2026robust,
  title   = {Robust parameter estimation for rational ordinary differential equations},
  author  = {Bassik, Oren and Berman, Yosef and Go, Soo and Hong, Hoon and Ilmer, Ilia and
             Ovchinnikov, Alexey and Rackauckas, Chris and Soto, Pedro and Yap, Chee},
  journal = {Applied Mathematics and Computation},
  volume  = {509},
  pages   = {129638},
  year    = {2026},
  doi     = {10.1016/j.amc.2025.129638}
}
```

## Feature Summary

  - Parameter estimation based on sample data
  - Estimated values are reported based on identifiability: local (finitely many), global (unique), unidentifiable.

"""
	Interpolant

A structure that stores information about the interpolation result.

# Fields
- `I`: the single variable differentiable function, result of interpolation.
- `dIdt`: the `TaylorSeries` derivative of `I`.
"""
struct Interpolant
	f::Any
end

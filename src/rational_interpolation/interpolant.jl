"""
    Interpolant

The result of interpolating the data of one measured quantity.

# Fields
- `f`: the interpolating function of one variable (time). Its derivatives are computed by automatic differentiation.
"""
struct Interpolant
    f::Any
end

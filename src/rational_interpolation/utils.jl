"""
    function eval_derivs(polynomial_system, interpolant::Interpolant,
                        y_function_name,
                        identifiability_result;
                        method = :homotopy)

This function evaluates the derivatives of the interpolated function `y_function_name` using the `interpolant` object.
    Derivatives are substituted into the polynomial system.
"""
function eval_derivs(polynomial_system, interpolant::Interpolant,
                     y_function_name,
                     identifiability_result;
                     method = :homotopy)
    if isequal(method, :homotopy)
        for (y_func, y_deriv_order) in pairs(identifiability_result["Y_eq"])
            if occursin(y_function_name, string(y_func))
                y_derivs_vals = Dict(ParameterEstimation.nemo2hc(y_func) => interpolant.dIdt[y_deriv_order] *
                                                                            factorial(y_deriv_order))
                polynomial_system = HomotopyContinuation.evaluate(ParameterEstimation.nemo2hc.(polynomial_system),
                                                                  y_derivs_vals)
            end
        end
    elseif isequal(method, :msolve)
        y_derivs = Vector{SIAN.Nemo.fmpq_mpoly}()
        y_vals = Vector{SIAN.Nemo.fmpq}()
        for (y_func, y_deriv_order) in pairs(identifiability_result["Y_eq"])
            if occursin(y_function_name, string(y_func))
                push!(y_derivs, y_func)
                push!(y_vals,
                      rationalize(Float64(interpolant.dIdt[y_deriv_order] *
                                          factorial(y_deriv_order))))
            end
        end
        polynomial_system = [SIAN.Nemo.evaluate(poly, y_derivs, y_vals)
                             for poly in polynomial_system]
    end
    return polynomial_system
end

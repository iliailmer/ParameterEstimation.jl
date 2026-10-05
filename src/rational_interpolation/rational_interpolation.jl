
"""
    interpolate(identifiability_result, data_sample, measured_quantities, inputs;
                interpolator, diff_order::Int = 1, at_t::Float = 0.0,
                method::Symbol = :homotopy)

Interpolate the data of each measured quantity and substitute the derivatives of the
interpolants at `at_t` into the polynomial system of `identifiability_result`.

# Arguments
- `identifiability_result`: the result of `check_identifiability`.
- `data_sample`: the data, see [`estimate`](@ref).
- `measured_quantities`: the measured quantities (outputs as equations of the form `y ~ x`).
- `inputs`: known input functions, as equations.
- `interpolator`: a pair `name => function`, for example `"AAA" => ParameterEstimation.aaad`.
- `diff_order::Int = 1`: the highest derivative order to compute.
- `at_t::Float = 0.0`: the time point at which the derivatives are computed.
- `method::Symbol = :homotopy`: the polynomial system solver the system is prepared for.

# Returns
- `interpolants`: a dictionary from each measured quantity to its [`ParameterEstimation.Interpolant`](@ref).
- `polynomial_system`: the polynomial system with the derivative values substituted, a `HomotopyContinuation.System` for `method = :homotopy`.
"""
function interpolate(identifiability_result, data_sample,
        measured_quantities, inputs; interpolator,
        diff_order::Int = 1, at_t::Float = 0.0,   #TODO(orebas)should we remove diff_order?
        method::Symbol = :homotopy)
    polynomial_system = identifiability_result["polynomial_system"]
    interpolants = Dict{Any, Interpolant}()
    sampling_times = data_sample["t"]
    for (key, sample) in pairs(data_sample)
        if key == "t"
            continue
        end
        y_function_name = map(x -> replace(string(x.lhs), "(t)" => ""),
            filter(x -> string(x.rhs) == string(key),
                measured_quantities))[1]
        interpolant = ParameterEstimation.interpolate(sampling_times, sample,
            interpolator,
            diff_order)
        interpolants[key] = interpolant
        err = sum(abs.(sample - interpolant.f.(sampling_times))) / length(sampling_times)
        @debug "Mean Absolute error in interpolation: $err interpolating $key"
        polynomial_system = eval_derivs(polynomial_system, interpolant, y_function_name,
            inputs, identifiability_result, at_time = at_t, method = method)
    end
    if isequal(method, :homotopy)
        try
            polynomial_system = HomotopyContinuation.System(polynomial_system)
        catch KeyError
            throw(ArgumentError("HomotopyContinuation threw a KeyError, it is possible that " *
                                "you are using Unicode characters in your input. Consider " *
                                "using ASCII characters instead."))
        end
    end
    identifiability_result["polynomial_system_to_solve"] = polynomial_system
    return interpolants, polynomial_system
end

"""
	interpolate(time, sample, numer_degree::Int, diff_order::Int = 1, at_t::Float = 0.0)

This function performs a rational interpolation of the data `sample` at the points `time` with numerator degree `numer_degree`.
It returns an `Interpolant` object that contains the interpolated function and its derivatives.
"""
function interpolate(time, sample, interpolator, diff_order::Int = 1)
    interpolated_function = ((interpolator.second))(time, sample)
    return Interpolant(interpolated_function)
end

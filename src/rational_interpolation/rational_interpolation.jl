

"""
	interpolate(identifiability_result, data_sample,
				measured_quantities; interpolation_degree::Int = 1,
				diff_order::Int = 1, at_t::Float = 0.0,
				method::Symbol = :homotopy)

This function performs the key step in parameter estimation.

	It interpolates the data in `data_sample` and computes the `TaylorSeries` expansion.
	These results are stored in the `Interpolant` object and are applied to the polynomial system in `identifiability_result`.

# Arguments
- `identifiability_result`: the result of the identifiability check.
- `data_sample`: a dictionary of the data samples. The keys are the symbols of the measured quantities and the values are the data samples.
- `measured_quantities`: the measured quantities (outputs as equations of the form `y ~ x`).
- `interpolation_degree::Int = 1`: the degree of the numerator of the rational interpolation.
- `diff_order::Int = 1`: the order of the derivative to be computed.
- `at_t::Float = 0.0`: the time point where the Taylor series expansion is computed.
- `method::Symbol = :homotopy`: the method used to solve the polynomial system. Can be one of :homotopy (recommended) and :msolve.

# Returns
- `System`: the polynomial system with the interpolated data applied. This system is compatible with `HomotopyContinuation` solving.
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
			identifiability_result["polynomial_system_to_solve"] = HomotopyContinuation.System(polynomial_system)
		catch KeyError
			throw(ArgumentError("HomotopyContinuation threw a KeyError, it is possible that " *
								"you are using Unicode characters in your input. Consider " *
								"using ASCII characters instead."))
		end
	else
		identifiability_result["polynomial_system_to_solve"] = polynomial_system
	end
	return interpolants
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


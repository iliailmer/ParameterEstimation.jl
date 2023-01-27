"""
    rational_interpolation_coefficients(x, y, n)

Perform a rational interpolation of the data `y` at the points `x` with numerator degree `n`.
This function only returns the coefficients of the numerator and denominator polynomials.

# Arguments
- `x`: the points where the data is sampled (e.g. time points).
- `y`: the data sample.
- `n`: the degree of the numerator.

# Returns
- `c`: the coefficients of the numerator polynomial.
- `d`: the coefficients of the denominator polynomial.
"""
function rational_interpolation_coefficients(x, y, n)
    N = length(x)
    m = N - n - 1
    A = zeros(N, N)
    if m > 0
        A_left_submatrix = reduce(hcat, [x .^ (i) for i in 0:(n)])
        A_right_submatrix = reduce(hcat, [x .^ (i) for i in 0:(m - 1)])
        A = hcat(A_left_submatrix, -y .* A_right_submatrix)
        b = y .* (x .^ m)
        try
            prob = LinearSolve.LinearProblem(A, b)
            c = LinearSolve.solve(prob)
            return c[1:(n + 1)], [c[(n + 2):end]; 1]
        catch SingularException
            lu_res = lu(A)
            y = lu_res.L \ lu_res.P * b
            c = lu_res.U \ y
            return c[1:(n + 1)], [c[(n + 2):end]; 1]
        end

    else
        A = reduce(hcat, [x .^ i for i in 0:n])
        b = y
        prob = LinearSolve.LinearProblem(A, b)
        c = LinearSolve.solve(prob)
        return c, [1]
    end
end

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
                     measured_quantities; interpolation_degree::Int = 1,
                     diff_order::Int = 1, at_t::Float = 0.0,
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
                                                      interpolation_degree,
                                                      diff_order, at_t)
        interpolants[key] = interpolant
        err = sum(abs.(sample - interpolant.I.(sampling_times))) / length(sampling_times)
        @debug "Mean Absolute error in interpolation: $err interpolating $key"
        polynomial_system = eval_derivs(polynomial_system, interpolant, y_function_name,
                                        identifiability_result, method = method)
    end
    if isequal(method, :homotopy)
        try
            identifiability_result["polynomial_system_to_solve"] = HomotopyContinuation.System(polynomial_system)
        catch KeyError
            throw(ArgumentError("HomotopyContinuation threw a KeyError, it is likely that " *
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
function interpolate(time, sample, numer_degree::Int, diff_order::Int = 1,
                     at_t::Float = 0.0)
    numer_coef, denom_coef = rational_interpolation_coefficients(time, sample,
                                                                 numer_degree)
    numer_function(t) = sum(numer_coef[i] * t^(i - 1) for i in eachindex(numer_coef))
    denom_function(t) = sum(denom_coef[i] * t^(i - 1) for i in eachindex(denom_coef))
    interpolated_function(t) = numer_function(t) / denom_function(t)
    return Interpolant(interpolated_function,
                       differentiate_interpolated(interpolated_function, diff_order, at_t))
end

function differentiate_interpolated(interpolated_function, diff_order::Int,
                                    at_t::Float = 0.0)
    τ = Taylor1(diff_order + 1)
    taylor_expantion = interpolated_function(τ - at_t)
    return taylor_expantion
end
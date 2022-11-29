function rational_interpolation_coefficients(x, y, n)
    # x and y are vectors of length N
    # n is numerator degree
    # returns the coefficients of the numerator and denominator polynomials
    N = length(x)
    m = N - n - 1
    A = zeros(N, N)
    if m > 0
        A_left_submatrix = reduce(hcat, [x .^ (i) for i in 0:(n)])
        A_right_submatrix = reduce(hcat, [x .^ (i) for i in 0:(m - 1)])
        A = hcat(A_left_submatrix, -y .* A_right_submatrix)
        b = y .* (x .^ m)
        prob = LinearSolve.LinearProblem(A, b)
        c = LinearSolve.solve(prob)
        return c[1:(n + 1)], [c[(n + 2):end]; 1]
    else
        A = reduce(hcat, [x .^ i for i in 0:n])
        b = y
        prob = LinearSolve.LinearProblem(A, b)
        c = LinearSolve.solve(prob)
        return c, [1]
    end
end

function interpolate(identifiability_result, data_sample, time_interval,
                     measured_quantities,
                     interpolation_degree::Int = 1,
                     diff_order::Int = 1)
    polynomial_system = identifiability_result["polynomial_system"]
    for (key, sample) in pairs(data_sample)
        y_function_name = map(x -> replace(string(x.lhs), "(t)" => ""),
                              filter(x -> string(x.rhs) == string(key),
                                     measured_quantities))[1]
        tsteps = range(time_interval[1], time_interval[2], length = length(sample))
        interpolant = ParameterEstimation.interpolate(tsteps, sample,
                                                      interpolation_degree,
                                                      diff_order)
        err = sum(abs.(sample - interpolant.I.(tsteps))) / length(tsteps)
        @info "Mean Absolute error in interpolation: $err interpolating $key"
        for (y_func, y_deriv_order) in pairs(identifiability_result["Y_eq"])
            if occursin(y_function_name, string(y_func))
                y_derivs_vals = Dict(ParameterEstimation.nemo2hc(y_func) => interpolant.dIdt[y_deriv_order] *
                                                                            factorial(y_deriv_order))
                polynomial_system = HomotopyContinuation.evaluate(ParameterEstimation.nemo2hc.(polynomial_system),
                                                                  y_derivs_vals)
            end
        end
    end
    return System(polynomial_system)
end

function interpolate(time, sample, numer_degree::Int, diff_order::Int = 1)
    numer_coef, denom_coef = rational_interpolation_coefficients(time, sample,
                                                                 numer_degree)
    numer_function(t) = sum(numer_coef[i] * t^(i - 1) for i in 1:length(numer_coef))
    denom_function(t) = sum(denom_coef[i] * t^(i - 1) for i in 1:length(denom_coef))
    interpolated_function(t) = numer_function(t) / denom_function(t)
    return Interpolant(interpolated_function,
                       differentiate_interpolated(interpolated_function, diff_order))
end

function differentiate_interpolated(interpolated_function, diff_order::Int)
    τ = Taylor1(diff_order + 1)
    taylor_expantion = interpolated_function(τ)
    return taylor_expantion
end
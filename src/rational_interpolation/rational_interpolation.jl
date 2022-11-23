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

function interpolate(time, sample, numer_degree::Int)
    numer_coef, denom_coef = rational_interpolation_coefficients(time, sample,
                                                                 numer_degree)
    numer_function(t) = sum(numer_coef[i] * t^(i - 1) for i in 1:length(numer_coef))
    denom_function(t) = sum(denom_coef[i] * t^(i - 1) for i in 1:length(denom_coef))
    interpolated_function(t) = numer_function(t) / denom_function(t)
    return interpolated_function
end
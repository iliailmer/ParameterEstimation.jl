function rational_interpolation_coefficients(x, y, n)
    # x and y are vectors of length N
    # n is numerator degree
    # returns the coefficients of the numerator and denominator polynomials
    N = length(x)
    m = N - n - 1
    A = zeros(N, n + m + 1)
    if m > 0
        A_left_submatrix = reduce(hcat, [(x .^ i) for i in 0:m-1])
        A_right_submatrix = reduce(hcat, [x .^ i for i in 0:n])
        A = hcat(y .* A_left_submatrix, A_right_submatrix)
        b = -y .* x .^ m
        prob = LinearProblem(A, b)
        c = solve(prob)
        return -c[m+1:end], push!(c[1:m], 1)
    else
        A = reduce(hcat, [x .^ i for i in 0:n])
        b = y
        prob = LinearProblem(A, b)
        c = solve(prob)
        return c, []
    end
end


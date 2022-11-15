using LinearSolve
using Plots
using Statistics

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

f(t) = sin(t)
N = 50
x = vcat(-10:(20/N):10)
y = f.(x)

errors = []
min_error = Inf
best_g = nothing
best_n = 0
plot(x, y)
for n ∈ 1:N-1
    a, b = rational_interpolation_coefficients(x, y, n)
    numerator(t) = sum(a[i] * t^(i - 1) for i in 1:length(a))
    if length(d) > 0
        denominator(t) = sum(b[i] * t^(i - 1) for i in 1:length(b))
    else
        denominator(t) = 1
    end
    g(t) = numerator(t) / denominator(t)
    err = maximum(abs.(f.(x) - g.(x)))
    push!(errors, err)
    if err < min_error
        min_error = err
        best_g = g
        best_n = n
    end
end
plot!(x, best_g.(x))

plot(errors)


println("Best n: $best_n")
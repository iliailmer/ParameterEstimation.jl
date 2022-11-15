# using Statistics, Plots

# f(t) = sin(t)
# N = 50
# x = vcat(-10:(20/N):10)
# y = f.(x)

# errors = []
# min_error = Inf
# best_g = nothing
# best_n = 0
# plot(x, y)
# for n ∈ 1:N-1
#     a, b = rational_interpolation_coefficients(x, y, n)
#     numerator(t) = sum(a[i] * t^(i - 1) for i in 1:length(a))
#     if length(d) > 0
#         denominator(t) = sum(b[i] * t^(i - 1) for i in 1:length(b))
#     else
#         denominator(t) = 1
#     end
#     g(t) = numerator(t) / denominator(t)
#     err = maximum(abs.(f.(x) - g.(x)))
#     push!(errors, err)
#     if err < min_error
#         min_error = err
#         best_g = g
#         best_n = n
#     end
# end
# plot!(x, best_g.(x))

# plot(errors)


# println("Best n: $best_n")
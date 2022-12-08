"""
    mean_abs_err(y_true, y_sample)

Compute the mean absolute error between `y_true` and `y_sample`.
"""
function mean_abs_err(y_true, y_sample)
    @assert length(y_true)==length(y_sample) "Unequal sample sizes: $(length(y_true)) and $(length(y_sample))"
    return sum(abs.(y_true .- y_sample)) / length(y_sample)
end
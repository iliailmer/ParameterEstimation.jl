function mean_abs_err(y_true, y_sample)
    @assert length(y_true)==length(y_sample) "Unequal sample sizes"
    return sum(abs.(y_true .- y_sample)) / length(y_sample)
end
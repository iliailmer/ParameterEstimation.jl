function differentiate_interpolated(interpolated_function, diff_order::Int)
    τ = Taylor1(diff_order + 1)
    taylor_expantion = interpolated_function(τ)
    return taylor_expantion
end

function squarify_system(poly_system::Vector{Expression})
    indets = HomotopyContinuation.variables(poly_system)
    M = randn(1, length(poly_system) - length(indets) + 1)
    return vcat(poly_system[1:(length(indets) - 1)], M * poly_system[length(indets):end])
end
function differentiate_interpolated(interpolated_function, diff_order::Int)
    τ = Taylor1(diff_order + 1)
    taylor_expation = interpolated_function(τ)
    return taylor_expantion
end

function construct_poly_equations(identifiability_result::Dict, interpolated_function,
                                  diff_order::Int)
    taylor_expantion = differentiate_interpolated(interpolated_function, diff_order)
    Et = identifiability_result["Et"]
    y_eq = OrderedDict()
    for (key, value) in pairs(identifiability_result["Y_eq"])
        println(value)
        y_eq[HomotopyContinuation.variables(Symbol(key))[1]] = taylor_expantion[Int(value)]
    end

    Et_hc = [HomotopyContinuation.subs(ParameterEstimation.nemo2hc(each), pairs(y_eq)...)
             for each in Et]
    return Et_hc
end

function squarify_system(poly_system::Vector{Expression})
    indets = HomotopyContinuation.variables(poly_system)
    M = randn(1, length(poly_system) - length(indets) + 1)
    return vcat(poly_system[1:(length(indets) - 1)], M * poly_system[length(indets):end])
end
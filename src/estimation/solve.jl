function solve_via_homotopy(identifiability_result, model; real_tol = 1e-10)
    parameters = ModelingToolkit.parameters(model)
    polynomial_system = identifiability_result["polynomial_system_to_solve"]

    @info "Solving $(length(polynomial_system)) polynomial equations in $(length(polynomial_system.variables)) variables"
    states = ModelingToolkit.states(model)
    results = HomotopyContinuation.solve(polynomial_system; show_progress = false)
    all_solutions = HomotopyContinuation.real_solutions(results)
    if length(all_solutions) == 0
        all_solutions = HomotopyContinuation.solutions(results)
        if length(all_solutions) == 0
            @debug "Interpolation numerator degree $(interpolation_degree): No solutions found"
            return Vector{EstimationResult}()
        end
    end
    all_solutions_ = Vector{Dict}([])
    state_param_map = merge(Dict(replace(string(x), "(t)" => "") => x for x in states),
                            Dict(string(x) => x for x in parameters))
    for sol in all_solutions
        tmp = Dict()
        sol = map(each -> to_exact(each; tol = real_tol), sol)
        for (idx, v) in enumerate(polynomial_system.variables)
            if endswith(string(v), "_0")
                tmp[state_param_map[string(v)[1:(end - 2)]]] = sol[idx]
            end
        end
        for (key, val) in identifiability_result.transcendence_basis_subs
            if endswith(string(key), "_0")
                tmp[state_param_map[string(key)[1:(end - 2)]]] = Int(Meta.parse("$val"))
            end
        end
        push!(all_solutions_, tmp)
    end
    return all_solutions_
end

function solve_via_msolve(identifiability_result, model; real_tol = 1e-10)
    polynomial_system = identifiability_result["polynomial_system_to_solve"]
    all_vars = reduce((x, y) -> union(x, y), vars(poly) for poly in polynomial_system)
    R, _ = Oscar.PolynomialRing(QQ, string.(all_vars))
    ps = [change_ring(poly, R) for poly in polynomial_system]
    i = Oscar.ideal(R, ps)
    solution = Oscar.msolve(i)
    throw(ArgumentError("Function solve_via_msolve is not implemented yet"))
end
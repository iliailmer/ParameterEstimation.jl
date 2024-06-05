function solve_via_homotopy(identifiability_result, model; real_tol = 1e-12)
	@debug "Solving $(length(identifiability_result["polynomial_system_to_solve"])) polynomial equations in $(length(identifiability_result["polynomial_system_to_solve"].variables)) variables"

	polynomial_system = identifiability_result["polynomial_system_to_solve"]
	state_param_map = merge(Dict(replace(string(x), "(t)" => "") => x
								 for x in ModelingToolkit.unknowns(model)),
		Dict(string(x) => x for x in ModelingToolkit.parameters(model)))
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
	for sol in all_solutions
		tmp = Dict()
		sol = map(each -> to_exact(each; tol = real_tol), sol)
		for (idx, v) in enumerate(polynomial_system.variables)
			if endswith(string(v), "_0")
				tmp[state_param_map[string(v)[1:(end-2)]]] = sol[idx]
			end
		end
		for (key, val) in identifiability_result.transcendence_basis_subs
			if endswith(string(key), "_0")
				tmp[state_param_map[string(key)[1:(end-2)]]] = Int(Meta.parse("$val"))
			end
		end
		push!(all_solutions_, tmp)
	end
	return all_solutions_
end

function solve_via_msolve(identifiability_result, model; real_tol = 1e-12)
    throw("Not implemented, sorry.")
    # polynomial_system = identifiability_result["polynomial_system_to_solve"]
    # state_param_map = merge(Dict(replace(string(x), "(t)" => "") => x
    # 							 for x in ModelingToolkit.states(model)),
    # 	Dict(string(x) => x for x in ModelingToolkit.parameters(model)))

    # all_vars = reduce((x, y) -> union(x, y), vars(poly) for poly in polynomial_system)
    # R, _ = Nemo.polynomial_ring(QQ, string.(all_vars); ordering = :degrevlex)
    # ps = [SIAN.parent_ring_change(poly, R) for poly in polynomial_system]
    # i = Oscar.ideal(R, ps)
    # all_solutions_ = Vector{Dict}([])

    # if Oscar.VERSION_NUMBER == v"0.11.3"
    # 	solutions, rat_param = Oscar.real_solutions(i)
    # elseif Oscar.VERSION_NUMBER == v"0.10.0"
    # 	rat_param, solutions = Oscar.msolve(i)
    # else
    # 	solutions, rat_param = Oscar.real_solutions(i)
    # end
    # for sol in solutions
    # 	tmp = Dict()
    # 	sol = map(each -> Float64(each), sol)
    # 	for (idx, v) in enumerate(all_vars)
    # 		if endswith(string(v), "_0")
    # 			tmp[state_param_map[string(v)[1:(end-2)]]] = sol[idx]
    # 		end
    # 	end
    # 	for (key, val) in identifiability_result.transcendence_basis_subs
    # 		if endswith(string(key), "_0")
    # 			tmp[state_param_map[string(key)[1:(end-2)]]] = Int(Meta.parse("$val"))
    # 		end
    # 	end
    # 	push!(all_solutions_, tmp)
    # end
    # return all_solutions_
end

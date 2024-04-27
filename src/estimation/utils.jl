function check_constraints(estimate, parameter_constraints, ic_constraints)
	sat = true
	if (!isnothing(parameter_constraints))
		#		println("Here: ", typeof(estimate))

		#		println("Here: ", estimate)

		#		println("Here: ", estimate.parameters)
		for (k, v) in pairs(estimate.parameters)
			if (haskey(parameter_constraints, k))
				if !(parameter_constraints[k][1] <= v && v <= parameter_constraints[k][2])
					sat = false
				end
			end
		end
	end

	if (!isnothing(ic_constraints))
		for (k, v) in pairs(estimate.states)
			if (haskey(ic_constraints, k))
				if !(ic_constraints[k][1] <= v && v <= ic_constraints[k][2])
					sat = false
				end
			end
		end
	end
	return sat
end


function compare_estimation_results(a, b)
	err = 0.0
	akeys = keys(a.states)
	bkeys = keys(b.states)
	for i in akeys
		if i in bkeys
			err = max(err, abs(a.states[i] - b.states[i]))
		else
			return 1.0e10
		end
	end
	akeysp = keys(a.parameters)
	bkeysp = keys(b.parameters)
	for j in akeysp
		if j in bkeysp
			err = max(err, abs(a.parameters[j] - b.parameters[j]))
		else
			return 1.0e10
		end
	end
	return err
end

function new_clustering(estimates, tol = 1e-3)
	new_estimates = []
	for e in estimates
		keep = true
		for f in new_estimates
			if (((abs(e.err - f.err)) < tol) && (compare_estimation_results(e, f) < tol))
				keep = false
			end
		end
		if (keep)
			push!(new_estimates, e)
		end

	end
	return new_estimates
end

function post_process(estimates, filtermode = :new, parameter_constraints = nothing, ic_constraints = nothing; threaded=false)
	if threaded
		estimates = filter(x -> !isnothing(x), estimates)
		estimates = vcat(estimates...)
	end
	#filter out the empty vectors
	estimates = filter(x -> length(x) > 0, estimates)
	estimates = filter(x -> x[1].return_code == ReturnCode.Success, estimates)

	if (filtermode == :new)
		estimates = collect(Iterators.flatten(estimates))
		estimates = filter(x -> check_constraints(x, parameter_constraints, ic_constraints), estimates)
		estimates = sort(estimates, by = x -> x.err)

		estimates_clustered_by_params = new_clustering(estimates, 1e-2)  #TODO(orebas): this is a magic number.
		estimates_clustered_by_params = filter(x -> (x.err < 1.0), estimates_clustered_by_params)  #also a magic number
		return estimates_clustered_by_params
	else
		best_solution = nothing
		for each in estimates
			if all(!isnothing(x.err) for x in each)
				if best_solution === nothing
					best_solution = each
				else
					if sum(x.err for x in each) < sum(x.err for x in best_solution)
						best_solution = each
					end
				end
			else
				best_solution = nothing
			end
		end
		if best_solution !== nothing
			@info "Best estimate found"
			return best_solution
		else
			@warn "No solution found"
			return estimates
		end
	end
end

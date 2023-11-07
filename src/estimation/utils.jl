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


function post_process(estimates, filtermode = :new, parameter_constraints = nothing, ic_constraints = nothing)
	if Threads.nthreads() > 1
		estimates = filter(x -> !isnothing(x), estimates)
		estimates = vcat(estimates...)
	end
	#filter out the empty vectors
	estimates = filter(x -> length(x) > 0, estimates)
	estimates = filter(x -> x[1].return_code == ReturnCode.Success, estimates)

	if (filtermode == :new)
		estimates = collect(Iterators.flatten(estimates))
		estimates = filter(x -> check_constraints(x, parameter_constraints, ic_constraints), estimates)
		return estimates
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

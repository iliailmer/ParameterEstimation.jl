function post_process(estimates, filtermode = :new)
	if Threads.nthreads() > 1
		estimates = filter(x -> !isnothing(x), estimates)
		estimates = vcat(estimates...)
	end
	#filter out the empty vectors
	estimates = filter(x -> length(x) > 0, estimates)
	estimates = filter(x -> x[1].return_code == ReturnCode.Success, estimates)

	if (filtermode == :new)
		return collect(Iterators.flatten(estimates))
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

function estimate_threaded(model::ModelingToolkit.ODESystem,
	 measured_quantities::Vector{ModelingToolkit.Equation},
	 inputs::Vector{ModelingToolkit.Equation},
	 data_sample::AbstractDict{Any, Vector{T}} = Dict{Any, Vector{T}}();
	 at_time::T, solver = Vern9(), interpolators = nothing, report_time = minimum(data_sample["t"]),
	 method = :homotopy, real_tol::Float64 = 1e-14,filtermode = :new, parameter_constraints = nothing, ic_constraints = nothing) where {T <: Float}
	@info "Using threaded estimation."
	check_inputs(measured_quantities, data_sample)
	datasize = length(first(values(data_sample)))
	
    if isnothing(interpolators)
        interpolators = default_interpolator(datasize)
    end

	id = ParameterEstimation.check_identifiability(model;
	measured_quantities = measured_quantities,
	inputs = [Num(each.lhs)
			  for each in inputs])
	@info "Estimating via the interpolators: $(keys(interpolators))"
	interpolators = collect(interpolators)
	estimates = Vector{Vector{EstimationResult}}(undef, length(interpolators))
	Threads.@threads for t in eachindex(interpolators)
		interp = interpolators[t]
		unfiltered = estimate_single_interpolator(model, measured_quantities, inputs, data_sample;
			identifiability_result = id,
			interpolator = interp, at_time = at_time, report_time,
			method = method,
			real_tol = real_tol)
		if length(unfiltered) > 0
			estimates[t] = filter_solutions(unfiltered, id, model, inputs,
				data_sample; solver = solver, filtermode)
		else
			estimates[t] = [
				EstimationResult(model, Dict(), interp, at_time,
					Dict{Any, Interpolant}(), ReturnCode.Failure,
					datasize, report_time),
			]
		end
	end
	return post_process(id, estimates, filtermode, parameter_constraints, ic_constraints)
end

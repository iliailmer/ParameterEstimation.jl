
function backsolve_initial_conditions(model, E, report_time, inputs::Vector{Equation}, data_sample;
	solver = Tsit5(), abstol = 1e-12, reltol = 1e-12)
	initial_conditions = [E[s] for s in ModelingToolkit.states(model)]
	parameter_values = [E[p] for p in ModelingToolkit.parameters(model)]
	tspan = (E.at_time, report_time)  #this is almost always backwards!

	ode_equations = ModelingToolkit.equations(model)
	ode_equations = substitute(ode_equations,
		Dict(each.lhs => Num(each.rhs) for each in inputs))
	t = ModelingToolkit.get_iv(model)
	@named new_model = ODESystem(ode_equations, t, ModelingToolkit.states(model),
		ModelingToolkit.parameters(model))
	prob = ODEProblem(new_model, initial_conditions, tspan, parameter_values)
	saveat = range(tspan[1], tspan[2], length = length(data_sample["t"]))

	#println("Starting backsolving.")
	#println("parameter values: ", parameter_values)
	#println("midpoint initial conditions ", initial_conditions)
	#println("saveat: ", saveat)
	ode_solution = ModelingToolkit.solve(prob, solver, p = parameter_values, saveat = saveat, abstol = abstol, reltol = reltol)

	state_param_map = (Dict(x => replace(string(x), "(t)" => "")
							for x in ModelingToolkit.states(model)))


	newstates = copy(E.states)

	for s in ModelingToolkit.states(model)
		#	println("BLAH")
		#	println(tspan)
		#	println("TEST ", state_param_map[s])
		#	println(" E[s] ", E[s])
		#	println(ode_solution[begin])
		#	println(ode_solution[Symbol(state_param_map[s])])

		#	println(ode_solution[Symbol(state_param_map[s])][end])
		temp = ode_solution[Symbol(state_param_map[s])][end]
		newstates[s] = temp
	end
	#println("startpoint initial conditions: ", newstates)
	ER = EstimationResult(E.parameters, newstates, E.degree, report_time,
		E.err, E.interpolants, E.return_code, E.datasize, report_time)
	#	println(ER.parameters)
	#	println("TEST")
	#	println(ER.states)
	#	println(ER)
	return ER

end






"""
	estimate_fixed_degree(model::ModelingToolkit.ODESystem,
						measured_quantities::Vector{ModelingToolkit.Equation},
						inputs::Vector{ModelingToolkit.Equation},
						data_sample::Dict{Any, Vector{T}} = Dict{Any, Vector{T}}();
						identifiability_result = Dict{String, Any}(),
						interpolation_degree::Int = 1,
						at_time::T = 0.0,
						method = :homotopy,
						real_tol = 1e-12) where {T <: Float}

Estimate the parameters of a model using the data sample `data_sample` and the
measured quantities `measured_quantities`.

# Arguments
- `model::ModelingToolkit.ODESystem`: the model with parameters and initial conditions to be estimated.
- `measured_quantities::Vector{ModelingToolkit.Equation}`: the measured quantities of the model. Used for identifiability assessment.
- `inputs::Vector{ModelingToolkit.Equation}`: the input equations of the model.
- `data_sample::Dict{Any, Vector{T}} = Dict{Any, Vector{T}}()`: the data sample used for estimation (same functions as `measured_quantities`).
																The keys of the dictionary are the measured quantities
																and the values are the corresponding data samples.
- `identifiability_result = Dict{String, Any}()`: the result of the identifiability assessment.
- `interpolation_degree::Int = 1`: the degree of the polynomial interpolation.
- `at_time::T = 0.0`: the time at which the parameters are estimated.
- `method = :homotopy`: the method used for solving the polynomial system. Can be one of :homotopy (recommended) and :msolve.
- `real_tol = 1e-12`: the tolerance for the real solutions of the polynomial system.

# Returns
- `EstimationResult`: the estimated parameters and initial conditions of the model.
"""
function estimate_single_interpolator(model::ModelingToolkit.ODESystem,
	measured_quantities::Vector{ModelingToolkit.Equation},
	inputs::Vector{ModelingToolkit.Equation},
	data_sample::AbstractDict{Any, Vector{T}} = Dict{Any,
		Vector{T}}();
	identifiability_result = Dict{String, Any}(),
	interpolator = ("AAA" => aaad),
	at_time::T = 0.0, report_time = minimum(data_sample["t"]),
	method = :homotopy,
	real_tol = 1e-12) where {T <: Float}
	time_interval = [minimum(data_sample["t"]), maximum(data_sample["t"])]  #TODO(orebas) will this break if key is missing?


	check_inputs(measured_quantities, data_sample)  #TODO(orebas): I took out checking the degree.  Do we want to check the interpolator otherwise?
	datasize = length(first(values(data_sample)))
	parameters = ModelingToolkit.parameters(model)
	states = ModelingToolkit.states(model)
	num_parameters = length(parameters) + length(states)
	@debug "Interpolating sample data via interpolation method $(interpolator.first)"
	if !haskey(data_sample, "t")
		@warn "No sampling time points found in data sample. Assuming uniform sampling t ∈ [$(time_interval[1]), $(time_interval[2])]."
		data_sample["t"] = range(time_interval[1], time_interval[2], length = datasize)
	end
	interpolants = ParameterEstimation.interpolate(identifiability_result,
		data_sample, measured_quantities, inputs;
		interpolator = interpolator,
		diff_order = num_parameters + 1,   #todo(OREBAS): is this always forcing num_parameters + 1 derivatives?
		at_t = at_time,
		method = method)

	if method == :homotopy
		all_solutions = solve_via_homotopy(identifiability_result, model;
			real_tol = real_tol)
	elseif method == :msolve
		all_solutions = solve_via_msolve(identifiability_result, model;
			real_tol = real_tol)
	else
		throw(ArgumentError("Method $method not supported"))
	end
	#println(all_solutions)
	#println("HERE")
	all_solutions = [EstimationResult(model, each, interpolator.first, at_time,
		interpolants, ReturnCode.Success, datasize, report_time)
					 for each in all_solutions]

	all_solutions_R = [backsolve_initial_conditions(model, each, report_time, inputs, data_sample)
					   for each in all_solutions]


	#println(all_solutions_R)
	return all_solutions_R
end

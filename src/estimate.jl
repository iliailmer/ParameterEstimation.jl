function estimate(model::ModelingToolkit.ODESystem,
                  measured_quantities::Vector{ModelingToolkit.Equation} = Vector{
                                                                                 ModelingToolkit.Equation
                                                                                 }([]),
                  data_sample = [],
                  time_interval = [],
                  interpolation_degree::Int = 1)
    if length(measured_quantities) == 0
        error("No measured states provided")
    end
    if length(data_sample) == 0
        error("No data sample provided")
    end
    if length(time_interval) > 2 || length(time_interval) == 1
        error("Time interval must be of the form [start, end]")
    end
    if length(time_interval) == 0
        error("No time interval provided")
    end
    parameters = ModelingToolkit.parameters(model)
    states = ModelingToolkit.states(model)
    num_parameters = length(parameters) + length(states)

    identifiability_result = ParameterEstimation.get_identifiability(model;
                                                                     measured_quantities = measured_quantities)
    # interpolate the data sample
    tsteps = range(time_interval[1], time_interval[2], length = length(data_sample))
    @info "Interpolating sample data"
    I = ParameterEstimation.interpolate(tsteps, data_sample, interpolation_degree)
    dIdt = ParameterEstimation.differentiate_interpolated(I, num_parameters + 1)
    err = sum(abs.(data_sample - I.(tsteps))) / length(tsteps)
    @info "Mean Absolute error in interpolation: $err"

    y_derivs_vals = Dict(ParameterEstimation.nemo2hc(key) => dIdt[val] * factorial(val)
                         for (key, val) in pairs(identifiability_result["Y_eq"]))
    polynomial_system = System(HomotopyContinuation.evaluate(ParameterEstimation.nemo2hc.(identifiability_result["Et"]),
                                                             y_derivs_vals))
    @info "HomotopyContinuations: solving $(length(polynomial_system)) polynomial equations in $(length(polynomial_system.variables)) variables"
    results = HomotopyContinuation.solve(polynomial_system)
    Base.show(results)

    all_solutions = HomotopyContinuation.real_solutions(results)
    if length(all_solutions) == 0
        all_solutions = HomotopyContinuation.solutions(results)
        if length(all_solutions) == 0
            @warn "No solutions found"
        end
    end
    all_solutions_dict = []
    state_str_map = Dict(string(x)[1:(end - 3)] => x for x in states)
    param_str_map = Dict(string(x) => x for x in parameters)
    state_param_map = merge(state_str_map, param_str_map)
    # remove unnecessary values, e.g. derivatives, etc.
    for sol in all_solutions
        tmp = Dict()
        for (idx, v) in enumerate(polynomial_system.variables)
            if endswith(string(v), "_0")
                tmp[state_param_map[string(v)[1:(end - 2)]]] = sol[idx]
            end
        end
        push!(all_solutions_dict, tmp)
    end
    return all_solutions_dict
end
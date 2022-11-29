function estimate(model::ModelingToolkit.ODESystem,
                  measured_quantities::Vector{ModelingToolkit.Equation} = Vector{
                                                                                 ModelingToolkit.Equation
                                                                                 }([]),
                  data_sample::Dict{Num, Vector{T}} = Dict{Num, Vector{T}}(),
                  time_interval = Vector{T}(),
                  identifiability_result = Dict{String, Any}(),
                  interpolation_degree::Int = 1) where {T <: Float}
    check_inputs(measured_quantities, data_sample, time_interval, interpolation_degree)

    parameters = ModelingToolkit.parameters(model)
    states = ModelingToolkit.states(model)
    num_parameters = length(parameters) + length(states)
    @info "Interpolating sample data"
    polynomial_system = ParameterEstimation.interpolate(identifiability_result,
                                                        data_sample,
                                                        time_interval,
                                                        measured_quantities,
                                                        interpolation_degree,
                                                        num_parameters + 1)

    @info "HomotopyContinuations: solving $(length(polynomial_system)) polynomial equations in $(length(polynomial_system.variables)) variables"
    results = HomotopyContinuation.solve(polynomial_system; show_progress = false)
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

function estimate_over_degrees(model::ModelingToolkit.ODESystem,
                               measured_quantities::Vector{ModelingToolkit.Equation} = Vector{
                                                                                              ModelingToolkit.Equation
                                                                                              }([]),
                               data_sample::Dict{Num, Vector{T}} = Dict{Num, Vector{T}}(),
                               time_interval = Vector{T}();
                               degree_range = nothing) where {T <: Float}
    check_inputs(measured_quantities, data_sample, time_interval)
    degree_range = 1:(length(data_sample[first(keys(data_sample))]) - 1)
    @info "Estimating over degrees between 1 and $(length(degree_range))"
    logger = ConsoleLogger(stdout, Logging.Error)
    estimates = nothing
    identifiability_result = ParameterEstimation.get_identifiability(model;
                                                                     measured_quantities = measured_quantities)
    with_logger(logger) do
        estimates = @showprogress pmap(x -> (x => filter_solutions(estimate(model,
                                                                            measured_quantities,
                                                                            data_sample,
                                                                            time_interval,
                                                                            identifiability_result,
                                                                            x), model,
                                                                   data_sample,
                                                                   time_interval)),
                                       degree_range)
    end
    best_solution = nothing
    for each in estimates
        if best_solution == nothing
            best_solution = each
        else
            if each[2][2] < best_solution[2][2]
                best_solution = each
            end
        end
    end
    @info "Best estimate $(best_solution[2][1]) found at degree $(best_solution[1]) with error $(best_solution[2][2])"
    return best_solution
end

function filter_solutions(results, model::ModelingToolkit.ODESystem,
                          data_sample::Dict{Num, Vector{T}} = Dict{Num, Vector{T}}(),
                          time_interval = Vector{T}()) where {T <: Float}
    @info "Filtering"
    min_error = 1e+10
    best_estimate = nothing
    tsteps = range(time_interval[1], time_interval[2],
                   length = length(first(values(data_sample))))
    @showprogress for (i, each_result) in enumerate(results)
        initial_conditions = [each_result[s] for s in ModelingToolkit.states(model)]
        parameter_values = [each_result[p] for p in ModelingToolkit.parameters(model)]
        prob = ModelingToolkit.ODEProblem(model, initial_conditions, time_interval,
                                          parameter_values)
        ode_solution = ModelingToolkit.solve(prob, Tsit5(),
                                             p = parameter_values,
                                             saveat = tsteps)
        if ode_solution.retcode == ReturnCode.Success
            err = 0
            for (key, sample) in data_sample
                err += ParameterEstimation.mean_abs_err(ode_solution[key], sample)
            end
            err /= length(data_sample)
            @info "\tSolution $i, error $err"
        else
            err = 1e+10
            @info "\tSolution $i, error ∞, probably unstable"
        end
        if err < min_error
            min_error = err
            best_estimate = each_result
        end
    end
    @info "Best estimate: $(best_estimate), error: $(min_error)"
    return best_estimate, min_error
end

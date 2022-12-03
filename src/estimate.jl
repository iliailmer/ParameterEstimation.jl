function estimate(model::ModelingToolkit.ODESystem,
                  measured_quantities::Vector{ModelingToolkit.Equation},
                  data_sample::Dict{Num, Vector{T}} = Dict{Num, Vector{T}}(),
                  time_interval = Vector{T}(),
                  identifiability_result = Dict{String, Any}(),
                  interpolation_degree::Int = 1; real_tol = 1e-10) where {T <: Float}
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
            @warn "Interpolation degree $(interpolation_degree): No solutions found"
            return []
        end
    end
    all_solutions_ = Vector{EstimationResult}([])
    state_str_map = Dict(replace(string(x), "(t)" => "") => x for x in states)
    param_str_map = Dict(string(x) => x for x in parameters)
    state_param_map = merge(state_str_map, param_str_map)
    # remove unnecessary values, e.g. derivatives, etc.
    for sol in all_solutions
        tmp = Dict()
        sol = map(each -> to_exact(each; tol = real_tol), sol)
        for (idx, v) in enumerate(polynomial_system.variables)
            if endswith(string(v), "_0")
                tmp[state_param_map[string(v)[1:(end - 2)]]] = sol[idx]
            end
        end
        param_est = EstimationResult(model, tmp, interpolation_degree, ReturnCode.Success)
        push!(all_solutions_, param_est)
    end

    return all_solutions_
end

function estimate_over_degrees(model::ModelingToolkit.ODESystem,
                               measured_quantities::Vector{ModelingToolkit.Equation},
                               data_sample::Dict{Num, Vector{T}} = Dict{Num, Vector{T}}(),
                               time_interval = Vector{T}();
                               degree_range = nothing, real_tol = 1e-10) where {T <: Float}
    check_inputs(measured_quantities, data_sample, time_interval)
    if degree_range == nothing
        degree_range = 1:(length(data_sample[first(keys(data_sample))]) - 1)
    end
    logger = ConsoleLogger(stdout, Logging.Warn)
    identifiability_result = ParameterEstimation.check_identifiability(model;
                                                                       measured_quantities = measured_quantities)
    estimates = []
    @info "Estimating over degrees between 1 and $(length(degree_range))"
    with_logger(logger) do
        @showprogress for deg in degree_range
            unfiltered = estimate(model,
                                  measured_quantities,
                                  data_sample,
                                  time_interval,
                                  identifiability_result,
                                  deg)
            if length(unfiltered) > 0
                filtered = filter_solutions(unfiltered, identifiability_result, model,
                                            data_sample, time_interval)
                push!(estimates, filtered)
            else
                push!(estimates,
                      EstimationResult(model, Dict(), deg, ReturnCode.Failure))
            end
        end
    end
    best_solution = nothing
    for est in estimates
        if !(est isa Vector)
            if est.return_code == ReturnCode.Success
                if best_solution == nothing
                    best_solution = est
                else
                    if est.err < best_solution.err
                        best_solution = est
                    end
                end
            end
        end
    end
    if best_solution != nothing
        @info "Best estimate found at degree $(best_solution.degree) with error $(best_solution.err)"
        return best_solution
    else
        @warn "No solution found"
        return estimates
    end
end

function solve_ode(model, estimate::EstimationResult, tsteps, data_sample; solver = Tsit5(),
                   return_ode = false)
    initial_conditions = [estimate[s] for s in ModelingToolkit.states(model)]
    parameter_values = [estimate[p] for p in ModelingToolkit.parameters(model)]
    prob = ModelingToolkit.ODEProblem(model, initial_conditions,
                                      (tsteps[1], tsteps[end]),
                                      parameter_values)
    ode_solution = ModelingToolkit.solve(prob, solver, saveat = tsteps)
    if ode_solution.retcode == ReturnCode.Success
        err = 0
        for (key, sample) in data_sample
            err += ParameterEstimation.mean_abs_err(ode_solution[key], sample)
        end
        err /= length(data_sample)
    else
        err = 1e+10
    end
    if return_ode
        return ode_solution,
               EstimationResult(estimate.parameters, estimate.states,
                                estimate.degree, err,
                                estimate.return_code)
    else
        return EstimationResult(estimate.parameters, estimate.states,
                                estimate.degree, err,
                                estimate.return_code)
    end
end

function solve_ode!(model, estimates::Vector{EstimationResult},
                    tsteps, data_sample; solver = Tsit5())
    estimates[:] = map(each -> solve_ode(model, each, tsteps, data_sample, solver = solver),
                       estimates)
end

function filter_solutions(results::Vector{EstimationResult},
                          identifiability_result::IdentifiabilityData,
                          model::ModelingToolkit.ODESystem,
                          data_sample::Dict{Num, Vector{T}} = Dict{Num, Vector{T}}(),
                          time_interval = Vector{T}(); topk = 1) where {T <: Float}
    @info "Filtering"
    if all(each -> each.return_code == ReturnCode.Failure, results)
        return results
    end
    min_error = 1e+10
    best_estimate = nothing
    tsteps = range(time_interval[1], time_interval[2],
                   length = length(first(values(data_sample))))
    if length(identifiability_result["identifiability"]["nonidentifiable"]) > 0
        @warn "The model contains non-identifiable parameters, no filtering was done."
        return results
    end
    solve_ode!(model, results, tsteps, data_sample)
    sorted = sort(results, by = x -> x.err)
    if topk == 1
        @info "Best estimate yelds ODE solution error $(sorted[1].err)"
        return sorted[1]
    else
        @info "Best $(topk) estimates yeld ODE solution errors $([s.err for s in sorted[1:topk]])"
        return sorted[1:topk]
    end
end

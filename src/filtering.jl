"""
    solve_ode(model, estimate::EstimationResult, sampling_times, data_sample; solver = Tsit5(),
              return_ode = false)

Solves the ODE system `model` with the parameters and initial conditions given by `estimate`.
Compute the error between the solution and the data sample. The error is recorded in the `EstimationResult`.

# Arguments
- `model`: the ODE system to be solved.
- `estimate::EstimationResult`: the parameters and initial conditions of the ODE system.
- `sampling_times`: the time steps of the ODE system. See `ModelingToolkit.solve`.
- `data_sample`: the data sample used for estimation (same functions as `measured_quantities`).
                 The keys of the dictionary are the measured quantities
                 and the values are the corresponding data samples.
- `solver = Tsit5()`: (optional) the solver used to solve the ODE system, see `DifferentialEquations` for available solvers.
- `return_ode = false`: (optional) whether to return the ODE solution.

# Returns
- ode_solution: the solution of the ODE system (if `return_ode` is set to `true`).
- `EstimationResult`: the estimated parameters and initial conditions of the model.
"""
function solve_ode(model, estimate::EstimationResult, data_sample;
                   solver = Tsit5(),
                   return_ode = false)
    initial_conditions = [estimate[s] for s in ModelingToolkit.states(model)]
    parameter_values = [estimate[p] for p in ModelingToolkit.parameters(model)]
    tspan = (estimate.at_time, data_sample["t"][end] + estimate.at_time)
    prob = ModelingToolkit.ODEProblem(model, initial_conditions,
                                      tspan, parameter_values)
    ode_solution = ModelingToolkit.solve(prob, solver,
                                         saveat = range(tspan[1], tspan[2],
                                                        length = length(data_sample["t"])))
    if ode_solution.retcode == ReturnCode.Success
        err = 0
        for (key, sample) in data_sample
            if isequal(key, "t")
                continue
            end
            err += ParameterEstimation.mean_abs_err(ode_solution(data_sample["t"])[key],
                                                    sample)
        end
        err /= length(data_sample)
    else
        err = 1e+10
    end
    if return_ode
        return ode_solution,
               EstimationResult(estimate.parameters, estimate.states,
                                estimate.degree, estimate.at_time, err,
                                estimate.interpolants, estimate.return_code,
                                estimate.datasize)
    else
        return EstimationResult(estimate.parameters, estimate.states,
                                estimate.degree, estimate.at_time, err,
                                estimate.interpolants, estimate.return_code,
                                estimate.datasize)
    end
end

"""
    solve_ode!(model, estimates::Vector{EstimationResult}, sampling_times, data_sample; solver = Tsit5())

Run solve_ode for multiple estimates and store the results (error between solution and sample) in each estimate.
This is done in-place.
"""
function solve_ode!(model, estimates::Vector{EstimationResult}, data_sample;
                    solver = Tsit5())
    estimates[:] = map(each -> solve_ode(model, each, data_sample, solver = solver),
                       estimates)
end

"""
    filter_solutions(results::Vector{EstimationResult},
                     identifiability_result::IdentifiabilityData,
                     model::ModelingToolkit.ODESystem,
                     data_sample::Dict{Num, Vector{T}} = Dict{Num, Vector{T}}(),
                     time_interval = Vector{T}(); topk = 1) where {T <: Float}

Filter estimation results stored in `results` vector based on ODE solving and checking against the sample.
In addition, takes into account global and local identifiability of parameters when filtering.

# Arguments
- `results::Vector{EstimationResult}`: the vector of estimation results.
- `identifiability_result::IdentifiabilityData`: the result of identifiability analysis.
- `model::ModelingToolkit.ODESystem`: the ODE system.
- `data_sample::Dict{Any, Any} = Dict{Any, Any}()`: the data sample used for estimation (same functions as `measured_quantities`).
                                                                The keys of the dictionary are the measured quantities
                                                                and the values are the corresponding data samples.
- `time_interval::Vector{T} = Vector{T}()`: the time interval of the ODE system.
- `topk = 1`: (optional) the number of best estimates to return.

# Returns
- `EstimationResult`: the best estimate (if `topk = 1`) or the vector of best estimates (if `topk > 1`).
"""
function filter_solutions(results::Vector{EstimationResult},
                          identifiability_result::IdentifiabilityData,
                          model::ModelingToolkit.ODESystem,
                          data_sample::Dict{Any, Any} = Dict{Any, Any}();
                          time_interval::Vector{T} = Vector{T}(),
                          solver = Tsit5(),
                          topk = 1) where {T <: Float}
    @info "Filtering"
    if length(results) == 0
        @warn "No results to filter."
        return results
    end
    if all(each -> each.return_code == ReturnCode.Failure, results)
        return results
    end
    try
        solve_ode!(model, results, data_sample; solver = solver) # this solves ODE with new parameters and computes err. between sample and solution
    catch InexactError
        @warn "InexactError when solving the ODE, no filtering was done."
        return results
    end
    if length(identifiability_result["identifiability"]["nonidentifiable"]) > 0
        @warn "The model contains non-identifiable parameters"
        filtered_results = Vector{ParameterEstimation.EstimationResult}()
        clustered = ParameterEstimation.cluster_estimates(model, results, data_sample,
                                                          solver = solver)
        if length(clustered) == 0
            @warn "No results to filter."
            return results
        end
        min_cluster_err, min_cluster_idx = findmin(sum(each.err for each in group) /
                                                   length(group)
                                                   for (id, group) in pairs(clustered))
        min_cluster = clustered[min_cluster_idx]
        @info "Best estimate yelds ODE solution error $(min_cluster_err)"
        filtered_results = min_cluster
        return results
    end
    if length(identifiability_result["identifiability"]["locally_not_globally"]) > 0
        filtered_results = Vector{ParameterEstimation.EstimationResult}()
        clustered = ParameterEstimation.cluster_estimates(model, results, data_sample,
                                                          solver = solver)
        if length(clustered) == 0
            @warn "No results to filter."
            return results
        end
        # find cluster with smallest error
        min_cluster_err, min_cluster_idx = findmin(sum(each.err for each in group) /
                                                   length(group)
                                                   for (id, group) in pairs(clustered))
        min_cluster = clustered[min_cluster_idx]
        @info "Best estimate yelds ODE solution error $(min_cluster_err)"
        filtered_results = min_cluster
    else
        sorted = sort(results, by = x -> x.err)
        if topk == 1
            filtered_results = Vector{EstimationResult}()
            @info "Best estimate yelds ODE solution error $(sorted[1].err)"
            push!(filtered_results, sorted[1])
        else
            filtered_results = Vector{Vector{ParameterEstimation.EstimationResult}}()
            @info "Best $(topk) estimates yeld ODE solution errors $([s.err for s in sorted[1:topk]])"
            push!(filtered_results, sorted[1:topk])
        end
    end
    return filtered_results
end

function cluster_estimates(model, res, data_sample; ε = 1e-6, solver = Tsit5())
    # clusrers the estimates by their error
    ParameterEstimation.solve_ode!(model, res, data_sample, solver = solver)
    clustered = Dict()
    #nearest neighbor search by err
    for i in eachindex(res)
        for j in (i + 1):length(res)
            if abs(res[i].err - res[j].err) < ε
                if !haskey(clustered, i)
                    clustered[i] = Vector{ParameterEstimation.EstimationResult}([])
                    push!(clustered[i], res[i])
                end
                push!(clustered[i], res[j])
            end
        end
    end
    # reset cluster keys
    clustered = Dict(i => clustered[key] for (i, key) in enumerate(keys(clustered)))
    return clustered
end


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

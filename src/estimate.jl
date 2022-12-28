"""
    estimate(model::ModelingToolkit.ODESystem,
             measured_quantities::Vector{ModelingToolkit.Equation},
             data_sample::Dict{Num, Vector{T}} = Dict{Num, Vector{T}}(),
             time_interval = Vector{T}(),
             identifiability_result = Dict{String, Any}(),
             interpolation_degree::Int = 1, at_time::T = 0.0;
             real_tol = 1e-10) where {T <: Float}

Estimate the parameters of a model using the data sample `data_sample` and the
measured quantities `measured_quantities`.

# Arguments
- `model::ModelingToolkit.ODESystem`: the model with parameters and initial conditions to be estimated.
- `measured_quantities::Vector{ModelingToolkit.Equation}`: the measured quantities of the model. Used for identifiability assessment.
- `data_sample::Dict{Num, Vector{T}} = Dict{Num, Vector{T}}()`: the data sample used for estimation (same functions as `measured_quantities`).
                                                                The keys of the dictionary are the measured quantities
                                                                and the values are the corresponding data samples.
- `time_interval = Vector{T}()`: the time interval of the data sample.
- `identifiability_result = Dict{String, Any}()`: the result of the identifiability assessment.
- `interpolation_degree::Int = 1`: the degree of the polynomial interpolation.
- `at_time::T = 0.0`: the time at which the parameters are estimated.
- `real_tol = 1e-10`: (optional) the tolerance for the real solutions of the polynomial system.

# Returns
- `EstimationResult`: the estimated parameters and initial conditions of the model.
"""
function estimate(model::ModelingToolkit.ODESystem,
                  measured_quantities::Vector{ModelingToolkit.Equation},
                  data_sample::Dict{Num, Vector{T}} = Dict{Num, Vector{T}}(),
                  time_interval = Vector{T}(),
                  identifiability_result = Dict{String, Any}(),
                  interpolation_degree::Int = 1, at_time::T = 0.0;
                  real_tol = 1e-10) where {T <: Float}
    check_inputs(measured_quantities, data_sample, time_interval, interpolation_degree)
    datasize = length(first(values(data_sample)))
    parameters = ModelingToolkit.parameters(model)
    states = ModelingToolkit.states(model)
    num_parameters = length(parameters) + length(states)
    @info "Interpolating sample data via rational interpolation"
    polynomial_system, interpolants = ParameterEstimation.interpolate(identifiability_result,
                                                                      data_sample,
                                                                      time_interval,
                                                                      measured_quantities,
                                                                      interpolation_degree,
                                                                      num_parameters + 1,
                                                                      at_time)

    @info "HomotopyContinuations: solving $(length(polynomial_system)) polynomial equations in $(length(polynomial_system.variables)) variables"
    results = HomotopyContinuation.solve(polynomial_system; show_progress = false)
    all_solutions = HomotopyContinuation.real_solutions(results)
    if length(all_solutions) == 0
        all_solutions = HomotopyContinuation.solutions(results)
        if length(all_solutions) == 0
            @warn "Interpolation numerator degree $(interpolation_degree): No solutions found"
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
        param_est = EstimationResult(model, tmp, interpolation_degree, at_time,
                                     interpolants, ReturnCode.Success, datasize)
        push!(all_solutions_, param_est)
    end
    return all_solutions_
end

# Path: src/estimate.jl
"""
    estimate_over_degrees(model::ModelingToolkit.ODESystem,
                          measured_quantities::Vector{ModelingToolkit.Equation},
                          data_sample::Dict{Num, Vector{T}} = Dict{Num, Vector{T}}(),
                          time_interval = Vector{T}(), at_time::T = 0.0;
                          solver = Tsit5(),
                          degree_range = nothing, real_tol = 1e-10) where {T <: Float}

Run estimation over a range of interpolation degrees. Return the best estimate according to a heuristic:
    - the best estimate is the one with the smallest error between sample data and ODE solution with current parameters (estimates);
"""
function estimate_over_degrees(model::ModelingToolkit.ODESystem,
                               measured_quantities::Vector{ModelingToolkit.Equation},
                               data_sample::Dict{Num, Vector{T}} = Dict{Num, Vector{T}}(),
                               time_interval = Vector{T}(), at_time::T = 0.0;
                               solver = Tsit5(),
                               degree_range = nothing,
                               real_tol = 1e-10) where {T <: Float}
    check_inputs(measured_quantities, data_sample, time_interval)
    datasize = length(first(values(data_sample)))
    if degree_range === nothing
        degree_range = 1:(length(data_sample[first(keys(data_sample))]) - 1)
    end
    logger = ConsoleLogger(stdout, Logging.Warn)
    identifiability_result = ParameterEstimation.check_identifiability(model;
                                                                       measured_quantities = measured_quantities)
    estimates = []
    @info "Estimating via rational interpolation with degrees between $(degree_range[1]) and $(degree_range[end])"
    with_logger(logger) do
        @showprogress for deg in degree_range
            unfiltered = estimate(model,
                                  measured_quantities,
                                  data_sample,
                                  time_interval,
                                  identifiability_result,
                                  deg, at_time)
            if length(unfiltered) > 0
                filtered = filter_solutions(unfiltered, identifiability_result, model,
                                            data_sample, time_interval; solver = solver)
                push!(estimates, filtered)
            else
                push!(estimates,
                      [
                          EstimationResult(model, Dict(), deg, at_time,
                                           Dict{Any, Interpolant}(),
                                           ReturnCode.Failure, datasize),
                      ])
            end
        end
    end
    #filter out the empty vectors
    estimates = filter(x -> length(x) > 0, estimates)

    #filter out the Failure results
    estimates = filter(x -> x[1].return_code == ReturnCode.Success, estimates)

    best_solution = nothing
    for each in estimates
        if best_solution === nothing
            best_solution = each
        else
            if sum(x.err for x in each) < sum(x.err for x in best_solution)
                best_solution = each
            end
        end
    end
    if best_solution !== nothing
        @info "Best estimate found" # error(s) $([x.err for x in best_solution])"
        return best_solution
    else
        @warn "No solution found"
        return estimates
    end
end

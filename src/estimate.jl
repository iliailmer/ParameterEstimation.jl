"""
estimate_fixed_degree(model::ModelingToolkit.ODESystem,
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
- `data_sample::Dict{Any, Vector{T}} = Dict{Any, Vector{T}}()`: the data sample used for estimation (same functions as `measured_quantities`).
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
function estimate_fixed_degree(model::ModelingToolkit.ODESystem,
                               measured_quantities::Vector{ModelingToolkit.Equation},
                               data_sample::Dict{Any, Vector{T}} = Dict{Any, Vector{T}}(),
                               identifiability_result = Dict{String, Any}(),
                               interpolation_degree::Int = 1, at_time::T = 0.0;
                               real_tol = 1e-10) where {T <: Float}
    time_interval = [minimum(data_sample["t"]), maximum(data_sample["t"])]
    check_inputs(measured_quantities, data_sample, interpolation_degree)
    datasize = length(first(values(data_sample)))
    parameters = ModelingToolkit.parameters(model)
    states = ModelingToolkit.states(model)
    num_parameters = length(parameters) + length(states)
    @info "Interpolating sample data via rational interpolation"
    if !haskey(data_sample, "t")
        @warn "No sampling time points found in data sample. Assuming uniform sampling t ∈ [$(time_interval[1]), $(time_interval[2])]."
        data_sample["t"] = range(time_interval[1], time_interval[2], length = datasize)
    end
    polynomial_system, interpolants = ParameterEstimation.interpolate(identifiability_result,
                                                                      data_sample,
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
            @debug "Interpolation numerator degree $(interpolation_degree): No solutions found"
            return Vector{EstimationResult}()
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
        for (key, val) in identifiability_result.transcendence_basis_subs
            if endswith(string(key), "_0")
                tmp[state_param_map[string(key)[1:(end - 2)]]] = Int(Meta.parse("$val"))
            end
        end
        param_est = EstimationResult(model, tmp, interpolation_degree, at_time,
                                     interpolants, ReturnCode.Success, datasize)
        push!(all_solutions_, param_est)
    end
    return all_solutions_
end

"""
estimate(model::ModelingToolkit.ODESystem,
                          measured_quantities::Vector{ModelingToolkit.Equation},
                          data_sample:::Dict{Any, Vector{T}} = Dict{Any, Vector{T}}(),
                          time_interval = Vector{T}(), at_time::T = 0.0;
                          solver = Tsit5(),
                          degree_range = nothing, real_tol = 1e-10) where {T <: Float}

Run estimation over a range of interpolation degrees. Return the best estimate according to a heuristic:
    - the best estimate is the one with the smallest error between sample data and ODE solution with current parameters (estimates);
"""
function estimate(model::ModelingToolkit.ODESystem,
                  measured_quantities::Vector{ModelingToolkit.Equation},
                  data_sample::Dict{Any, Vector{T}} = Dict{Any, Vector{T}}(),
                  at_time::T = 0.0;
                  solver = Tsit5(),
                  degree_range = nothing,
                  real_tol = 1e-10,
                  threaded = Threads.nthreads() > 1) where {T <: Float}
    if threaded
        result = estimate_threaded(model, measured_quantities, data_sample,
                                   at_time; solver = solver,
                                   degree_range = degree_range,
                                   real_tol = real_tol)
        display(result)
        return result
    else
        result = estimate_serial(model, measured_quantities, data_sample,
                                 at_time; solver = solver,
                                 degree_range = degree_range,
                                 real_tol = real_tol)
        display(result)
        return result
    end
end

function estimate_threaded(model, measured_quantities, data_sample,
                           at_time; solver = solver,
                           degree_range = degree_range,
                           real_tol = real_tol)
    @warn "Using threaded estimation."
    check_inputs(measured_quantities, data_sample)
    datasize = length(first(values(data_sample)))
    if degree_range === nothing
        degree_range = 1:(length(data_sample[first(keys(data_sample))]) - 1)
    end

    logger = ConsoleLogger(stdout, Logging.Warn)
    identifiability_result = ParameterEstimation.check_identifiability(model;
                                                                       measured_quantities = measured_quantities)
    n_threads = Threads.nthreads()
    N = length(degree_range)
    estimates = Vector{Any}(nothing, n_threads)
    @info "Estimating via rational interpolation with degrees between $(degree_range[1]) and $(degree_range[end]) using $n_threads threads"

    with_logger(logger) do
        Threads.@threads for t in 1:N
            deg = degree_range[t]
            id = Threads.threadid()
            if isnothing(estimates[id])
                estimates[id] = []
            end
            unfiltered = estimate_fixed_degree(model, measured_quantities,
                                               data_sample, identifiability_result,
                                               deg, at_time)
            if length(unfiltered) > 0
                filtered = filter_solutions(unfiltered, identifiability_result, model,
                                            data_sample; solver = solver)
                push!(estimates[id], filtered)

            else
                push!(estimates[id],
                      [
                          EstimationResult(model, Dict(), deg,
                                           at_time,
                                           Dict{Any, Interpolant}(),
                                           ReturnCode.Failure,
                                           datasize),
                      ])
            end
        end
    end
    #remove undef
    estimates = filter(x -> x !== nothing, estimates)
    estimates = vcat(estimates...)
    estimates = filter(x -> length(x) > 0, estimates)

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

function estimate_serial(model::ModelingToolkit.ODESystem,
                         measured_quantities::Vector{ModelingToolkit.Equation},
                         data_sample::Dict{Any, Vector{T}} = Dict{Any, Vector{T}}(),
                         at_time::T = 0.0; solver = Tsit5(), degree_range = nothing,
                         real_tol::Float64 = 1e-10) where {T <: Float}
    check_inputs(measured_quantities, data_sample)
    datasize = length(first(values(data_sample)))
    if degree_range === nothing
        degree_range = 1:(length(data_sample[first(keys(data_sample))]) - 1)
    end
    logger = ConsoleLogger(stdout, Logging.Warn)
    identifiability_result = ParameterEstimation.check_identifiability(model;
                                                                       measured_quantities = measured_quantities)
    estimates = Vector{Vector{ParameterEstimation.EstimationResult}}()
    @info "Estimating via rational interpolation with degrees between $(degree_range[1]) and $(degree_range[end])"
    with_logger(logger) do
        @showprogress for deg in degree_range
            unfiltered = estimate_fixed_degree(model, measured_quantities, data_sample,
                                               identifiability_result, deg, at_time)
            if length(unfiltered) > 0
                filtered = filter_solutions(unfiltered, identifiability_result,
                                            model, data_sample; solver = solver)
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
"""
    estimate_fixed_degree(model::ModelingToolkit.ODESystem,
                        measured_quantities::Vector{ModelingToolkit.Equation},
                        data_sample::Dict{Any, Vector{T}} = Dict{Any, Vector{T}}(),
                        identifiability_result = Dict{String, Any}(),
                        interpolation_degree::Int = 1, at_time::T = 0.0;
                        method = :homotopy_continuation,
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
                               method = :homotopy,
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
    interpolants = ParameterEstimation.interpolate(identifiability_result,
                                                   data_sample, measured_quantities,
                                                   interpolation_degree = interpolation_degree,
                                                   diff_order = num_parameters + 1,
                                                   at_t = at_time,
                                                   method = method)

    if method == :homotopy
        all_solutions = solve_via_homotopy(identifiability_result, model;
                                           real_tol = real_tol)
        all_solutions = [EstimationResult(model, each, interpolation_degree, at_time,
                                          interpolants, ReturnCode.Success, datasize)
                         for each in all_solutions]
    elseif method == :msolve
        all_solutions = solve_via_msolve(identifiability_result, model;
                                         real_tol = real_tol)

    else
        throw(ArgumentError("Method $method not supported"))
    end
    return all_solutions
end

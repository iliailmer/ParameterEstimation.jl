function estimate_serial(model::ModelingToolkit.ODESystem,
                         measured_quantities::Vector{ModelingToolkit.Equation},
                         data_sample::Dict{Any, Vector{T}} = Dict{Any, Vector{T}}();
                         at_time::T = 0.0, solver = Tsit5(), degree_range = nothing,
                         method = :homotopy,
                         real_tol::Float64 = 1e-10) where {T <: Float}
    check_inputs(measured_quantities, data_sample)
    datasize = length(first(values(data_sample)))
    if degree_range === nothing
        degree_range = 1:(length(data_sample[first(keys(data_sample))]) - 1)
    end
    id = ParameterEstimation.check_identifiability(model;
                                                   measured_quantities = measured_quantities)
    estimates = Vector{Vector{ParameterEstimation.EstimationResult}}()
    @info "Estimating via rational interpolation with degrees between $(degree_range[1]) and $(degree_range[end])"
    @showprogress for deg in degree_range
        unfiltered = estimate_fixed_degree(model, measured_quantities, data_sample;
                                           identifiability_result = id,
                                           interpolation_degree = deg, at_time = at_time,
                                           method = method, real_tol = real_tol)
        if length(unfiltered) > 0
            filtered = filter_solutions(unfiltered, id, model, data_sample; solver = solver)
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
    return post_process(estimates)
end

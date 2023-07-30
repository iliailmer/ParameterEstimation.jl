include("bary_derivs.jl")

function estimate_OB(model::ModelingToolkit.ODESystem,
    measured_quantities::Vector{ModelingToolkit.Equation},
    inputs::Vector{ModelingToolkit.Equation},
    data_sample::AbstractDict{Any,Vector{T}}=Dict{Any,Vector{T}}();
    at_time::T=0.0, solver=Tsit5(), degree_range=nothing,
    method=:homotopy,
    real_tol::Float64=1e-10) where {T<:Float}

    check_inputs(measured_quantities, data_sample)
    datasize = length(first(values(data_sample)))
    if degree_range === nothing
        degree_range = 1:(length(data_sample[first(keys(data_sample))])-1)
    end
    #skip identifiability check.  Add it back later!
    estimates = Vector{Vector{ParameterEstimation.EstimationResult}}()

    #start actual work here
    time_interval = [minimum(data_sample["t"]), maximum(data_sample["t"])]
    parameters = ModelingToolkit.parameters(model)
    states = ModelingToolkit.states(model)
    num_parameters = length(parameters) + length(states)

    if !haskey(data_sample, "t")
        @warn "No sampling time points found in data sample. Assuming uniform sampling t ∈ [$(time_interval[1]), $(time_interval[2])]."
        data_sample["t"] = range(time_interval[1], time_interval[2], length=datasize)
    end
    #polynomial_system = identifiability_result["polynomial_system"]
    interpolants = Dict{Any,Interpolant}()
    sampling_times = data_sample["t"]
    for (key, sample) in pairs(data_sample)  #copied from interpolate() but what is it doing?
        if key == "t"
            continue
        end
        y_function_name = map(x -> replace(string(x.lhs), "(t)" => ""),
            filter(x -> string(x.rhs) == string(key),
                measured_quantities))[1]
        aaa_interpolant = aaad(sampling_times, sample)
        interpolants[key] = aaa_interpolant
        err = sum(abs.(sample - aaa_interpolant.(sampling_times))) / length(sampling_times)
        @debug "Mean Absolute error in interpolation: $err interpolating $key"
        polynomial_system = eval_derivs(polynomial_system, interpolant, y_function_name,
            inputs, identifiability_result, method=method)
    end
    return post_process(estimates)


end



function estimate_serial(model::ModelingToolkit.ODESystem,
    measured_quantities::Vector{ModelingToolkit.Equation},
    inputs::Vector{ModelingToolkit.Equation},
    data_sample::AbstractDict{Any,Vector{T}}=Dict{Any,Vector{T}}();
    at_time::T=0.0, solver=Tsit5(), degree_range=nothing,
    method=:homotopy,
    real_tol::Float64=1e-10) where {T<:Float}

    check_inputs(measured_quantities, data_sample)
    #estimate_OB(model, measured_quantities, inputs, data_sample)


    datasize = length(first(values(data_sample)))
    if degree_range === nothing
        degree_range = 1:(length(data_sample[first(keys(data_sample))])-1)
    end
    id = ParameterEstimation.check_identifiability(model;
        measured_quantities=measured_quantities,
        inputs=[Num(each.lhs)
                for each in inputs])
    estimates = Vector{Vector{ParameterEstimation.EstimationResult}}()
    @info "Estimating via rational interpolation with degrees between $(degree_range[1]) and $(degree_range[end])"
    @showprogress for deg in degree_range
        unfiltered = estimate_fixed_degree(model, measured_quantities, inputs, data_sample;
            identifiability_result=id,
            interpolation_degree=deg, at_time=at_time,
            method=method, real_tol=real_tol)
        if length(unfiltered) > 0
            filtered = filter_solutions(unfiltered, id, model, inputs, data_sample;
                solver=solver)
            push!(estimates, filtered)
        else
            push!(estimates,
                [
                    EstimationResult(model, Dict(), deg, at_time,
                        Dict{Any,Interpolant}(),
                        ReturnCode.Failure, datasize),
                ])
        end
    end
    return post_process(estimates)
end


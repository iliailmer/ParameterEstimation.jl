function estimate_threaded(model, measured_quantities, inputs, data_sample;
                           at_time::Float = 0.0, solver = solver,
                           degree_range = degree_range,
                           method = :homotopy, real_tol::Float64 = 1e-10)
    @warn "Using threaded estimation."
    check_inputs(measured_quantities, data_sample)
    datasize = length(first(values(data_sample)))
    if isnothing(degree_range)
        degree_range = 1:(length(data_sample[first(keys(data_sample))]) - 1)
    end

    identifiability_result = ParameterEstimation.check_identifiability(model;
                                                                       measured_quantities = measured_quantities,
                                                                       inputs = [Num(each.lhs)
                                                                                 for each in inputs])
    n_threads = Threads.nthreads()
    N = length(degree_range)
    estimates = Vector{Any}(nothing, n_threads)
    @info "Estimating via rational interpolation with degrees between $(degree_range[1]) and $(degree_range[end]) using $n_threads threads"
    Threads.@threads for t in 1:N
        deg = degree_range[t]
        id = Threads.threadid()
        if isnothing(estimates[id])
            estimates[id] = []
        end
        unfiltered = estimate_fixed_degree(model, measured_quantities, inputs, data_sample;
                                           identifiability_result = identifiability_result,
                                           interpolation_degree = deg, at_time = at_time,
                                           method = method,
                                           real_tol = real_tol)
        if length(unfiltered) > 0
            filtered = filter_solutions(unfiltered, identifiability_result, model, inputs,
                                        data_sample; solver = solver)
            push!(estimates[id], filtered)

        else
            push!(estimates[id],
                  [
                      EstimationResult(model, Dict(), deg, at_time,
                                       Dict{Any, Interpolant}(), ReturnCode.Failure,
                                       datasize),
                  ])
        end
    end
    return post_process(estimates)
end

"""
    estimate(model::ModelingToolkit.ODESystem,
            measured_quantities::Vector{ModelingToolkit.Equation},
            data_sample::Dict{Any, Vector{T}} = Dict{Any, Vector{T}}(),
            at_time::T = 0.0; method = :homotopy, solver = Tsit5(),
            degree_range = nothing, real_tol = 1e-10,
            threaded = Threads.nthreads() > 1) where {T <: Float}

Run estimation over a range of interpolation degrees. Return the best estimate according to a heuristic:
    - the best estimate is the one with the smallest error between sample data and ODE solution with current parameters (estimates);
"""
function estimate(model::ModelingToolkit.ODESystem,
                  measured_quantities::Vector{ModelingToolkit.Equation},
                  data_sample::Dict{Any, Vector{T}} = Dict{Any, Vector{T}}();
                  at_time::T = 0.0, method = :homotopy, solver = Tsit5(),
                  degree_range = nothing, real_tol = 1e-10,
                  threaded = Threads.nthreads() > 1) where {T <: Float}
    if !(method in [:homotopy, :msolve])
        throw(ArgumentError("Method $method is not supported, must be one of :homotopy or :msolve."))
    end
    if threaded
        result = estimate_threaded(model, measured_quantities, data_sample;
                                   at_time = at_time, solver = solver,
                                   degree_range = degree_range,
                                   method = method,
                                   real_tol = real_tol)
    else
        result = estimate_serial(model, measured_quantities, data_sample;
                                 solver = solver, at_time = at_time,
                                 degree_range = degree_range, method = method,
                                 real_tol = real_tol)
    end
    for each in result
        display(each)
    end
    return result
end

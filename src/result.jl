"""
    EstimationResult

A container for the results of an estimation.
Contains the estimated parameters and initial conditions (state values at a given time), the degree of the rational interpolation used,
the error between the estimated ODE solution and the sample data, and the return code.

# Fields
- `parameters::OrderedDict`: The estimated parameters.
- `states::OrderedDict`: The estimated initial conditions.
- `degree::Int64`: The degree of the rational interpolation used.
- `at_time::Float64`: The time at which the initial conditions are estimated.
- `err::Union{Nothing, Float64}`: The error between the estimated ODE solution and the sample data.
- `interpolants::Union{Nothing, Dict{Any, Interpolant}}`: The rational interpolants used to estimate the parameters and initial conditions.
- `return_code::Any`: The return code of the estimation.
"""
struct EstimationResult
    parameters::OrderedDict
    states::OrderedDict
    degree::Int64
    at_time::Float64
    err::Union{Nothing, Float64}
    interpolants::Union{Nothing, Dict{Any, Interpolant}}
    return_code::Any
    datasize::Int64
    function EstimationResult(model::ModelingToolkit.ODESystem,
                              poly_sol::Union{Dict, OrderedDict}, degree::Int64,
                              at_time::Float64, interpolants::Dict{Any, Interpolant},
                              return_code, datasize)
        parameters = OrderedDict{Any, Any}()
        states = OrderedDict{Any, Any}()
        for p in ModelingToolkit.parameters(model)
            parameters[ModelingToolkit.Num(p)] = get(poly_sol, p, nothing)
        end
        for s in ModelingToolkit.states(model)
            states[ModelingToolkit.Num(s)] = get(poly_sol, s, nothing)
        end
        new(parameters, states, degree, at_time, nothing, interpolants, return_code,
            datasize)
    end
    function EstimationResult(parameters::OrderedDict, states::OrderedDict, degree::Int64,
                              at_time::Float64, err, interpolants, return_code, datasize)
        new(parameters, states, degree, at_time, err, interpolants, return_code, datasize)
    end
end

Base.get(sol::EstimationResult, k) = get(sol.parameters, k, get(sol.states, k, nothing))

function Base.getindex(sol::EstimationResult, k)
    if haskey(sol.parameters, k)
        return sol.parameters[k]
    elseif haskey(sol.states, k)
        return sol.states[k]
    else
        throw(KeyError(k))
    end
end

function Base.show(io::IO, e::EstimationResult)
    if any(isnothing.(values(e.parameters)))
        println(io, "Parameter(s)        :\t",
                join([@sprintf("%3s = %3s", k, v) for (k, v) in pairs(e.parameters)],
                     ", "))
        println(io, "Initial Condition(s):\t",
                join([@sprintf("%3s = %3s", k, v) for (k, v) in pairs(e.states)], ", "))
    elseif !all(isreal.(values(e.parameters)))
        println(io, "Parameter(s)        :\t",
                join([@sprintf("%3s = %.3f+%.3fim", k, real(v), imag(v))
                      for (k, v) in pairs(e.parameters)],
                     ", "))
        println(io, "Initial Condition(s):\t",
                join([@sprintf("%3s = %.3f+%.3fim", k, real(v), imag(v))
                      for (k, v) in pairs(e.states)], ", "))
    else
        println(io, "Parameter(s)        :\t",
                join([@sprintf("%3s = %.3f", k, v) for (k, v) in pairs(e.parameters)],
                     ", "))
        println(io, "Initial Condition(s):\t",
                join([@sprintf("%3s = %.3f", k, v) for (k, v) in pairs(e.states)], ", "))
    end
    # println(io, "Interpolation Degree (numerator): ", e.degree)
    # println(io, "Interpolation Degree (denominator): ", e.datasize - e.degree - 1)
    # if isnothing(e.err)
    # println(io, "Error: Not yet calculated")
    # else
    # println(io, "Error: ", @sprintf("%.4e", e.err))
    # end
    # println(io, "Return Code: ", e.return_code)
end

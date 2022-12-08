"""
    EstimationResult

A container for the results of an estimation.
Contains the estimated parameters and initial conditions (state values at a given time), the degree of the rational interpolation used,
the error between the estimated ODE solution and the sample data, and the return code.
"""
struct EstimationResult
    parameters::OrderedDict
    states::OrderedDict
    degree::Int
    err::Union{Nothing, Float64}
    return_code::Any
    function EstimationResult(model, poly_sol, degree, return_code)
        parameters = OrderedDict{Any, Any}()
        states = OrderedDict{Any, Any}()
        for p in ModelingToolkit.parameters(model)
            parameters[ModelingToolkit.Num(p)] = get(poly_sol, p, nothing)
        end
        for s in ModelingToolkit.states(model)
            states[ModelingToolkit.Num(s)] = get(poly_sol, s, nothing)
        end
        new(parameters, states, degree, nothing, return_code)
    end
    function EstimationResult(parameters, states, degree, err, return_code)
        new(parameters, states, degree, err, return_code)
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
        println(io, "Parameter Estimates:\n\t",
                join([@sprintf("%s = %s", k, v) for (k, v) in pairs(e.parameters)],
                     ", "))
        println(io, "Initial Condition Estimates:\n\t",
                join([@sprintf("%s = %s", k, v) for (k, v) in pairs(e.states)], ", "))
    elseif !all(isreal.(values(e.parameters)))
        println(io, "Parameter Estimates:\n\t",
                join([@sprintf("%s = %.6f+%.6fim", k, real(v), imag(v))
                      for (k, v) in pairs(e.parameters)],
                     ", "))
        println(io, "Initial Condition Estimates:\n\t",
                join([@sprintf("%s = %.6f+%.6fim", k, real(v), imag(v))
                      for (k, v) in pairs(e.states)], ", "))
    else
        println(io, "Parameter Estimates:\n\t",
                join([@sprintf("%s = %.6f", k, v) for (k, v) in pairs(e.parameters)],
                     ", "))
        println(io, "Initial Condition Estimates:\n\t",
                join([@sprintf("%s = %.6f", k, v) for (k, v) in pairs(e.states)], ", "))
    end
    println(io, "Degree: ", e.degree)
    if isnothing(e.err)
        println(io, "Error: Not yet calculated")
    else
        println(io, "Error: ", @sprintf("%.4e", e.err))
    end
    println(io, "Return Code: ", e.return_code)
end

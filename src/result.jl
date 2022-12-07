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
    println(io, "Parameter Estimates:\n\t",
            join([@sprintf("%s = %.6f", k, v) for (k, v) in pairs(e.parameters)],
                 ", "))
    println(io, "Initial Condition Estimates:\n\t",
            join([@sprintf("%s = %.6f", k, v) for (k, v) in pairs(e.states)], ", "))
    println(io, "Degree: ", e.degree)
    if isequal(e.err, nothing)
        println(io, "Error: Not yet calculated")
    else
        println(io, "Error: ", @sprintf("%.4e", e.err))
    end
    println(io, "Return Code: ", e.return_code)
end

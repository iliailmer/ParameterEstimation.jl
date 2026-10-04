"""
    EstimationResult

A container for one estimate. Indexing with a parameter or a state returns its estimated
value, for example `result[mu]`.

# Fields
- `parameters::AbstractDict`: the estimated parameters.
- `states::AbstractDict`: the estimated states (initial conditions) at time `at_time`.
- `degree`: the name of the interpolator used, for example `"AAA"`.
- `at_time::Float64`: the time at which `states` are given. In the results of `estimate` it is equal to `report_time`.
- `err::Union{Nothing, Float64}`: the mean absolute error between the ODE solution with the estimated values and the data.
- `interpolants::Union{Nothing, AbstractDict{Any, Interpolant}}`: the interpolants of the measured quantities.
- `return_code`: `ReturnCode.Success` or `ReturnCode.Failure`.
- `datasize::Int64`: the number of data points.
- `report_time`: the time at which the states are reported (by default the first time point).
"""
struct EstimationResult
    parameters::AbstractDict
    states::AbstractDict
    degree::Any
    at_time::Float64
    err::Union{Nothing, Float64}
    interpolants::Union{Nothing, AbstractDict{Any, Interpolant}}
    return_code::Any
    datasize::Int64
    report_time::Any
    function EstimationResult(model::ModelingToolkit.System,
            poly_sol::AbstractDict, degree,
            at_time::Float64,
            interpolants::AbstractDict{Any, Interpolant},
            return_code, datasize, report_time)
        parameters = OrderedDict{Any, Any}()
        states = OrderedDict{Any, Any}()
        for p in ModelingToolkit.parameters(model)
            parameters[ModelingToolkit.Num(p)] = get(poly_sol, p, nothing)
        end
        for s in ModelingToolkit.unknowns(model)
            states[ModelingToolkit.Num(s)] = get(poly_sol, s, nothing)
        end
        new(parameters, states, degree, at_time, nothing, interpolants, return_code,
            datasize, report_time)
    end
    function EstimationResult(parameters::AbstractDict, states::AbstractDict, degree,
            at_time::Float64, err, interpolants, return_code, datasize, report_time::Float64)
        new(parameters, states, degree, at_time, err,
            interpolants, return_code, datasize, report_time)
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
    if (!isnothing(e.report_time))
        report_time_string = @sprintf(", where t = %.3f", e.report_time)
    else
        report_time_string = ""
    end

    if any(isnothing.(values(e.parameters)))
        println(io, "Parameter(s)        :\t",
            join([@sprintf("%3s = %3s", k, v) for (k, v) in pairs(e.parameters)],
                ", "))
        println(io, "Initial Condition(s):\t",
            join([@sprintf("%3s = %3s", k, v) for (k, v) in pairs(e.states)], ", "), report_time_string)
    elseif !all(isreal.(values(e.parameters)))
        println(io, "Parameter(s)        :\t",
            join(
                [@sprintf("%3s = %.3f+%.3fim", k, real(v), imag(v))
                 for (k, v) in pairs(e.parameters)],
                ", "))
        println(io, "Initial Condition(s):\t",
            join(
                [@sprintf("%3s = %.3f+%.3fim", k, real(v), imag(v))
                 for (k, v) in pairs(e.states)], ", "), report_time_string)
    else
        println(io, "Parameter(s)        :\t",
            join([@sprintf("%3s = %.3f", k, v) for (k, v) in pairs(e.parameters)],
                ", "))
        println(io, "Initial Condition(s):\t",
            join([@sprintf("%3s = %.3f", k, v) for (k, v) in pairs(e.states)], ", "), report_time_string)
    end
    # println(io, "Interpolation Degree (numerator): ", e.degree)
    # println(io, "Interpolation Degree (denominator): ", e.datasize - e.degree - 1)
    if isnothing(e.err)
        println(io, "Error: Not yet calculated")
    else
        println(io, "Error: ", @sprintf("%.4e", e.err))
    end
    #	if isnothing(e.at_time)
    #		println(io, "Time: Not specified")
    #	else
    #		println(io, "Time: ", @sprintf("%.4e", e.at_time))
    #	end
    # end
    # println(io, "Return Code: ", e.return_code)
end

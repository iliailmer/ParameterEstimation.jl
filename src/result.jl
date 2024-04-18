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
- `datasize::Int64`: The number of data points used in the estimation.
- `report_time::Any`: The time at which the initial conditions are reported (usually the leftmost point in the time span).

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
	function EstimationResult(model::ModelingToolkit.ODESystem,
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
		new(parameters, states, degree, at_time, err, interpolants, return_code, datasize, report_time)
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
			join([@sprintf("%3s = %.3f+%.3fim", k, real(v), imag(v))
				  for (k, v) in pairs(e.parameters)],
				", "))
		println(io, "Initial Condition(s):\t",
			join([@sprintf("%3s = %.3f+%.3fim", k, real(v), imag(v))
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

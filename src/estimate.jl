function estimate(model::ModelingToolkit.ODESystem, data_sample=[], time_interval=[], measured_states::Vector{Equation}=[], interpolation_degree::Int=1)
    if length(data_sample) == 0
        error("No data sample provided")
    end
    if length(measured_states) == 0
        error("No measured states provided")
    end
    if length(time_interval) > 2 || length(time_interval) == 1
        error("Time interval must be of the form [start, end]")
    end
    if length(time_interval) == 0
        error("No time interval provided")
    end

    identifiability_result = get_identifiability(model; measured_quantities=measured_states)
    parameters = ModelingToolkit.parameters(model)
    states = ModelingToolkit.states(model)
    num_parameters = length(parameters) + length(states)

    # interpolate the data sample
    tsteps = time_interval[1]
    data_interpolator = interpolate(tsteps, data_sample, interpolation_degree)
end


function interpolate(identifiability_result::Dict, num_parameters::Int)

end
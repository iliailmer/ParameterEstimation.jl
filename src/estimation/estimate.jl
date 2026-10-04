"""
    estimate(model::ModelingToolkit.System,
             measured_quantities::Vector{ModelingToolkit.Equation},
             data_sample::AbstractDict{Any, Vector{T}};
             inputs = ModelingToolkit.Equation[],
             at_time = data_sample["t"][length(data_sample["t"]) ÷ 2],
             method = :homotopy, solver = Vern9(),
             report_time = minimum(data_sample["t"]),
             interpolators = nothing, real_tol = 1e-14, seed = 42,
             threaded = Threads.nthreads() > 1, filtermode = :new,
             parameter_constraints = nothing,
             ic_constraints = nothing) where {T <: Float}

Estimate the parameters and initial conditions of `model` from `data_sample`.

The data of each measured quantity is interpolated. The derivatives of the interpolants at
`at_time` are substituted into a polynomial system obtained from identifiability analysis,
and the system is solved. Each solution is then checked: the ODE is solved with the
estimated values and compared with the data. This is repeated for every interpolator.

# Arguments
- `model::ModelingToolkit.System`: the ODE model;
- `measured_quantities::Vector{ModelingToolkit.Equation}`: the measured quantities (outputs), as equations such as `y ~ x^2 + x`;
- `data_sample`: the data, a dictionary with the time points under the key `"t"` and the samples of each measured quantity under the right-hand side of its equation (`x^2 + x` above).

# Keyword arguments
- `inputs`: known input functions, as equations such as `u ~ sin(t)`;
- `at_time`: the time point at which the derivatives are computed. Default: a point near the middle of the data;
- `report_time`: the time at which the estimated states (initial conditions) are reported. Default: the first time point;
- `method = :homotopy`: the polynomial system solver. Only `:homotopy` is implemented;
- `solver = Vern9()`: the ODE solver used to check the estimates;
- `interpolators = nothing`: a dictionary `name => function` of interpolators. If `nothing`, AAA and Floater-Hormann interpolation are used;
- `real_tol = 1e-14`: real and imaginary parts of a solution smaller than this are set to zero;
- `seed = 42`: the random seed for identifiability analysis and polynomial system solving. Results are reproducible for a fixed seed; the caller's random number state is not changed. If expected solutions are missing, try another seed;
- `threaded = Threads.nthreads() > 1`: run the interpolators in parallel threads;
- `parameter_constraints = nothing`: a dictionary `parameter => (lower, upper)`. Estimates outside the bounds are dropped;
- `ic_constraints = nothing`: the same for the states;
- `filtermode = :new`: how the estimates are selected. `:new` returns all estimates that fit the data.

# Returns
- `Vector{EstimationResult}`: the estimates, sorted by `err` (best first). The vector is empty if no estimate fits the data. Models that are only locally identifiable can return several estimates that fit the data equally well.
"""
function estimate(model::ModelingToolkit.System,
        measured_quantities::Vector{ModelingToolkit.Equation},
        data_sample::AbstractDict{Any, Vector{T}} = Dict{Any, Vector{T}}();
        inputs::Vector{ModelingToolkit.Equation} = Vector{ModelingToolkit.Equation}(),
        at_time::T = data_sample["t"][fld(length((data_sample["t"])), 2)],  #uses something akin to a midpoint by default
        method = :homotopy, solver = Vern9(),
        report_time = minimum(data_sample["t"]),
        interpolators = nothing, real_tol = 1e-14, seed = 42,
        threaded = Threads.nthreads() > 1, filtermode = :new, parameter_constraints = nothing,
        ic_constraints = nothing) where {T <: Float}
    if !(method in [:homotopy, :msolve])
        throw(ArgumentError("Method $method is not supported, must be one of :homotopy or :msolve."))
    end
    rng_state = copy(Random.default_rng())
    Random.seed!(seed)
    result = try
        if threaded
            estimate_threaded(model, measured_quantities, inputs, data_sample;
                at_time = at_time, solver = solver, report_time,
                interpolators = interpolators, method = method,
                real_tol = real_tol, filtermode, parameter_constraints = parameter_constraints,
                ic_constraints = ic_constraints)
        else
            estimate_serial(model, measured_quantities,
                inputs,
                data_sample;
                solver = solver, at_time = at_time, report_time,
                interpolators = interpolators, method = method,
                real_tol = real_tol, filtermode, parameter_constraints = parameter_constraints,
                ic_constraints = ic_constraints)
        end
    finally
        copy!(Random.default_rng(), rng_state)
    end
    println("Final Results:")
    for each in result
        display(each)
    end
    return result
end

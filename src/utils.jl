"""
    nemo2hc(expr_tree::Union{Expr, Symbol})

Converts a symbolic expression from Nemo to HomotopyContinuation format.
"""
function nemo2hc(expr_tree::Union{Expr, Symbol})
    #traverse expr_tree
    if typeof(expr_tree) == Symbol
        return HomotopyContinuation.Expression(HomotopyContinuation.variables(expr_tree)[1])
    end
    if typeof(expr_tree) == Expr
        if expr_tree.head == :call
            if expr_tree.args[1] in [:+, :-, :*, :/, :^, ://]
                if length(expr_tree.args) == 2
                    return eval(expr_tree.args[1])(nemo2hc(expr_tree.args[2]))
                else
                    return reduce(eval(expr_tree.args[1]),
                                  map(nemo2hc, expr_tree.args[2:end]))
                end
            end
        end
    end
end

function nemo2hc(expr_tree::SIAN.Nemo.fmpq_mpoly)
    return nemo2hc(Meta.parse(string(expr_tree)))
end

function nemo2hc(expr_tree::Number)
    return expr_tree
end

"""
    squarify_system(poly_system::Vector{Expression})

Given a non-square polynomial system in `n` variables, takes first `n - 1` equations and
adds a random linear combination of the remaining equations to the system.
"""
function squarify_system(poly_system::Vector{Expression})
    indets = HomotopyContinuation.variables(poly_system)
    M = randn(1, length(poly_system) - length(indets) + 1)
    return vcat(poly_system[1:(length(indets) - 1)], M * poly_system[length(indets):end])
end

"""
    check_inputs(measured_quantities::Vector{ModelingToolkit.Equation} = Vector{ModelingToolkit.Equation}([]),
                 data_sample::Dict{Num, Vector{T}} = Dict{Num, Vector{T}}(),
                 time_interval = Vector{T}(),
                 interpolation_degree::Int = 1) where {T <: Float}

Checks that the inputs to `estimate` are valid.
"""
function check_inputs(measured_quantities::Vector{ModelingToolkit.Equation} = Vector{
                                                                                     ModelingToolkit.Equation
                                                                                     }([]),
                      data_sample::Dict{Num, Vector{T}} = Dict{Num, Vector{T}}(),
                      time_interval = Vector{T}(),
                      interpolation_degree::Int = 1) where {T <: Float}
    if length(measured_quantities) == 0
        error("No measured states provided")
    end
    if length(data_sample) == 0
        error("No data sample provided")
    end
    if length(time_interval) > 2 || length(time_interval) == 1
        error("Time interval must be of the form [start, end]")
    end
    if length(time_interval) == 0
        error("No time interval provided")
    end
    if interpolation_degree < 1
        error("Interpolation degree must be ≥ 1")
    end
end

to_real(x::Number; tol = 1e-10) = abs(imag(x)) < tol ? real(x) : x
function to_exact(x::Number; tol = 1e-10)
    r = abs(real(x)) < tol ? 0 : real(x)
    i = abs(imag(x)) < tol ? 0 : imag(x)
    if i == 0
        return r
    else
        return r + i * im
    end
end

function sample_data(model::ModelingToolkit.ODESystem,
                     measured_data::Vector{ModelingToolkit.Equation},
                     time_interval::Vector{T},
                     p_true::Vector{T},
                     u0::Vector{T},
                     num_points::Int;
                     solver = Tsit5(), inject_noise = false, mean_noise = 0,
                     stddev_noise = 1) where {T <: Float}
    tsteps = range(time_interval[1], time_interval[2], length = num_points)
    problem = ODEProblem(model, u0, time_interval, p_true)
    solution_true = ModelingToolkit.solve(problem, solver, p = p_true, saveat = tsteps;
                                          abstol = 1e-10, reltol = 1e-10)
    data_sample = Dict(Num(v.rhs) => solution_true[Num(v.rhs)] for v in measured_data)
    if inject_noise
        for (key, sample) in data_sample
            data_sample[key] = sample + randn(num_points) .* stddev_noise .+ mean_noise
        end
    end

    return data_sample
end

function write_sample(data_sample; filename = "sample_data.txt")
    open(filename, "w") do io
        idx_iter = eachindex(first(values(data_sample)))
        # print keys
        for (key, sample) in data_sample
            print(io, key, " ")
        end
        println(io)
        for i in idx_iter
            for (key, sample) in data_sample
                print(io, sample[i], " ")
            end
            println(io)
        end
    end
end
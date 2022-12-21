using ModelingToolkit, Flux, DiffEqFlux, DifferentialEquations, Plots
using Distributions, Random
solver = Tsit5()

@parameters k1 k2 k3 k4 k5 k6
@variables t x1(t) x2(t) x3(t) x4(t) x5(t) x6(t) y1(t) y2(t)
D = Differential(t)
states = [x1, x2, x3, x4, x5, x6]
parameters = [k1, k2, k3, k4, k5, k6]

@named model = ODESystem([
                             D(x1) ~ -k1 * x1 * x2 + k2 * x4 + k4 * x6,
                             D(x2) ~ -k1 * x1 * x2 + k2 * x4 + k3 * x4,
                             D(x3) ~ k3 * x4 + k5 * x6 - k6 * x3 * x5,
                             D(x4) ~ k1 * x1 * x2 - k2 * x4 - k3 * x4,
                             D(x5) ~ k4 * x6 + k5 * x6 - k6 * x3 * x5,
                             D(x6) ~ -k4 * x6 - k5 * x6 + k6 * x3 * x5,
                         ], t, states, parameters)
measured_quantities = [y1 ~ x3, y2 ~ x2]

ic = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
p_true = [0.03, 0.02, 0.05, 0.03, 0.02, 0.05] # True Parameters

time_interval = [0.0, 30.0]
datasize = 10
tsteps = range(time_interval[1], time_interval[2], length = datasize)

prob_true = ODEProblem(model, ic, time_interval, p_true)
solution_true = solve(prob_true, solver, p = p_true, saveat = tsteps)

# Initial Parameter Vector
p = randn!(MersenneTwister(123), zeros(length(p_true) + length(ic)))
_params = Flux.params(p)

prob = ODEProblem(model, [p[1], p[2], p[3], p[4], p[5], p[6]], time_interval, p[7:end])

function predict_rd() # Our 1-layer "neural network"
    sol = solve(remake(prob, u0 = [p[1], p[2], p[3], p[4], p[5], p[6]]), solver,
                p = p[7:end], saveat = tsteps) # override with new parameters
    return [sol[3, :] sol[2, :]]
end

function loss_rd()
    y_pred = predict_rd()
    y_true = [solution_true[3, :] solution_true[2, :]]
    return sum(abs2, y_true .- y_pred)
end # loss function
data = Iterators.repeated((), 3000)
opt = ADAM(0.01)
cb = function () #callback function to observe training
    display(loss_rd())
    # using `remake` to re-create our `prob` with current parameters `p`
    # display(plot(solve(remake(prob, p=p), AutoTsit5(Rosenbrock23()), saveat=tsteps), ylim=(0, 6)))
end

# Display the ODE with the initial parameter values.
cb()

Flux.train!(loss_rd, _params, data, opt, cb = cb)

println("Parameters: ", p)
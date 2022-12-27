using ModelingToolkit, Flux, DiffEqFlux, DifferentialEquations
using Distributions, Random
solver = Tsit5()

@parameters a12 a13 a21 a31 a01
@variables t x1(t) x2(t) x3(t) y1(t) y2(t)
D = Differential(t)

ic = [1.0, 2.0, 1.0]
time_interval = [0.0, 1.0]
datasize = 10
tsteps = range(time_interval[1], time_interval[2], length = datasize)
p_true = [0.2, 0.3, 0.5, 0.6, -0.2] # True Parameters

states = [x1, x2, x3]
parameters = [a12, a13, a21, a31, a01]
@named model = ODESystem([D(x1) ~ -(a21 + a31 + a01) * x1 + a12 * x2 + a13 * x3,
                             D(x2) ~ a21 * x1 - a12 * x2,
                             D(x3) ~ a31 * x1 - a13 * x3],
                         t, states, parameters)
measured_quantities = [y1 ~ x1, y2 ~ x2]
prob_true = ODEProblem(model, ic, time_interval, p_true)
solution_true = solve(prob_true, solver, p = p_true, saveat = tsteps)

# Initial Parameter Vector
p = randn!(MersenneTwister(123), zeros(length(p_true) + length(ic)))
_params = Flux.params(p)

prob = ODEProblem(model, [p[1], p[2], p[3]], time_interval, p[4:end])

function predict_rd() # Our 1-layer "neural network"
    sol = solve(remake(prob, u0 = [p[1], p[2], p[3]]), solver, p = p[4:end],
                saveat = tsteps) # override with new parameters
    return [sol[1, :] sol[2, :]]
end

function loss_rd()
    y_pred = predict_rd()
    y_true = [solution_true[1, :] solution_true[2, :]]
    return sum(abs, y_true .- y_pred)
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
println("parameters:", p)
using ModelingToolkit, Flux, DiffEqFlux, DifferentialEquations, Plots
using Distributions, Random
solver = Tsit5()

@parameters k01 k12 k13 k14 k21 k31 k41
@variables t x1(t) x2(t) x3(t) x4(t) y1(t) y2(t) y3(t) y4(t)
D = Differential(t)

ic = [1.0, 2.0, 1.0, -1.0]
time_interval = [0.0, 10.0]
datasize = 20
tsteps = range(time_interval[1], time_interval[2], length = datasize)
p_true = [0.2, 0.3, 0.5, 0.6, -0.2, 1.1, 0.02] # True Parameters

states = [x1, x2, x3, x4]
parameters = [k01, k12, k13, k14, k21, k31, k41]
@named model = ODESystem([
                             D(x1) ~ -k01 * x1 + k12 * x2 + k13 * x3 + k14 * x4 - k21 * x1 -
                                     k31 * x1 - k41 * x1,
                             D(x2) ~ -k12 * x2 + k21 * x1,
                             D(x3) ~ -k13 * x3 + k31 * x1,
                             D(x4) ~ -k14 * x4 + k41 * x1],
                         t, states, parameters)
measured_quantities = [y1 ~ x1 + x3, y2 ~ x2 + x4, y3 ~ x1 + x2, y4 ~ x3 + x4]
prob_true = ODEProblem(model, ic, time_interval, p_true)
solution_true = solve(prob_true, solver, p = p_true, saveat = tsteps)

# Initial Parameter Vector
p = randn!(MersenneTwister(123), zeros(length(p_true) + length(ic)))
_params = Flux.params(p)

prob = ODEProblem(model, [p[1], p[2], p[3], p[4]], time_interval, p[5:end])

function predict_rd() # Our 1-layer "neural network"
    sol = solve(remake(prob, u0 = [p[1], p[2], p[3], p[4]]), solver, p = p[5:end],
                saveat = tsteps) # override with new parameters
    return [sol[1, :] + sol[3, :] sol[2, :] + sol[4, :] sol[1, :] + sol[2, :] sol[3, :] +
                                                                              sol[4, :]]
end

function loss_rd()
    y_pred = predict_rd()
    y_true = [solution_true[1, :] + solution_true[2, :] solution_true[3, :] +
                                                        solution_true[4, :] solution_true[1,
                                                                                          :] +
                                                                            solution_true[3,
                                                                                          :] solution_true[2,
                                                                                                           :] +
                                                                                             solution_true[4,
                                                                                                           :]]
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
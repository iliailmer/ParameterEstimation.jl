using ModelingToolkit, Flux, DiffEqFlux, DifferentialEquations, Plots
using Distributions, Random
solver = Tsit5()

@parameters lm d beta a k u c q b h
@variables t x(t) y(t) v(t) w(t) z(t) y1(t) y2(t) y3(t) y4(t)
D = Differential(t)
states = [x, y, v, w, z]
parameters = [lm, d, beta, a, k, u, c, q, b, h]

@named model = ODESystem([
                             D(x) ~ lm - d * x - beta * x * v,
                             D(y) ~ beta * x * v - a * y,
                             D(v) ~ k * y - u * v,
                             D(w) ~ c * x * y * w - c * q * y * w - b * w,
                             D(z) ~ c * q * y * w - h * z,
                         ], t, states, parameters)
measured_quantities = [y1 ~ w, y2 ~ z, y3 ~ x, y4 ~ y + v]

ic = [1.0, 1.0, 1.0, 1.0, 1.0]
time_interval = [0.0, 10.0]
datasize = 100
tsteps = range(time_interval[1], time_interval[2], length = datasize)
p_true = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]

prob_true = ODEProblem(model, ic, time_interval, p_true)
solution_true = solve(prob_true, solver, p = p_true, saveat = tsteps; abstol = 1e-12,
                      reltol = 1e-12)

# Initial Parameter Vector
p = randn!(MersenneTwister(123), zeros(length(p_true) + length(ic)))
_params = Flux.params(p)

prob = ODEProblem(model, [p[1], p[2], p[3], p[4], p[5]], time_interval, p[6:end])

function predict_rd() # Our 1-layer "neural network"
    sol = solve(remake(prob, u0 = [p[1], p[2], p[3], p[4], p[5]]), solver,
                p = p[6:end], saveat = tsteps) # override with new parameters
    return [sol[4, :] sol[5, :] sol[1, :] sol[2, :] .+ sol[3, :]]
end

function loss_rd()
    y_pred = predict_rd()
    y_true = [solution_true[4, :] solution_true[5, :] solution_true[1, :] solution_true[2,
                                                                                        :] .+
                                                                          solution_true[3,
                                                                                        :]]
    return sum(abs, y_true .- y_pred) / datasize
end # loss function
data = Iterators.repeated((), 3000)
opt = ADAM(0.05)
cb = function () #callback function to observe training
    display(loss_rd())
    # using `remake` to re-create our `prob` with current parameters `p`
    # display(plot(solve(remake(prob, p=p), AutoTsit5(Rosenbrock23()), saveat=tsteps), ylim=(0, 6)))
end

# Display the ODE with the initial parameter values.
cb()

Flux.train!(loss_rd, _params, data, opt, cb = cb)

println("Final Parameter Values: ")
println(p)
# plot(solution_true[4, :], label = "True Solution: w")
# plot!(predict_rd()[:, 1], label = "Predicted Solution: w")

# plot(solution_true[5, :], label = "True Solution: z")
# plot!(predict_rd()[:, 2], label = "Predicted Solution: z")

# plot(solution_true[1, :], label = "True Solution: x")
# plot!(predict_rd()[:, 3], label = "Predicted Solution: x")

# plot(solution_true[2, :] .+ solution_true[3, :], label = "True Solution: y + v")
# plot!(predict_rd()[:, 4], label = "Predicted Solution: y + v")
using Pkg
Pkg.activate(; temp = true)
Pkg.add("ModelingToolkit, DifferentialEquations, Gadfly, DiffEqParamEstim, Optim, RecursiveArrayTools")

using ModelingToolkit, DifferentialEquations, Gadfly
using DiffEqParamEstim, Optim, RecursiveArrayTools

@parameters mu bi bw al g dz k
@variables t s(t) i(t) w(t) r(t) y1(t) y2(t)
D = Differential(t)
solver = AutoTsit5(Rosenbrock23())

states = [s, i, w, r]
parameters = [mu, bi, bw, al, g, dz, k]
function cholera(du, u, p, t)
    mu, bi, bw, al, g, dz, k = p
    s, i, w, r = u
    du[1] = mu - bi * s * i - bw * s * w - mu * s + al * r
    du[2] = bw * s * w + bi * s * i - g * i - mu * i
    du[3] = dz * (i - w)
    du[4] = g * i - mu * r - al * r
end
@named model = ODESystem([
                             D(s) ~ mu - bi * s * i - bw * s * w - mu * s + al * r,
                             D(i) ~ bw * s * w + bi * s * i - g * i - mu * i,
                             D(w) ~ dz * (i - w),
                             D(r) ~ g * i - mu * r - al * r,
                         ], t, states, parameters)

ic = [1.0, -1.0, 1.0, -1.0]
time_interval = [0.0, 1.0]
datasize = 10
p_true = [1, 1.3, 1.1, 1.2, 1.1, 1, 0.5] # True Parameters
prob = ODEProblem(cholera, ic, time_interval, p_true)
sol = solve(prob, solver,
            saveat = range(time_interval[1], time_interval[2], length = datasize))

data_sample = Dict(1 => sol[2, :],
                   2 => sol[2, :] + sol[4, :] + sol[1, :])

p = plot(x = sol.t, y = data_sample[1], Geom.line, Guide.xlabel("t"), Guide.ylabel("data"))
push!(p,
      layer(x = sol.t, y = data_sample[2], Geom.line,
            Theme(default_color = colorant"red")))
τ = collect(range(time_interval[1], time_interval[2], length = datasize))
cost_function = build_loss_objective(prob, solver,
                                     L2Loss(τ, [data_sample[1]; data_sample[2]]),
                                     maxiters = 10000, verbose = false)

start = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
result = Optim.optimize(cost_function, start, BFGS())
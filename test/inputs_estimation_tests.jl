@testset "Estimation with known inputs" begin
    using ParameterEstimation
    using ModelingToolkit

    @parameters a b
    @independent_variables t
    @variables x1(t) x2(t) u(t) y1(t) y2(t)
    D = Differential(t)
    eqs = [D(x1) ~ -a * x1 + u, D(x2) ~ b * x1 - x2]
    @named model = System(eqs, t, [x1, x2], [a, b])
    outs = [y1 ~ x1, y2 ~ x2]
    inputs = [u ~ sin(t)]
    data = ParameterEstimation.sample_data(model, outs, [0.0, 1.0], [0.5, 2.0], [1.0, 0.5],
        21; inputs)

    res = ParameterEstimation.estimate(model, outs, data; inputs)
    @test !isempty(res)
    isempty(res) && return
    @test isapprox(res[1].parameters[a], 0.5, rtol = 1e-6)
    @test isapprox(res[1].parameters[b], 2.0, rtol = 1e-6)
    @test isapprox(res[1].states[x1], 1.0, rtol = 1e-6)
    @test isapprox(res[1].states[x2], 0.5, rtol = 1e-6)
end

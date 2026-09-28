@testset "Unicode parameter and state names" begin
    using ParameterEstimation
    using ModelingToolkit

    @parameters α β
    @variables t θ(t) x2(t) y1(t) y2(t)
    D = Differential(t)
    @named model = ODESystem([D(θ) ~ -α * θ, D(x2) ~ β * θ - x2], t, [θ, x2], [α, β])
    outs = [y1 ~ θ, y2 ~ x2]
    p_true = [0.5, 2.0]
    ic = [1.0, 0.5]
    data = ParameterEstimation.sample_data(model, outs, [0.0, 1.0], p_true, ic, 11)

    res = ParameterEstimation.estimate(model, outs, data)
    @test !isempty(res)
    isempty(res) && return
    @test isapprox(res[1].parameters[α], 0.5, rtol = 1e-6)
    @test isapprox(res[1].parameters[β], 2.0, rtol = 1e-6)
    @test isapprox(res[1].states[θ], 1.0, rtol = 1e-6)
    @test isapprox(res[1].states[x2], 0.5, rtol = 1e-6)
end

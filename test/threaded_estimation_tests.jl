@testset "Threaded and serial estimation agree" begin
    using ParameterEstimation
    using ModelingToolkit

    @parameters mu
    @variables t x1(t) y1(t)
    D = Differential(t)
    @named model = ODESystem([D(x1) ~ -mu * x1], t, [x1], [mu])
    outs = [y1 ~ x1 + x1^2]
    data = Dict{Any, Vector{Float64}}("t" => [0.0, 1 / 3, 2 / 3, 1.0],
                                      x1 + x1^2 => [2.0, 1.56301, 1.22995, 0.97441])

    for threaded in (false, true)
        res = ParameterEstimation.estimate(model, outs, data; threaded)
        @test Symbol(res[1].return_code) == :Success
        @test isapprox(res[1].parameters[mu], 0.5, atol = 1e-3)
        @test isapprox(res[1].states[x1], 1.0, atol = 1e-3)

        constrained = ParameterEstimation.estimate(model, outs, data; threaded,
                                                   parameter_constraints = Dict(mu => (0.6, 1.0)))
        @test isempty(constrained)
    end
end

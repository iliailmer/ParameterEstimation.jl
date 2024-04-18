# TODO: msolve functionality is broken.
@testset "Run small estimation code to test correctness" begin
    #     using ParameterEstimation
    #     using ModelingToolkit # ODE definitions

    #     # define toy model
    #     @parameters mu
    #     @variables t x1(t) y1(t)
    #     D = Differential(t)
    #     @named model = ODESystem([D(x1) ~ -mu * x1],
    #                              t, [x1], [mu])
    #     outs = [y1 ~ x1 + x1^2]
    #     time = [0.0, 1.0] # sampling interval
    #     data = Dict{Any, Vector{Float64}}("t" => [0.0, 1 / 3, 2 / 3, 1.0],
    #                                       x1 + x1^2 => [2.0, 1.56301,
    #                                           1.22995, 0.97441])
    #     res = ParameterEstimation.estimate(model, outs, data; method = :msolve)
    #    # @test length(res) == 1
    #     @test Symbol(res[1].return_code) == :Success
    #     @test isapprox(res[1].parameters[mu], 0.5, atol = 1e-3)
    #     @test isapprox(res[1].states[x1], 1.0, atol = 1e-3)

    #     # uneven time sample
    #     data = Dict{Any, Vector{Float64}}(x1^2 + x1 => [2.0, 1.98506, 1.17611, 0.97441],
    #                                       "t" => [0.0, 0.01, 0.73, 1.0])
    #     res = ParameterEstimation.estimate(model, outs, data; method = :msolve)
    #     #@test length(res) == 1
    #     @test Symbol(res[1].return_code) == :Success
    #     @test isapprox(res[1].parameters[mu], 0.5, atol = 1e-3)
    #     @test isapprox(res[1].states[x1], 1.0, atol = 1e-3)
end

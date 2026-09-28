@testset "Fixed seed gives reproducible results" begin
    using ParameterEstimation
    using ModelingToolkit
    using Random

    @parameters a b c
    @variables t x1(t) x2(t) y1(t)
    D = Differential(t)
    @named model = ODESystem([D(x1) ~ a * x1 - b * x1 * x2, D(x2) ~ -c * x2 + x1 * x2],
        t, [x1, x2], [a, b, c])
    outs = [y1 ~ x1]
    data = ParameterEstimation.sample_data(model, outs, [0.0, 1.0], [0.4, 0.8, 0.3],
        [1.0, 0.5], 21)

    Random.seed!(1)
    first_run = ParameterEstimation.estimate(model, outs, data)
    Random.seed!(2)
    second_run = ParameterEstimation.estimate(model, outs, data)
    @test length(first_run) == length(second_run)
    for (r1, r2) in zip(first_run, second_run)
        @test r1.parameters == r2.parameters
        @test r1.states == r2.states
    end

    Random.seed!(1)
    ParameterEstimation.estimate(model, outs, data)
    after_first = rand()
    ParameterEstimation.estimate(model, outs, data)
    after_second = rand()
    @test after_first != after_second
end

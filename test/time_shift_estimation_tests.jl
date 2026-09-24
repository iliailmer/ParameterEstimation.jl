@testset "Estimates do not depend on where the time interval starts" begin
    using ParameterEstimation
    using ModelingToolkit

    @parameters a b c d
    @variables t x1(t) x2(t) x3(t) x4(t) y1(t) y2(t) y3(t) y4(t)
    D = Differential(t)
    @named model = ODESystem([D(x1) ~ a + x2, D(x2) ~ b + x3, D(x3) ~ c + x4, D(x4) ~ d],
                             t, [x1, x2, x3, x4], [a, b, c, d])
    outs = [y1 ~ x1, y2 ~ x2, y3 ~ x3, y4 ~ x4]
    p_true = [2.0, 3.0, 4.0, 5.0]
    ic = [0.0, 0.0, 0.0, 0.0]

    for time_interval in ([0.0, 8.0], [-4.0, 4.0])
        data = ParameterEstimation.sample_data(model, outs, time_interval, p_true, ic, 9)
        res = ParameterEstimation.estimate(model, outs, data)
        @test !isempty(res)
        isempty(res) && continue
        for (p, v) in zip([a, b, c, d], p_true)
            @test isapprox(res[1].parameters[p], v, rtol = 1e-6)
        end
        for s in [x1, x2, x3, x4]
            @test isapprox(res[1].states[s], 0.0, atol = 1e-6)
        end
    end
end

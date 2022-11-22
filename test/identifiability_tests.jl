@testset "Check identifiability of `ODESystem` object" begin
    @parameters a01 a21 a12
    @variables t x0(t) x1(t) y1(t)
    D = Differential(t)

    eqs = [D(x0) ~ -(a01 + a21) * x0 + a12 * x1, D(x1) ~ a21 * x0 - a12 * x1]
    measured_quantities = [y1 ~ x0]
    de = ODESystem(eqs, t, name = :Test)
    correct = Dict{Any, Any}("nonidentifiable" => [a12, a21, substitute(x1, t => 0), a01],
                             "locally_not_globally" => Num[],
                             "globally" => [substitute(x0, t => 0)])
    output = get_identifiability(de; measured_quantities = measured_quantities)["identifiability"]

    @test isequal(Set(correct["nonidentifiable"]), Set(output["nonidentifiable"]))
    @test isequal(Set(correct["locally_not_globally"]), Set(output["locally_not_globally"]))
    @test isequal(Set(correct["globally"]), Set(output["globally"]))

    #--------------------------------------------------------------------------
    @parameters a01 a21 a12
    @variables t x0(t) x1(t) y1(t)
    D = Differential(t)

    eqs = [D(x0) ~ -(a01 + a21) * x0 + a12 * x1, D(x1) ~ a21 * x0 - a12 * x1]
    measured_quantities = [y1 ~ x0]

    de = ODESystem(eqs, t, name = :Test)
    correct = Dict{Any, Any}("nonidentifiable" => [a12, a21, substitute(x1, t => 0), a01],
                             "locally_not_globally" => Num[],
                             "globally" => [substitute(x0, t => 0)])
    output = get_identifiability(de; measured_quantities = measured_quantities)["identifiability"]

    @test isequal(Set(correct["nonidentifiable"]), Set(output["nonidentifiable"]))
    @test isequal(Set(correct["locally_not_globally"]), Set(output["locally_not_globally"]))
    @test isequal(Set(correct["globally"]), Set(output["globally"]))
end
@testset "Check size of transcendence_basis of polynomial systems" begin
    R, (a, b, c) = Nemo.polynomial_ring(Nemo.QQ, ["a", "b", "c"];
        internal_ordering = :degrevlex)
    Et = [a + b, a + b + c]
    vals = [1, 1, 1]
    _, tr_basis = ParameterEstimation.algebraic_independence(Et, [a, b, c], vals)
    @test length(tr_basis) == 1

    Et = [a + b, a^2 + b + c]
    vals = [1, 1, 1]
    _, tr_basis = ParameterEstimation.algebraic_independence(Et, [a, b, c], vals)
    @test length(tr_basis) == 1

    Et = [a + b + c]
    vals = [1, 1, 1]
    _, tr_basis = ParameterEstimation.algebraic_independence(Et, [a, b, c], vals)
    @test length(tr_basis) == 2

    Et = [a + b, a^2 + b + c, a^3 + b + c]
    vals = [1, 1, 1]
    _, tr_basis = ParameterEstimation.algebraic_independence(Et, [a, b, c], vals)
    @test length(tr_basis) == 0

    R, (a, b, c, d) = Nemo.polynomial_ring(Nemo.QQ, ["a", "b", "c", "d"];
        internal_ordering = :degrevlex)
    vals = [1, 1, 1, 1]
    Et = [a + d, b + c, a + a * b * c + a^2 * d]
    _, tr_basis = ParameterEstimation.algebraic_independence(Et, [a, b, d, c], vals)
    @test length(tr_basis) == 1
end

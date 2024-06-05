@testset "Convert Nemo polynomial systems into HomotopyContinuation type" begin
    R, (X, Y, Z) = Nemo.polynomial_ring(Nemo.QQ, ["x", "y", "z"])
    @var x y z
    Et = [X + Y + Z - 1, X + Y - 2, Z * X - 1]
    E = HomotopyContinuation.System([ParameterEstimation.nemo2hc(e) for e in Et])
    E_native = HomotopyContinuation.System([x + y + z - 1, x + y - 2, x * z - 1])
    @info "Source System: $Et"
    @info "Converted System: $E"
    @test E_native == E

    Et = [X^2 + 1 // 3 * Y + Z * X - 1, -X + Y^3 * X^2 - 2, Z * X - 1]
    E = HomotopyContinuation.System([ParameterEstimation.nemo2hc(e) for e in Et])
    E_native = HomotopyContinuation.System([x^2 + 1 // 3 * y + z * x - 1,
        -x + y^3 * x^2 - 2,
        x * z - 1])
    @info "Source System: $Et"
    @info "Converted System: $E"
    @test E_native == E

    Et = [X * Y * Z * X // 10 - 1, -X - Y - 2, Z * X - 1]
    E = HomotopyContinuation.System([ParameterEstimation.nemo2hc(e) for e in Et])
    E_native = HomotopyContinuation.System([
        x * y * z * x // 10 - 1,
        -x - y - 2,
        z * x - 1
    ])
    @info "Source System: $Et"
    @info "Converted System: $E"
    @test E_native == E
end

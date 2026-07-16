using Test, TestSetExtensions
using ModelingToolkit, SIAN, HomotopyContinuation
using Nemo
using ParameterEstimation

@info "Testing started"

@testset "All the tests" begin
    tests = isempty(ARGS) ?
            filter(f -> endswith(f, ".jl") && f != basename(@__FILE__), readdir(@__DIR__)) :
            [string(t, ".jl") for t in ARGS]
    for test in tests
        println(splitext(test)[1], ": ")
        include(test)
        println()
    end
end

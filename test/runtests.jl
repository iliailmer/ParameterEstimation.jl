using Test, TestSetExtensions
using ModelingToolkit, SIAN, HomotopyContinuation
using ParameterEstimation, Oscar

@info "Testing started"

@testset "All the tests" begin @includetests ARGS end

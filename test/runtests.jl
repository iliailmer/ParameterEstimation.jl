using Test, TestSetExtensions
using ModelingToolkit, Nemo, HomotopyContinuation
using ParameterEstimation

@info "Testing started"

@testset "All the tests" begin @includetests ARGS end

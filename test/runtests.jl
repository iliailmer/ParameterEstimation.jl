using Test, TestSetExtensions
using ModelingToolkit, SIAN, HomotopyContinuation
using Nemo
using ParameterEstimation

@info "Testing started"

@testset "All the tests" begin @includetests ARGS end

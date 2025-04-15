using TestParticles
using BifrostTools
using Test
using LinearAlgebra
using DifferentialEquations
using DataFrames
using Interpolations
using JLD2
using Random

if !isdefined(Main, :verbose)
    verbose = 3
end

include("testfields.jl")

if "fast_tests" in ARGS
    include("fast_tests.jl")
else
    @testset verbose = verbose ≥ 1 "All tests" begin
        include("fast_tests.jl")
        include("slow_tests.jl")
   end # testset all test
end
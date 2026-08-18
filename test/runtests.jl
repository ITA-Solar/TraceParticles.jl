using BifrostTools
using DataFrames
using DiffEqCallbacks
using Interpolations
using JLD2
using LinearAlgebra
using Logging
using OrdinaryDiffEq
using Random
using StaticArrays
using Statistics
using Test
using TraceParticles
import TraceParticles as TP

include("test_utils.jl")

if !isdefined(Main, :verbose)
    verbose = 3
end

if "debug" in ARGS
    loglevel = Logging.Debug
elseif "info" in ARGS
    loglevel = Logging.Info
else
    loglevel = Logging.Warn
end
logfilepath = joinpath(mktempdir(), "runtests.log.txt")
logger = TraceParticles.logger2fileandstderr(loglevel, logfilepath)

with_logger(logger) do
    if "fast_tests" in ARGS
        @testset verbose = verbose ≥ 1 "Tests" begin
            include("fast_tests.jl")
            @test isfile(logfilepath)
        end
    elseif "blank" in ARGS
        nothing
    else
        @testset verbose = verbose ≥ 1 "Tests" begin
            include("fast_tests.jl")
            include("slow_tests.jl")
            @test isfile(logfilepath)
        end # testset all test
    end
end

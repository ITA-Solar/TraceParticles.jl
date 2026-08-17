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

logfilepath = "runtests.log.txt"
if "debug" in ARGS
    loglevel = Logging.Debug
elseif "info" in ARGS
    loglevel = Logging.Info
else
    loglevel = Logging.Warn
end
logger = TraceParticles.logger2fileandstderr(loglevel, logfilepath)

if "fast_tests" in ARGS
    @testset verbose = verbose ≥ 1 "Tests" begin
        with_logger(logger) do
            include("fast_tests.jl")
        end
        targetfile = joinpath(@__DIR__,"fast_tests.log.txt")
        include("test_log.jl")
        testlog(loglevel, logfilepath, targetfile, verbose)
    end
elseif "blank" in ARGS
    nothing
else
    @testset verbose = verbose ≥ 1 "Tests" begin
        with_logger(logger) do
            include("fast_tests.jl")
            include("slow_tests.jl")
        end
        targetfile = joinpath(@__DIR__,"all_tests.log.txt")
        include("test_log.jl")
        testlog(loglevel, logfilepath, targetfile, verbose)
    end # testset all test
end

if loglevel == Logging.Info rm(logfilepath) end

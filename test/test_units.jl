#-------------------------------------------------------------------------------
# Created 10.01.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
# Script for running fast tests.
#-------------------------------------------------------------------------------

using Test

if !isdefined(Main, :verbose)
    verbose = 4
end

@testset verbose = verbose ≥ 2 "Fast tests" begin

    @testset verbose = verbose ≥ 3 "Unit tests" begin
        include("unit_tests/test_physics.jl")
        include("unit_tests/test_callbacks.jl")
        include("unit_tests/test_callback_specifiers.jl")
        include("unit_tests/test_traceparticlesparameters.jl")
        include("unit_tests/test_traceparticlesproblem.jl")
        include("unit_tests/test_rerun.jl")
    end
end # testset All test

#-------------------------------------------------------------------------------
# Created 10.01.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
# Script for running tests checked by Jenkins after it registrers activity on 
# repository.
#-------------------------------------------------------------------------------
#

using LinearAlgebra
using DifferentialEquations
using Interpolations


include("test_schemes.jl")
#include("test_interpolations_tp.jl")
#include("test_solvers.jl") # This file declear global variables whcich can't
# be overridden.

verbose = 4

@testset verbose = verbose ≥ 1 "Quick tests" begin
    @testset verbose = verbose ≥ 2 "Unit tests" begin
        #@testset verbose = verbose ≥ 3 "Interpolations_tp" begin
        #    testtrilinear(verbose ≥ 4)
        #    testlocateCell(verbose ≥ 4)
        #end
        @testset verbose = verbose ≥ 3 "Schemes" begin
            testeuler(verbose ≥ 4)
            testeulerCromer(verbose ≥ 4)
        end
        @testset verbose = verbose ≥ 3 "Solvers" begin
            #    testfullOrbit(verbose ≥ 4)
            #    testrelfullOrbitExplLeapFrog(verbose ≥ 4)
        end
    end # testset Unit tests

    @testset verbose = verbose ≥ 3 "I/O" begin
        #@testset verbose = true "tp" begin
        #    test_tp_save_and_tp_load(verbose >= 4)
        #end
    end

    @testset verbose = verbose ≥ 3 "Experiments" begin
        @testset verbose = true "ExB-drift" begin
            include("experiments/ExBdrift.jl")
        end # testset ExB-drift
        @testset verbose = true "∇B-drift" begin
            include("experiments/gradBdrift.jl")
        end # testset ∇B-drift
        @testset verbose = true "Mirroring" begin
            include("experiments/mirroring.jl")
        end # testset mirroring
        @testset verbose = true "Speiser (1965)" begin
            include("experiments/speiser1965.jl")
        end # testset Speiser 1965
    end # testset Experiments

end # testset All test

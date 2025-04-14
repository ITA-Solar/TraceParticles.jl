#-------------------------------------------------------------------------------
# Created 10.01.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
# Script for running tests checked by Jenkins after it registrers activity on 
# repository.
#-------------------------------------------------------------------------------
#

using TestParticles
using BifrostTools
using Test
using LinearAlgebra
using DifferentialEquations
using DataFrames
using Interpolations
using JLD2
using Random

include("testfields.jl")

if !isdefined(Main, :verbose)
    verbose = 3
end
@testset verbose = verbose ≥ 1 "Fast tests" begin

    @testset verbose = verbose ≥ 3 "Unit tests" begin
        include("test_physics.jl")
    end

    @testset verbose = verbose ≥ 3 "Regression tests" begin
        @testset verbose = verbose ≥ 4 "ExB-drift" begin
            include("experiments/ExBdrift.jl")
        end # testset ExB-drift
        @testset verbose = verbose ≥ 4 "∇B-drift" begin
            include("experiments/gradBdrift.jl")
        end # testset ∇B-drift
        @testset verbose = verbose ≥ 4 "Mirroring" begin
            include("experiments/mirroring.jl")
        end # testset mirroring
        @testset verbose = verbose ≥ 4 "Speiser (1965)" begin
            include("experiments/speiser1965.jl")
        end # testset Speiser 1965
        @testset verbose = verbose ≥ 4 "2D coronal jet" begin
            include("experiments/2Dcoronaljet/2Dcoronaljet.jl")
        end
    end # testset Experiments

end # testset All test
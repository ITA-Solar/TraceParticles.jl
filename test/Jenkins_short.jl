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


include("testfields.jl")

verbose = 4

@testset verbose = verbose ≥ 1 "Quick tests" begin

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

#-------------------------------------------------------------------------------
# Created 10.01.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
# Script for running tests checked by Jenkins after it registrers activity on 
# repository.
#-------------------------------------------------------------------------------

using TestParticles
using Test
using LinearAlgebra
using DifferentialEquations
using Interpolations

verbose = 4

include("testfields.jl")

include("Jenkins_short.jl")
@testset verbose = verbose ≥ 1 "Long tests" begin
    @testset verbose = verbose ≥ 3 "Experiments" begin
        @testset verbose = true "Dipole" begin
            include("experiments/dipoleloop.jl")
        end # testset dipole
    end
end # testset All test
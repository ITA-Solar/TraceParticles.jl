# Created 21.04.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 dipoleloop.jl
#
#-------------------------------------------------------------------------------
# Runs a series of magnetic dipole runs using 'dipole.jl' to reproduce fig. 19
# in Ripperda et al., 2018. The tests uses a previous result to compare
# with. The previous result was empirically validated by making a plot of
# 'phisNew' and comparing it with fig 19. 
#-------------------------------------------------------------------------------
using Interpolations
using LinearAlgebra
using OrdinaryDiffEq
using StaticArrays
using Test
using TraceParticles:
    get_fullorbit,
    get_guidingcentre!,
    guidingcentreapproximation!,
    lorentzforce!,
    magneticmoment

"""
    magneticdipolefield(x, y, z, M)
Return the magnetic field at (`x`, `y`, `z`) in a magnetic dipole field with
mangetic moment `M`.

B⃗(x, y, z) = 3αzx x̂  + 3αzy ŷ  + α(2z² - x² - y²)ẑ,
where
α = M/(x² + y² + z²)^(2.5).
"""
function magneticdipolefield(x, y, z, M)
    a = M / (x^2 + y^2 + z^2)^(5 / 2)
    return SVector(3a * z * x, 3a * z * y, a * (2z^2 - x^2 - y^2))
end
function T_dipole(q, m, M, R0, v, α)
    return @. 2π * q * M / (m * v^2 * R0) * (1 - 1 / 3 * sin(α)^0.62)
end # function T_dipole

include("../utils.jl")
include("dipole.jl")

qMmvec = collect(1:10) * 10.0
phis_fo = zeros(length(qMmvec))
phis_gca = zeros(length(qMmvec))
phis_analytical = zeros(length(qMmvec))

for i in eachindex(qMmvec)
    qMm = qMmvec[i]
    global ϕ
    global ϕ_FO
    global ϕ_GCA 
    ϕ, ϕ_FO, ϕ_GCA, _, _ = dipolesimulation(qMm)
    # The calculation of the angle phi is by evaluating atan(y/x).
    # Rotations more than π/2 degrees needs to be adjusted because atan(x/y)
    # jumps to -π/2 in the second circle quadrant.
    # Angles in (0, π/2) needs no addition
    # Angles in (π/2, 3π/2) needs +π 
    # Angles in (3/2, 5π/2) needs +2π
    # etc...
    # I use the analytical angle as reference for how much to add to the
    # numerical angles.
    n = div(ϕ, pi/2)
    extra_rads = (div(n, 2) + n % 2) * pi
    phis_fo[i] = ϕ_FO + extra_rads
    phis_gca[i] = ϕ_GCA + extra_rads
    phis_analytical[i] = ϕ
end

@testset verbose = true "GCA vs approx" begin
    @test all(@. isapprox(phis_gca[3:end], phis_analytical[3:end], rtol=1e-2))
end
@testset verbose = true "Full orbit vs approx:" begin
    @test all(@. isapprox(phis_fo[3:end], phis_analytical[3:end], rtol=1e-2))
end
@testset verbose = true "GCA vs Full orbit" begin
    @test all(@. isapprox(phis_fo[3:end], phis_gca[3:end], rtol=1e-2))
end

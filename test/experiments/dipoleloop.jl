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

qMmvec = collect(1:10)*10.0
phis = zeros(3, length(qMmvec))

for i = 1:length(qMmvec)
    global qMm = qMmvec[i]
    include("dipole.jl")
    phis[1,i] = ϕ # Analytical approximation
    phis[2,i] = ϕ_FO
    phis[3,i] = ϕ_GCA
end

# The calculation of the angle phi is by evaluating atan(y/x).
# Rotations more than π/2 degrees thus needs to be adjusted.
# The adjustments above are the ones that make the result most
# equal to Ripperda et al., 2018. All of them are pretty certain
# except maybe the first angle (qMm = 10) where the analytical
# approximation fail.
phisNew = copy(phis)
phisNew[3,4:end-1] .= phis[3,4:end-1] .+ pi
phisNew[2,4] = phis[2,4] + 2pi
phisNew[2:3,3] .= phis[2:3,3] .+ 2pi
phisNew[2:3,2] .= phis[2:3,2] .+ 3pi
phisNew[2:3,1] .= phis[2:3,1] .+ 4pi

@testset verbose = true "GCA vs approx" begin
    @test all(@. isapprox(phisNew[3,3:end], phisNew[1,3:end], rtol=1e-2))
end
@testset verbose = true "Full orbit vs approx:" begin
    @test all(@. isapprox(phisNew[2,3:end], phisNew[1,3:end], rtol=1e-2))
end
@testset verbose = true "GCA vs Full orbit" begin
    @test all(@. isapprox(phisNew[2,2:end], phisNew[3,2:end], rtol=1e-2))
end

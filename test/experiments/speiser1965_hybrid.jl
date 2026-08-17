# Created 13.09.25
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
"""
                speiser1965_hybrid.jl

A reproduction of the speiser test, using Adam Stanier's PhD thesis as
reference for parameter values, this time using the hybrid solver.

The particle should have a damped oscillation in x-direction (around zero)
while being accelerated in the z-direction. The damping should be proportional
to ∝ t^(-1/4), and the acceleration should be such that the z-position equal
 0.5q/m*Ez*t^2.

With a nonzero η the particle will eventually be ejected from its damped
oscillation. Should yield an eject time of πm/(qηb) (≈ 1.3e-4 sec with current
parameters).
"""

using OrdinaryDiffEq
using Interpolations
using StaticArrays
using Test
using TraceParticles
import TraceParticles as TP

include("../test_utils.jl")

#-------------------------------------------------------------------------------
#                            EXPERIMENTAL PARAMETERS
#
#-------------------------------------------------------------------------------
# NUMBER OF PARTICLES, SIMULATION DURATION, TIMESTEP   |
#......................................................|
numparticles = 1  # Number of particles to simulate    |
dt = 0.2e-8        # Time step [s]                     |
t0 = 0.0
tf = 2e-4 #n=10   # End time of simulation [s]         |
tspan = (t0, tf)   # timespan                           | 
#......................................................|

#...............................................
# PARTICLE TYPE
mass = TP.m_p
charge = TP.e

#...............................................
# INITIAL CONDITIONS
L0 = 1e4 # m
vel0 = [0.0, 0.0, 0.0] #478941.6577971818]
pos0 = [1e-6, 1e-10, 1e-10] * L0


#...............................................
# SPATIAL PARAMETERS (x, y, z)
numdims = 3
# Lower bounds of the three spatial axes
xi0 = (-3e-6 * L0, -0.2 * L0, -0.2 * L0)
# Upper bound of the three spatial axes
xif = (3e-6 * L0, 0.2 * L0, 0.2 * L0)
# Grid resolution of the axes
#n = (10, 10, 2)
ni = (100, 100, 100)

#...............................................
# MAGNETIC FIELD PARAMETERS
η = 0.025  # guide field?
d = 1e-4 # Current sheet width
b = 1e-2  # Characteristic field strength
if !isdefined(Main, :speiserBfield)
    function speiserBfield(
        x::Float64,
    )
        return [η, -x / d, 0.0] * b
    end
end

t_eject = π * mass / (charge * η * b)
tf = 1.1 * t_eject

#...............................................
# ELECTRIC FIELD PARAMETERS
v0 = 1e7
a = 1e-2 * v0 * b
Ez = -a

#-------------------------------------------------------------------------------
#                            RUN EXPERIMENT
#
#-------------------------------------------------------------------------------
# COMPUTING THE AXES, MAGNETIC FIELD AND ELECTRIC FIELD
xx, yy, zz, dx, dy, dz = create3Daxes(xi0, xif, ni)
Bfield = zeros(Float64, numdims, ni[1], ni[2], ni[3])
Efield = zeros(size(Bfield))
Efield[3, :, :, :] .= Ez
Bfield = stack(discretise(
    (x, y, z) -> speiserBfield(x),
    xx,
    yy,
    zz)
)
# Create interpolation objects
emfields = eachslice(vcat(Bfield, Efield), dims=(2, 3, 4))
emfields_itp = linear_interpolation((xx, yy, zz), emfields,
    extrapolation_bc=Flat()
)
emfields_itp = ElectromagneticFieldInterpolator(
    StaticInterpolation(emfields_itp)
)

#-------------------------------------------------------------------------------
# SIMULATION DURATION
#
numsteps = trunc(Int64, tf / dt)   # Number of timesteps in the simulation

#-------------------------------------------------------------------------------
# CREATE PROBLEM
# Set initial position and velocities
prob = ODEProblem(
    lorentzforce!,
    [pos0; vel0],
    tspan,
    (charge=charge, mass=mass, electromagneticfield=emfields_itp)
)
hprob = ODEProblem(
    hybridgcafo!,
    [pos0; vel0],
    tspan,
    HybridParamsWithDetection(
        charge=charge,
        mass=mass,
        electromagneticfield=emfields_itp,
        getphase=(integrator) -> π / 2,
    )
)
R, vparal, mu = get_guidingcentre(
    pos0,
    vel0,
    emfields_itp(pos0..., t0)[2],
    emfields_itp(pos0..., t0)[1],
    charge,
    mass
)
gcaprob = ODEProblem(
    guidingcentreapproximation!,
    [R[1], R[2], R[3], vparal],
    tspan,
    GCAParams(
        charge=charge,
        mass=mass,
        electromagneticfield=emfields_itp,
        magneticmoment=mu,
    )
)

hcb = DiscreteCallback(
    TP.HybridSwitchCondition(1e-2, 1e-2, Inf, 0.9),
    TP.hybridswitchaffect_withdetection!
)
#...............................................................................
# RUN SIMULATION
fosol = solve(prob)
hsol = solve(hprob, callback=hcb)
gcasol = solve(gcaprob)
#...............................................................................


#-------------------------------------------------------------------------------
# TESTING
times = collect(range(0.0, step=dt, length=numsteps + 1))
indx = @. isapprox(fosol.t, t_eject, rtol=2e-5)
velysimτ = fosol.u[indx][1][5]
velyτ = -0.5 * 3.0 * a / (η * b)
hindx = @. isapprox(hsol.t, t_eject, rtol=3e-5)
hvelysimτ = hsol.u[hindx][1][5]
foendpos = fosol.u[end][1:3]
hsolendR = hsol.u[end][1:3]
foendR, foendvparal, foendmu = get_guidingcentre(
    foendpos,
    fosol.u[end][4:6],
    emfields_itp(foendpos..., tf)[2],
    emfields_itp(foendpos..., tf)[1],
    charge,
    mass
)
foef = get_observable(fosol, :energy, tf, EoM="FO")
hef = get_observable(hsol, :energy, tf)

@testset "Speiser 1965 - Hybrid solver" begin
    @testset verbose = true "Full orbit: Ejection-time" begin
        @test isapprox(velysimτ, velyτ, atol=abs(5 * velyτ))
    end # testset Full orbit: RK4

    @testset verbose = true "Hybrid: Ejection time and Nof. switches" begin
        @test isapprox(hvelysimτ, velyτ, atol=abs(5 * velyτ))
        @test hsol.prob.p.nswitches == 3
    end
    @testset verbose = true "Hybrid vs full orbit" begin
        @test isapprox(hsolendR[1], foendR[1], rtol=1e-3)
        @test isapprox(hsolendR[2], foendR[2], rtol=1e-5)
        @test isapprox(hsolendR[3], foendR[3], rtol=1e-4)
        @test isapprox(hsol.u[end][4], foendvparal, rtol=1e-2)
        @test isapprox(hsol.prob.p.magneticmoment, foendmu, rtol=5)
        @test isapprox(hef, foef, rtol=1e-5)
    end
end

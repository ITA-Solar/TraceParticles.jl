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
using TraceParticles:
    EoMID,
    get_guidingcentre,
    get_observable,
    guidingcentreapproximation!,
    GCAParams,
    hybridgcafo!,
    HybridSwitchCondition,
    HybridParamsWithDetection,
    hybridswitchaffect_withdetection!,
    lorentzforce!,
    proton_mass,
    proton_charge

include("../utils.jl")

# MAGNETIC FIELD PARAMETERS
η = 0.025 # Field component perpendicular to sheet
d = 1e-4 # Current sheet width
b = 1e-2  # Characteristic field strength
function speiserBfield(x)
    return SVector(b*η, -x * b / d, 0.0)
end

# ELECTRIC FIELD PARAMETERS
v0 = 1e7
a = 1e-2 * v0 * b
Ez = -a

# PARTICLE Parameters
mass = proton_mass
charge = proton_charge
# Speiser's analytical ejection time
t_eject = π * mass / (charge * η * b)
# End time of simulation [s]
tf = 1.1 * t_eject
tspan = (0, tf)   # timespan

# INITIAL CONDITIONS
L0 = 1e4 # m
vel0 = [0.0, 0.0, 0.0] #478941.6577971818]
pos0 = [1e-6, 1e-10, 1e-10] * L0

# SPATIAL PARAMETERS (x, y, z)
numdims = 3
# Lower bounds of the three spatial axes
xi0 = (-3e-6 * L0, -0.2 * L0, -0.2 * L0)
# Upper bound of the three spatial axes
xif = (3e-6 * L0, 0.2 * L0, 0.2 * L0)
# Grid resolution of the axes
ni = (100, 100, 100)

# COMPUTING THE AXES, MAGNETIC FIELD AND ELECTRIC FIELD
xx, yy, zz, dx, dy, dz = create3Daxes(xi0, xif, ni)
bfield = discretise3D((x, y, z) -> speiserBfield(x), xx, yy,zz)
# Create interpolation object
bfield_itp = linear_interpolation((xx, yy, zz), bfield, extrapolation_bc=Flat())
electromagneticfield(x, y, z, _) = SVector(0,0,Ez), bfield_itp(x, y, z)

# CREATE PROBLEM
# Set initial position and velocities
prob = ODEProblem(
    lorentzforce!,
    [pos0; vel0],
    tspan,
    (charge=charge, mass=mass, electromagneticfield=electromagneticfield)
)
hprob = ODEProblem(
    hybridgcafo!,
    [pos0; vel0],
    tspan,
    HybridParamsWithDetection(
        charge=charge,
        mass=mass,
        electromagneticfield=electromagneticfield,
        getphase=(integrator) -> π / 2,
        initialeomid=EoMID.FullOrbit
    )
)
R, vparal, mu = get_guidingcentre(
    pos0,
    vel0,
    electromagneticfield(pos0..., 0)[2],
    electromagneticfield(pos0..., 0)[1],
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
        electromagneticfield=electromagneticfield,
        magneticmoment=mu,
    )
)

hcb = DiscreteCallback(
    HybridSwitchCondition(1e-2, 1e-2, Inf, 0.9),
    hybridswitchaffect_withdetection!
)

# RUN SIMULATION
fosol = solve(prob, Vern9(); reltol=1e-7)
hsol = solve(hprob, Vern9(lazy=false); reltol=1e-7, callback=hcb)
gcasol = solve(gcaprob)

# TESTING
velyτ = -0.5 * 3.0 * a / (η * b)
velysimτ = fosol(t_eject)[5]
hvelysimτ = hsol(t_eject)[5]
foendpos = fosol.u[end][1:3]
hsolendR = hsol.u[end][1:3]
foendR, foendvparal, foendmu = get_guidingcentre(
    foendpos,
    fosol.u[end][4:6],
    electromagneticfield(foendpos..., tf)[2],
    electromagneticfield(foendpos..., tf)[1],
    charge,
    mass
)
foef = get_observable(fosol, :energy, tf, EoM="FO")
hef = get_observable(hsol, :energy, tf)

@testset "Speiser 1965 - Hybrid solver" begin
    @testset verbose = true "Full orbit: Ejection-time" begin
        @test isapprox(velysimτ, velyτ, atol=abs(1 * velyτ))
    end # testset Full orbit: RK4

    @testset verbose = true "Hybrid: Ejection time and Nof. switches" begin
        @test isapprox(hvelysimτ, velyτ, atol=abs(1 * velyτ))
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

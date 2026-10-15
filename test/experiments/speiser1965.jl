# Created 21.04.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
"""
                speiser.jl

A reproduction of the speiser test, using Adam Stanier's PhD thesis as
reference  for parameter values.

The particle should have a damped oscillation in x-direction (around zero)
while being accelerated in the z-direction. The damping should be proportional
to ∝ t^(-1/4), and the acceleration should be such that the z-position equal
 0.5q/m*Ez*t^2.

With a nonzero η the particle will eventually be ejected from its damped
oscillation. Should yield an eject time of πm/(qηb) (≈ 1.3e-4 sec with current
parameters).
"""

using Interpolations
using OrdinaryDiffEq
using StaticArrays
using Test
using TraceParticles:
    proton_mass,
    proton_charge,
    lorentzforce!

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
prob = ODEProblem(
    lorentzforce!,
    [pos0; vel0],
    tspan,
    (charge=charge, mass=mass, electromagneticfield=electromagneticfield)
)

# RUN SIMULATION
sol = solve(prob, Vern9(), reltol=1e-7)

# TESTING
velysimτ = sol(t_eject)[5]
velyτ = -0.5 * 3.0 * a / (η * b)
@testset verbose = true "Y-velocity at ejection" begin
    @test isapprox(velysimτ, velyτ, atol=abs(1 * velyτ))
end

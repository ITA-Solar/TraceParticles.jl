# Created 31.03.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
"""
    gradBdrift.jl

This experiment tests the ∇B-drift term in the GCA against the full orbit
solution in a static magnetic field where the gradient is perpendicular to the
magnetic field direction. The parameters yield a Larmor radius that is much
smaller than the characteristic length scale of the magnetic field
gradient, such that the GCA is valid. The solutions is also compared to the
approximate gradient drift presented in e.g. Chen (2016) or Aschwanden (2006).

To be able to compare the GCA and full orbit, the simulation duration has to
be a multiple of the gyration period.
"""

#-------------------------------------------------------------------------------
#                            EXPERIMENTAL PARAMETERS
#
#-------------------------------------------------------------------------------
# NUMBER OF PARTICLES, SIMULATION DURATION, TIMESTEP   |
#......................................................|
numparticles = 1  # Number of particles to simulate    |
dt = 0.01         # Time step [s]                      |
# Use a final time equal to 10 gyrations for the full  |
# orbit, according to the parameters set (mass, charge,|
# vel0, pos0, B0, a)                                   |
tf = 10 * 1 / 1.7507  # End time of simulation [s]     |
tspan = (0.0, tf)     # simulation time span           |
#tf = 2.3 #n=10   # End time of simulation [s]         |
#......................................................|

#...............................................
# PARTICLE TYPE
mass = 1
charge = 1

#...............................................
# INITIAL CONDITIONS
pos0 = [0.5, 0.5, 0.5]
vel0 = [0.0, 0.1, 0.0]

#...............................................
# SPATIAL PARAMETERS (x, y, z)
numdims = 3
# Lower bounds of the three spatial axes
xi0 = (0.0, 0.0, 0.0)
# Upper bound of the three spatial axes
xif = (1.0, 1.0, 1.0)
# Grid resolution of the axes
#n = (10, 10, 2)
ni = (100, 100, 2)

#...............................................
# MAGNETIC FIELD PARAMETERS
a = 10.0 # gradient in magnetic field
B0 = 6.0 # Additional constant
function gradBfield(
    y::Float64,
)
    return [0.0, 0.0, a * y + B0]
end

#...............................................
# ELECTRIC FIELD PARAMETERS
# It will not exist
Ex = 0.0
Ey = 0.0

#...............................................
# SOLVER CONDITIONS
# RK4 for GCA? Euler-Cromer for full orbit?

#-------------------------------------------------------------------------------
#                            RUN EXPERIMENT
#
#-------------------------------------------------------------------------------
# COMPUTING THE AXES, MAGNETIC FIELD AND ELECTRIC FIELD
xx, yy, zz, dx, dy, dz = create3Daxes(xi0, xif, ni)
Bfield = zeros(Float64, numdims, ni[1], ni[2], ni[3])
Efield = zeros(size(Bfield))
Efield[1, :, :, :] .= Ex
Efield[2, :, :, :] .= Ey
#Efield[1,:,1:30,:] .= 0.0
Bfield = stack(discretise((x, y, z) -> gradBfield(y), xx, yy, zz))
# Create interpolation objects
emfields = eachslice(vcat(Bfield, Efield), dims=(2, 3, 4))
emfields_itp = linear_interpolation((xx, yy, zz), emfields,
    extrapolation_bc=Flat()
)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# SIMULATION DURATION
#
numsteps = trunc(Int64, tf / dt)   # Number of timesteps in the simulation
#println("Number of time steps = $numsteps.")

#-------------------------------------------------------------------------------
# PARTICLE CREATION
# Set initial position and velocities
#
B⃗ = emfields_itp(pos0...)[1:3]
E⃗ = zeros(3)
B = norm(B⃗)
b̂ = B⃗ / B
v = norm(vel0)
vparal = vel0 ⋅ b̂
vperp = √(v^2 - vparal^2)
μ = mass * vperp^2 / (2B)

#-------------------------------------------------------------------------------
# CREATE PROBLEM
# Non-relativisitc Euler-Cromer ?
prob_FO = ODEProblem(
    lorentzforce!,
    [pos0; vel0],
    tspan,
    (charge=charge, mass=mass, fields=emfields_itp),
)
prob_GCA = ODEProblem(
    guidingcentreapproximation!,
    [pos0; vparal],
    tspan,
    (charge=charge, mass=mass, magneticmoment=μ, fields=emfields_itp)
)

#...............................................................................
# RUN SIMULATION
sol_FO = DifferentialEquations.solve(prob_FO)
sol_GCA = DifferentialEquations.solve(prob_GCA)

#...............................................................................


#-------------------------------------------------------------------------------
# TESTING
# Initial gyrofrequency
B = norm(gradBfield(pos0[2]))
vperp = √(vel0[1]^2 + vel0[2]^2)
ff = charge * B / (mass * 2π)
rL = mass * vperp / (charge * B)
L = a / B
vdrift = [-a * mass * vel0[2]^2 / (2charge * B^2), 0.0, 0.0]
posf_anal = tf * vdrift + pos0
@testset verbose = true "GCA: Defualt alg." begin
    @test isapprox(posf_anal, sol_GCA.u[end][1:3], rtol=0.0001)
    @test isapprox(sol_FO.u[end][1:3], sol_GCA.u[end][1:3], rtol=0.001)
end # testset GCA: Euler
@testset verbose = true "Full orbit: Defualt alg." begin
    @test isapprox(posf_anal, sol_FO.u[end][1:3], rtol=0.001)
end # testset GCA: Euler


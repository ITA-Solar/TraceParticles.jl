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

#-------------------------------------------------------------------------------
#                            EXPERIMENTAL PARAMETERS
#
#-------------------------------------------------------------------------------
# NUMBER OF PARTICLES, SIMULATION DURATION, TIMESTEP   |
#......................................................|
numparticles = 1  # Number of particles to simulate    |
dt = 0.2e-8        # Time step [s]                     |
tf = 2e-4 #n=10   # End time of simulation [s]         |
tspan = (0, tf)   # timespan                           | 
#......................................................|

#...............................................
# PARTICLE TYPE
mass = tp.m_p
charge = tp.e

#...............................................
# INITIAL CONDITIONS
L0 = 1e4 # m
vel0 = [0.0, 0.0, 0.0] #478941.6577971818]
pos0 = [1e-6, 1e-10, 1e-10]*L0


#...............................................
# SPATIAL PARAMETERS (x, y, z)
numdims = 3
# Lower bounds of the three spatial axes
xi0 = (-3e-6*L0, -0.2*L0, -0.2*L0)
# Upper bound of the three spatial axes
xif = (3e-6*L0, 0.2*L0, 0.2*L0)
# Grid resolution of the axes
#n = (10, 10, 2)
ni = (100, 100, 100)

#...............................................
# MAGNETIC FIELD PARAMETERS
η = 0.025  # guide field?
d = 1e-4 # Current sheet width
b = 1e-2  # Characteristic field strength
function speiserBfield(
    x::Float64,
    )
    return [η, -x/d, 0.0]*b
end

t_eject = π*mass/(charge*η*b)
tf = 1.1*t_eject

#...............................................
# ELECTRIC FIELD PARAMETERS
v0 = 1e7 
a = 1e-2 * v0*b
Ez = -a

#-------------------------------------------------------------------------------
#                            RUN EXPERIMENT
#
#-------------------------------------------------------------------------------
# COMPUTING THE AXES, MAGNETIC FIELD AND ELECTRIC FIELD
xx, yy, zz, dx, dy, dz = create3Daxes(xi0, xif, ni)
Bfield = zeros(Float64, numdims, ni[1], ni[2], ni[3])
Efield = zeros(size(Bfield))
Efield[3, :,:,:] .= Ez
Bfield = stack(discretise(
    (x,y,z) -> speiserBfield(x),
    xx,
    yy,
    zz)
)
# Create interpolation objects
emfields = eachslice(vcat(Bfield, Efield), dims=(2,3,4))
emfields_itp = linear_interpolation((xx, yy, zz), emfields, 
    extrapolation_bc=Flat()
    )

#-------------------------------------------------------------------------------
# SIMULATION DURATION
#
numsteps = trunc(Int64, tf/dt)   # Number of timesteps in the simulation

#-------------------------------------------------------------------------------
# CREATE PROBLEM
# Set initial position and velocities
prob = ODEProblem(
    lorentzforce!,
    [pos0; vel0],
    tspan,
    (charge, mass, emfields_itp)
    )

#...............................................................................
# RUN SIMULATION
sol = DifferentialEquations.solve(prob)
#...............................................................................


#-------------------------------------------------------------------------------
# TESTING
times = collect(range(0.0, step=dt, length=numsteps+1))

indx = @. isapprox(sol.t, t_eject, rtol=2e-5)
velysimτ = sol.u[indx][1][5]
velyτ = -0.5*3.0*a/(η*b)
@testset verbose = true "Full orbit: Default alg." begin
    @test isapprox(velysimτ, velyτ, atol=abs(5*velyτ))
end # testset Full orbit: RK4


#-------------------------------------------------------------------------------
# Created 19.12.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
"""
    testExBdrift.jl

    This experiment tests numerical particle paths produced by the tp-code 
against analytical solutions.
    An electron and proton in a static and homogenous electromagnetic field.
The electric field is normal to the magnetic field such that the particles will 
experience an ExB-drift. The resulting numerical paths are compared to the
analytic paths, which root mean square errors entails the test assertions. 

    In practice, this experiment also tests dependent type-constructors and
methods, such as:
    ParticleSoA(pos::Matrix, vel::Matrix, species::Vector, numSteps::Integer)
    Mesh(bField::Array{4},
         eField::Array{4},
         xx::Vector,
         yy::Vector,
         zz::Vector)
    run!(patch:Patch)
and the chosen numerical solver, scheme and interpolation chosen.
"""

using OrdinaryDiffEq
using Interpolations
using Test
using TraceParticles


#-------------------------------------------------------------------------------
# EXPERIMENT PARAMTERS
#
numDims = 3       # Number of spatial dimensions
numParticles = 2  # Number of particles to simulae
#   (electron, proton)
tf = 2 * 2π / (5.93096958e7) # s. End time of simulation
tspan = (0, tf)
#   multiple of electron (1eV) gyrofreq. at Bz=3.3e-4 G
numSteps = 1000   # Number of timesteps in the simulation
dt = tf / numSteps  # Size of timestep
# The electromagnetic field is static and homogeneous. Hence the mesh only 
# needs one cell, i.e. the boundary.
Bx = 0.0  # magnetic field component only in the positive z-direction
By = 0.0
Bz = 3.37213e-4 # G (so that electron larmor radius is 1 cm)
Ex = 0.0
Ey = 0.02Bz / tf  # Electric field component only in the positive y-direction
Ez = 0.0

# Particle conditions
qe = -TraceParticles.e  # Electron charge
qp = TraceParticles.e   # Proton charge
me = TraceParticles.m_e # Electron mass
mp = TraceParticles.m_p # Proton mass
# Set initial position and velocities
vdriftx = Ey / Bz
# Electron initial velocity
vxe = 0.0#593_096.95848 # m/s (1eV electron
vye = 0.0
vze = 0.0
# Proton initial velocity
vxp = 13_841.122177 # m/s (1eV proton)
vyp = 0.0
vzp = 0.0
# Electron initial position
x0e = 0.0
y0e = 0.0
z0e = 0.0
# Proton initial position
x0p = 0.0
y0p = 0.0
z0p = 0.0

#-------------------------------------------------------------------------------
# ANALYTICAL SOLUTION 
times = LinRange(0, tf, numSteps + 1)

# Gyrofrequency
function ω(q, B, m)
    return q * B / m
end
# Larmor radius
function larmorRadius(ω, vperp)
    return vperp / ω
end
# x-position
function x(t, x0, ωt, rL, Ey, Bz)
    return @. rL * sin(ωt * t) + Ey / Bz * t
end
# y-positoin
function y(t, y0, ωt, rL, Ex, Bz)
    return @. -rL + rL * cos(ωt * t) - Ex / Bz * t
end
# z-position
function z(t, z0, vz0, ω, Ez, Bz)
    return @. z0 + 0.5ω / Bz * Ez * t^2 + vz0 * t
end

# Define gyrofrequency and Larmor radius for both particles
ωe = ω(qe, Bz, me)
ωp = ω(qp, Bz, mp)
vperpe = vxe - vdriftx # Subtract drift velocity to find vperp
vperpp = vxp - vdriftx
rLe = larmorRadius(ωe, vperpe)
rLp = larmorRadius(ωp, vperpp)

# Compute time evolution of position of both particles
pose = zeros(Float64, numDims, numSteps + 1)
posp = zeros(Float64, numDims, numSteps + 1)
pose[1, :] .= x(times, x0e, ωe, rLe, Ey, Bz)
pose[2, :] .= y(times, y0e, ωe, rLe, Ex, Bz)
pose[3, :] .= z(times, z0e, vze, ωe, Ez, Bz)
posp[1, :] .= x(times, x0p, ωp, rLp, Ey, Bz)
posp[2, :] .= y(times, y0p, ωp, rLp, Ex, Bz)
posp[3, :] .= z(times, z0p, vzp, ωp, Ez, Bz)

analytic_ex = linear_interpolation(times, x(times, x0e, ωe, rLe, Ey, Bz))
analytic_ey = linear_interpolation(times, y(times, y0e, ωe, rLe, Ex, Bz))
analytic_ez = linear_interpolation(times, z(times, z0e, vze, ωe, Ez, Bz))
analytic_px = linear_interpolation(times, x(times, x0p, ωp, rLp, Ey, Bz))
analytic_py = linear_interpolation(times, y(times, y0p, ωp, rLp, Ex, Bz))
analytic_pz = linear_interpolation(times, z(times, z0p, vzp, ωp, Ez, Bz))

#-------------------------------------------------------------------------------
# NUMERICAL SOLUTION 
#-------------------------------------------------------------------------------
# MESH CREATION
B = zeros(Float64, numDims, 2, 2, 2)
E = zeros(Float64, numDims, 2, 2, 2)
B[3, :, :, :] .= Bz
E[1, :, :, :] .= Ex
E[2, :, :, :] .= Ey
E[3, :, :, :] .= Ez
# Create coordinates of cell
xx = [-0.1, 1.1]
yy = [-0.1, 1.1]
zz = [-0.1, 1.1]
# Create interpolation objects
emfields = eachslice(vcat(B, E), dims=(2, 3, 4))
emfields_itp = linear_interpolation((xx, yy, zz), emfields,
    extrapolation_bc=Flat()
)
emfields_itp = ElectromagneticFieldInterpolator(
    StaticInterpolation(emfields_itp)
)

#-------------------------------------------------------------------------------
# PARTICLE CREATION
# Full orbit particles
# Set initial position and velocities
pos = zeros(Float64, numDims, numParticles) # Position -> origin
vel = zeros(Float64, numDims, numParticles) # Velocity
#  Create DifferentialEquations.jl particle type
ic_FO = [[x0e, y0e, z0e, vxe, vye, vze], [x0p, y0p, z0p, vxp, vyp, vzp]]
params_FO = (
    [-TraceParticles.e, TraceParticles.e],
    [TraceParticles.m_e, TraceParticles.m_p],
    emfields_itp
)

# GCA particles
# Set initial position and velocities
R = zeros(Float64, numDims, numParticles) # Position -> origin
velGCA = zeros(Float64, numDims + 3, numParticles) # Velocity
vperpe = √(vxe^2 + vye^2)
vperpp = √(vxp^2 + vyp^2)
# Magnetic moments of particles
μₑ = me * vperpe^2 / 2Bz
μₚ = mp * vperpp^2 / 2Bz

#  Create initial condition vector  and params
ic_GCA = [[x0e, y0e, z0e, vze], [x0p, y0p, z0p, vzp]]
params_GCA = (
    [-TraceParticles.e, TraceParticles.e],
    [TraceParticles.m_e, TraceParticles.m_p],
    [μₑ, μₚ],
    emfields_itp
)

#-------------------------------------------------------------------------------
# PUSHERS TO TEST:
# Non-relativisitc Euler-Cromer
# Relativistic Vay pusher
# Relativistic Boris pusher
# non-relativistic GCA

#-------------------------------------------------------------------------------
# RUN SIMULATION

# Problem functions
prob_func_FO(prob, ctx) = remake(prob, u0=ic_FO[ctx.sim_id], p=(
    charge=params_FO[1][ctx.sim_id],
    mass=params_FO[2][ctx.sim_id],
    electromagneticfield=params_FO[3],
))
prob_func_GCA(prob, ctx) = remake(prob, u0=ic_GCA[ctx.sim_id], p=(
    charge=params_GCA[1][ctx.sim_id],
    mass=params_GCA[2][ctx.sim_id],
    magneticmoment=params_GCA[3][ctx.sim_id],
    electromagneticfield=params_GCA[4],
))

# ODEProblems
prob_FO = ODEProblem(lorentzforce!, zeros(6), tspan)
prob_GCA = ODEProblem(guidingcentreapproximation!, zeros(4), tspan)
# EnsembleProblems
eprob_FO = EnsembleProblem(prob_FO, prob_func=prob_func_FO)
eprob_GCA = EnsembleProblem(prob_GCA, prob_func=prob_func_GCA)

# SOLVE
sol_FO = solve(eprob_FO, Tsit5, reltol=1.1e-4, trajectories=2)
sol_GCA = solve(eprob_GCA, Tsit5, reltol=1.1e-4, trajectories=2)

#-------------------------------------------------------------------------------
# CALCULATE DEVIATIONS FROM ANALYTICAL SOLUTION
# Calculate the root mean squared error between numerical and analytical 
#   trajectories for both electron and proton for all position components.

# Euler-Cromer
# Electron
rmse_ex = 0.0
rmse_ey = 0.0
rmse_ez = 0.0
rmse_px = 0.0
rmse_py = 0.0
rmse_pz = 0.0
nsteps_e = length(sol_FO.u[1].t)
nsteps_p = length(sol_FO.u[2].t)
for i = 1:nsteps_e
    global rmse_ex += (sol_FO.u[1].u[i][1] - analytic_ex(sol_FO.u[1].t[i]))^2
    global rmse_ey += (sol_FO.u[1].u[i][2] - analytic_ey(sol_FO.u[1].t[i]))^2
    global rmse_ez += (sol_FO.u[1].u[i][3] - analytic_ez(sol_FO.u[1].t[i]))^2
end
for i = 1:nsteps_p
    global rmse_px += (sol_FO.u[2].u[i][1] - analytic_px(sol_FO.u[2].t[i]))^2
    global rmse_py += (sol_FO.u[2].u[i][2] - analytic_py(sol_FO.u[2].t[i]))^2
    global rmse_pz += (sol_FO.u[2].u[i][3] - analytic_pz(sol_FO.u[2].t[i]))^2
end
rmse_ex = sqrt(rmse_ex / nsteps_e)
rmse_ey = sqrt(rmse_ey / nsteps_e)
rmse_ez = sqrt(rmse_ez / nsteps_e)
rmse_px = sqrt(rmse_px / nsteps_p)
rmse_py = sqrt(rmse_py / nsteps_p)
rmse_pz = sqrt(rmse_pz / nsteps_p)

#-------------------------------------------------------------------------------
# TEST RESULTS
@testset verbose = true "Full orbit: Default alg." begin
    @test isapprox(rmse_ex, 1e-7, rtol=1e-1)
    @test isapprox(rmse_ey, 1e-7, rtol=1e-1)
    @test isapprox(rmse_ez, 0.0, rtol=1e-1)
    @test isapprox(rmse_px, 5e-15, rtol=1e-1)
    @test isapprox(rmse_py, 1e-11, rtol=1e-0)
    @test isapprox(rmse_pz, 0.0, rtol=1e-6)
end # testset Vay
# GCA using DifferentialEquations.jl
@testset verbose = true "GCA: Default alg." begin
    # Electron
    @test isapprox(sol_GCA.u[1].u[2][1], 0.02, rtol=1e-13)
    @test isapprox(sol_GCA.u[1].u[2][2], 0.00, atol=1e-15)
    @test isapprox(sol_GCA.u[1].u[2][3], 0.00, atol=1e-15)
end



using OrdinaryDiffEq
using Interpolations
using LinearAlgebra
using Random
using StaticArrays
using Statistics
using Test
using TraceParticles
import DiffEqCallbacks as CB

#-------------------------------------------------------------------------------
# EXPERIMENT PARAMTERS
#
numDims = 3       # Number of spatial dimensions
numParticles = 1  # Number of particles to simulae
#   (electron, proton)
#tf = 2*2π/(5.93096958e7) # s. End time of simulation
tf = 10 * 2 * 2π / (5.93096958e7) # s. End time of simulation
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
#Ey = 0.02Bz/tf  # Electric field component only in the positive y-direction
Ey = 0.02Bz / (tf / 10)
Ez = 0.0

# Particle conditions
q = -TraceParticles.e  # Electron charge
m = TraceParticles.m_e # Electron mass
# Set initial position and velocities
Evec = [Ex, Ey, Ez]
Bvec = [Bx, By, Bz]
vdrift = cross(Evec, Bvec) / norm(Bvec)^2
# Electron initial velocity
vx = 0.0#593_096.95848 # m/s (1eV electron
vy = 0.0
vz = 0.0
vel = [vx, vy, vz]
# Electron initial position
x0 = 0.0
y0 = 0.0
z0 = 0.0

#-------------------------------------------------------------------------------
# ANALYTICAL SOLUTION 
times = LinRange(0, tf, numSteps + 1)

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
# Particle starting with GCA needs magnetic moment
gyrationvel = vel .- vdrift
vperp = norm(gyrationvel)
vparal = dot(vel, Bvec)
# Magnetic moments of particles
μ = m * vperp^2 / (2norm(Bvec))
q = -TraceParticles.e  # Electron charge
m = TraceParticles.m_e # Electron mass
# Set initial position and velocities
E, B = emfields_itp(x0, y0, z0)
R, vparal, mu = get_guidingcentre([x0, y0, z0], [vx, vy, vz], B, E, q, m)
ic = [
    [x0, y0, z0, vx, vy, vz],
    [R..., vparal, 0, 0],
    [x0, y0, z0, vx, vy, vz],
    [R..., vparal, 0, 0],
]
eoms = [
    hybridgcafo!,
    hybridgcafo!,
    lorentzforce!,
    guidingcentreapproximation!,
]
params = (
    charge=q,
    mass=m,
    electromagneticfield=emfields_itp,
    magneticmoment=ones(2) * μ,
    rng=[Xoshiro(1), Xoshiro(2)],
    getphase=(integrator) -> π / 2,
    initialeomid=[
                  TraceParticles.EoMID.FullOrbit,
                  TraceParticles.EoMID.GuidingCentreApproximation,
                ]
)

#-------------------------------------------------------------------------------
# CALLBACKS
ω_c = q * norm(B) / m
freq = ω_c / (2π)  # Gyrofrequency in Hz
period = abs(1 / freq)  # Gyroperiod in seconds
cb = CB.PresetTimeCallback(
    [
        4period,
        8period,
        12period,
        16period,
    ],
    hybridswitchaffect!
)
#cb = CB.PresetTimeCallback(
#    [
#        1 / 6 * tf,
#        2 / 6 * tf,
#        3 / 6 * tf,
#        4 / 6 * tf,
#        5 / 6 * tf,
#    ],
#    TraceParticles.switch_affect!
#)
#-------------------------------------------------------------------------------
# RUN SIMULATION


# Problem functions
prob_func(prob, ctx) = begin
    remake(
        prob,
        f=eoms[ctx.sim_id],
        u0=ic[ctx.sim_id],
        p=HybridParams(;
            charge=prob.p.charge,
            mass=prob.p.mass,
            electromagneticfield=prob.p.electromagneticfield,
            magneticmoment=params.magneticmoment[ctx.sim_id],
            terminationcode=TraceParticles.TerminationCode.NotTerminated,
            rng=params.rng[ctx.sim_id],
            initialeomid=params.initialeomid[ctx.sim_id],
            getphase=prob.p.getphase,
        )
    )
end
prob = ODEProblem(hybridgcafo!, zeros(6), tspan, params)

foprob = remake(
    prob,
    f=eoms[3],
    u0=ic[3],
    p=FullOrbitParams(;
        charge=prob.p.charge,
        mass=prob.p.mass,
        electromagneticfield=prob.p.electromagneticfield,
        terminationcode=TraceParticles.TerminationCode.NotTerminated,
    )
)
gcaprob = remake(
    prob,
    f=eoms[4],
    u0=ic[4],
    p=GCAParams(;
        charge=prob.p.charge,
        mass=prob.p.mass,
        electromagneticfield=prob.p.electromagneticfield,
        terminationcode=TraceParticles.TerminationCode.NotTerminated,
        magneticmoment=μ
    )
)

# ODEProblems
# EnsembleProblems
eprob = EnsembleProblem(prob, prob_func=prob_func)

# SOLVE
sol = solve(
    eprob,
    Tsit5();
    #reltol=1e-8,
    trajectories=2,
    callback=cb,
    adaptive=true,
    #dt=tf / 1000
);
fosol = solve(
    foprob,
    Tsit5();
    #reltol=1e-8,
    adaptive=true,
    #dt=tf / 1000
);
gcasol = solve(
    gcaprob,
    Tsit5();
    #reltol=1e-8,
    adaptive=true,
    #dt=tf / 1000
);
times = range(16period, tf, length=1000)
switch1 = [sol.u[1](t, idxs=1:3) for t in times]
switch2 = [sol.u[2](t, idxs=1:3) for t in times]
fo = [fosol(t, idxs=1:3) for t in times]
gca = [gcasol(t, idxs=1:3) for t in times]
rmse1 = sqrt(mean(norm.(switch1 .- fo) .^ 2))
rmse2 = sqrt(mean(norm.(switch2 .- gca) .^ 2))
@testset "Hybrid switching" begin
    @test isapprox(rmse1, 0.0; atol=1e-4)
    @test isapprox(rmse2, 0.0; atol=1e-4)
end

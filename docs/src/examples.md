# Examples

## Basic use of with analytical fields 
This example demonstrates the most basic use of `TraceParticles`. The
electromagnetic fields are analytically prescribed, and the initial conditions
of the particle ensemble are pre-determined. In this case, the only tool that
is needed from `TraceParticles` are the equations of motion for the particles.
The rest is setup for solving an `EnsembleProblem` from `DifferentialEquations`.

```julia
using TraceParticles: lorentzforce!
using DifferentialEquations

"""
We will in this example simulate two charged particles embedded in an
electromagnetic field.
"""
eom = lorentzforce! # The equations of motion of charge particles
                    # in a collisionless electromagntic field is the
                    # Lorentz force.
# Here we define the initial conditions of our two particles. The initial
# conditions are 3 position- and 3 velocity-components, which could be
# represented in a 6-component state vector.
ic = [[0.0, 1.0, 0.0, 1.0, 0.0, 0.0], [0.0, 0.5, 0.0, -0.5, 0.0, 0.0]]
# The parameters of the Lorentz equation is the particle's charge, mass and
# functor a giving the electromagnetic fields at the particle position.
mass = 1
charge = [1, -1]
struct StaticUniformElectromagneticField
    bx::Real
    by::Real
    bz::Real
    ex::Real
    ey::Real
    ez::Real
end
function (field::StaticUniformElectromagneticField)(x,y,z)
    return [field.bx, field.by, field.bz, field.ex, field.ey, field.ez]
end
emfield = StaticUniformElectromagneticField(0,0,1,0,1,0)
tspan = (0,10) # Define the simulation time-span
npart = 2     # Define the number of particles in the simulation
prob = ODEProblem(eom, ic[1], tspan) # Create an initial ODE-problem

# Define a problem function, which will define a new ODE-problem for each
# particle.
function prob_func(prob, i, repeat)
    remake(prob, u0=ic[i], p=(mass, charge[i], emfield))
end
# Define an output function - what to be saved from each particle solution
# The ODESolution type contains a lot of metadata, including alogorithm
# choice, number of algorithm switches, etc. In a large ensemble, this
# overhead will quickly consume a lot of memory.
#
# In this output function we save only the initial and final state of the
# particle, in addition to its timesteps. We set the required "rerun" return
# argument to false.
function output_func(sol,i)
    return ([first(sol), last(sol)], sol.t), false
end
# Define an Ensamble-problem, to simulate the ensemble of particles
ensamble_prob = EnsembleProblem(
    prob # Initial problem
    ;
    prob_func=prob_func,
    output_func=output_func
)

# Run the simulation
sol = DifferentialEquations.solve(
    ensamble_prob, 
    EnsembleSerial()
    ;
    trajectories=npart
)
```

## Gridded fields and statistical sampling
The following example demonstrates a more complex test particle simulation,
using
- interpolated electromagnetic field,
- initial positions drawn from a proposal probability density distribution.
- importance weights corresponding to a target distribution,
- initial velocities drawn from Maxwell-Boltzmann distributions,
- the guiding centre approximation (GCA) as equations of motion,
- callbacks terminating the particle if the GCA is invalid or the particle is
  out of bounds,
- a reduction and output function that stores the initial and final state of the
  solution to an HDF5-file, along with some observables,
- a large particle ensemble.

This script assumes you have the external electromangetic field and temperature of the system, and the  as interpolation objects stored in JLD2 files. `target_itp` and `proposal_itp` are probability density distributions
# used to draw new particle initial conditions.
n 

```julia
using TraceParticles
using DifferentialEquations
using Interpolations
using JLD2

JLD2.@load "electromagneticfield.jld2" em_fields_itp
JLD2.@load "temperature.jld2" temperature_itp
JLD2.@load "target_distribution" target_itp
JLD2.@load "proposal_distribution" proposal_itp

tspan = (0, 1)                 # Temporal bounds of the simulation [s]
xbounds = (0.0, 31.958047e6)   # Spatial bounds [m]
zbounds = (-32e6, 0.00)        # This is a 2.5D simulation

# Particle parameters
charge = -TraceParticles.e      # Charge of particles
mass = TraceParticles.m_e       # Mass of particles
eom = TraceParticles.gca_2Dxz!  # Equation of motion
npart = 10_000_000
rng = MersenneTwister()        # Random number generator
alg = Rosenbrock23()           # Numerical integration algorithm

# Parameters for saving data
expname = "test"
datadir = "data"

# Size of a "particle batch"
batch_size = 1_000_000
nbatches = ceil(Int, npart / batch_size)

# Maximum value of the proposal distribution is needed for the rejection 
# algorithm
maxval = maximum(propdistr.itp.itp.coefs)
samplingdomain = [xbounds, zbounds]

#
# Define Callbacks
#
# Terminate particle when it is outside bounds
out_cb = DiscreteCallback(
    tp.OutOfDomainCondition_2Dxz(xbounds, zbounds),
    tp.outofdomainaffect!
)
# Terminate the particle when EoM is invalid
gca_cb = tp.GCABreakDownCB("lowmem", "2Dxz", tolerance=0.001)
cbset = CallbackSet(out_cb, gca_cb)

# Define the problem function that draws particle initial conditions.
# `DposMBvel` draws position from the proposal distribition and velocity from
# a Maxwell-Boltzmann distribution at the position. It also gives the particle
# a statistical weight according to the ratio between the proposal distribution
# and the target distribution (importance weighing).
prob_func = DposMBvel(
    rng, proposal_itp, target_itp, samplingdomain, maxval, temperature_itp
)
output_func = output_func_max_lightweight
# The output function determines what to with each particle solution.
# `output_func_max_lightweight` keeps only the initial and final state of the
# solution, together with a crude estimate of the largest experienced Larmor
# radius, smallest experienced characteristic length scale of the magnetic
# field, and largest experiences ration between them.

# The reduciton determines what to with each particle batch.
# `SaveBatchAsHDF5` saves each batch in seperate directories in an HDF5-file.
# It also prints some summary statistics for each batch.
reduction = SaveBatchAsHDF5(datadir, expname, batch_size, nbatches)

#
# Define the ODE- and ensemble problem and solve
#
# See the documentation of DifferentialEquations.jl for more keyword arguments
# and solver options
#
prob = ODEProblem(
    eom,
    zeros(4), # Dummy initial condition
    ( # Parameters
        charge=charge,
        mass=mass,
        fields=em_fields_itp,
        magneticmoment=nothing # Dummy magnetic moment
    )
)

ensemble_prob = EnsembleProblem(
    prob;
    prob_func=prob_func,
    output_func=output_func,
    reduction=reduction,
    safetycopy=false
)

sol = solve(
    ensemble_prob,
    alg,
    EnsembleSerial();
    reltol=3e-8,
    abstol=1e0,
    trajectories=npart,
    maxiters=1e6,
    callback=cbset,
    batch_size=batch_size,
)

# Calculate and save `GCAState`s for the initial and final solution for each
# particle. A `GCAState` includes more information than the solution, including
# drifts, Larmor radius, electromagnetic fields etc. See the documentation for
# more information.
save_gcastates(get_filename(reduction), em_fields_itp)
```

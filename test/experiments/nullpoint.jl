# Nullpoint simulation with 
#  Uniform distribution in position 
#  Maxwellian velocity distributions
using tp
using Random
using DifferentialEquations
using Interpolations


# SET EXPEIRMENT PARAMETERS
# -----------------------------------------------------------------------------
# Particle parameters
dtsnap = 1e-2 # in code units (hs)
tspan = (0.0, dtsnap*100)   # Simulation time-span
npart = convert(Int64, 1e4) # Number of particles
xbounds = (14.527344e6, 18.824219e6) # Bounds of of initial position in x
ybounds = (1e6, 1e6)                 # Bounds of of initial position in x
zbounds = (-7.78315e6,-3.8759463e6)  # Bounds of of initial position in x
charge = -tp.e              # Charge of particles
mass = tp.m_e               # Mass of particles
seed = 0                    # Random number generator seed
rng = MersenneTwister(seed) # Random number generator
wp_part=Float64             # Working precision of particles
units = "SI"               # Choice of units for Bifrost snapshot

# Bifrost-input parameters
br_expname="nullpoint"                        # Bifrost experiment name
br_expdir="/Users/eilifo/code/repos/data/8K"  # Bifrost experiment directory
br_snap=1490          # Bifrost experiment snapshot
wp_snap=Float32  # Working precision of Bifrost snapshot
itp_bc = (Flat(), Flat()) # Interpolation boundary conditions

# Parameters for saving data
tp_expname="nullpointexp"                    # Simulation experiment name
tp_expdir="/Users/eilifo/code/repos/data/8K" # Where to save results
 
# Other options
eom = gca_2Dxz! # Equation of motion


# DEFINE PARTICLE INITIAL CONDITIONS  AND ODE-PARAMETERS
# -----------------------------------------------------------------------------
# Background field interpolation object
fields_itp = get_br_emfield_vecof_interpolators(
    br_expname,
    br_snap,
    br_expdir
    ;
    itp_bc = itp_bc
    );

# Initial position of particles
posx = rand(rng, wp_part, xbounds..., npart)
posy = rand(rng, wp_part, ybounds..., npart)
posz = rand(rng, wp_part, zbounds..., npart)

# Initial velocity of particles
# Draw from a Maxwell-distribution with temperature corresponding to the
# initial position of the particles.
br_temp_itp = get_br_var_interpolator(
    br_expname,
    br_snap,
    br_expdir,
    "tg";
    units=units,
    itp_bc = itp_bc
    );
vel = maxwellianvelocitysample(rng, br_temp_itp, mass, posx, posz)
#
# Next, find the guiding centre, velocity magnitude parallel to the magnetic
# field and the magnetic moment of the particle. This is necessary for the
# guiding centre approximation.
#
# First, find EM-fields at initial position
fields_at_pos = [itp.(posx, posz) for itp in fields_itp] # Has type Vector{Vector}
B_at_pos = [[fields_at_pos[i][j] for i = 1:3] for j = 1:npart]
E_at_pos = [[fields_at_pos[i][j] for i = 4:6] for j = 1:npart]
#B_at_pos = [el[1:3] for el in fields_at_pos]
#E_at_pos = [el[4:6] for el in fields_at_pos]

# Then, calculate the 
# and cosine of pitch angle
ic = Vector{Vector{wp_part}}(undef, npart)
magneticmoments = Vector{wp_part}(undef, npart)
for i in eachindex(ic)
    R, vparal, magneticmoment = get_guidingcentre(
    [posx[i], posy[i], posz[i]],
        vel[i],
        B_at_pos[i],
        E_at_pos[i],
        charge,
        mass
        )
    ic[i] = [R; vparal]
    magneticmoments[i] = magneticmoment
end

# DEFINE PROBLEMS
# -----------------------------------------------------------------------------
# Define the problem for the first particle
prob = ODEProblem(
    eom,
    ic[1], 
    tspan, 
    (charge, mass, magneticmoments[1], fields_itp)
    )

# Define the problem function. The function that remakes the ODE-problem for
# every new particle. Here we change the initial conditions and/or the
# parameters.
function prob_func(prob, i, _)
    remake(prob, u0=ic[i], p=(charge, mass, magneticmoments[i], fields_itp))
end

# Define the ensamble problem for the whole lot
@time ensamble_prob = EnsembleProblem(
        prob, 
        prob_func=prob_func
        ;
        safetycopy=false,
    )


# SOLVE PROBLEM
# -----------------------------------------------------------------------------
@time sim = DifferentialEquations.solve(
        ensamble_prob, 
        EnsembleThreads()
        ;
        trajectories=npart,
        progress=true,
        maxiters=1000
    )


# SAVE THE RESULTS
# -----------------------------------------------------------------------------

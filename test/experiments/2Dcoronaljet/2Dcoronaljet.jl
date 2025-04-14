# Create interpolation objects
# To create electron density interpolators demands EoS-tables so we avoid that.
expdir = Base.source_dir()
brxp = BifrostExperiment("testsnap", expdir)
@info "Creating interpolation objects..."
create_bifrost_itps(
    brxp,
    1,
    BSpline(Cubic()),
    Flat(),
    expdir*"/cs",
    auxvariables= ["tg"],
    normalise=[false],
    xzinterpolation=true,
)

# EXPEIRMENT PARAMETERS
#===============================================================================#
# Concerning field, duration and bounds
@info "Loading interpolation objects..."
fields_itp = load_object(expdir*"/cs_BE.jld2")
tg_itp = load_object(expdir*"/cs_tg_normfalse.jld2")
ne_itp = load_object(expdir*"/cs_ne_normfalse.jld2")
ne_norm_itp = load_object(expdir*"/cs_ne_normtrue.jld2")
dtsnap = 1e-2 # in code units (hs)
tspan = (0.0, dtsnap * 100)   # Simulation time-span
#xbounds = (14.527344e6, 18.824219e6) # Bounds of of initial position in x
# x bounds for the whole simulation snapshot
# The whole snapshot except a 40km edge on both sides.
xbounds = (1.5822270393371582e7, 1.6318359375e7)
ybounds = (1e6, 1e6)                 # Bounds of of initial position in y
zbounds = (-6.598462104797363e6, -6.102306842803955e6)

#==============================================================================#
# Particle parameters
charge = -TestParticles.e              # Charge of particles
mass = TestParticles.m_e               # Mass of particles
eom = gca_2Dxz!
npart = Int(1e2)
seed = 5
rng = Xoshiro(seed) # Random number generator
# Parameters for saving data
expname = "testresults"
fname = joinpath(expdir, expname, expname*".h5")
# Where to save the results
datadir = expdir
# Define Callbacks
out_cb = DiscreteCallback(
    OutOfDomainCondition_2Dxz(xbounds, zbounds),
    outofdomainaffect!
)
rel_cb = DiscreteCallback(
    RelativisticConditionGCA_2Dxz(;mass=mass, fraction=0.02),
    relativisticaffect!
)
gca_cb = DiscreteCallback(
    GCABreakDownCondition_2Dxz(0.001),
    gcabreakdownaffect!,
)
# Size of the batch of particles to be saved in a single file
batchsize = Int(1e1)
# Set batchsize to npart if the npart < batchsize
batchsize = npart > batchsize ? batchsize : npart
nbatches = ceil(Int, npart / batchsize)
# Max-value needed for the rejection algorithm
maxval = maximum(ne_itp.itp.itp.itp.coefs)
domain = [xbounds, zbounds] # Domain of sampling.

alg = Rosenbrock23()
"""
Keyword arguments in the creation of the DifferentialEquaions.EnsembleProblem.
"""
ensembleprob_kwargs = Dict(
    # The data footprint and format of output data through an output function
    :output_func => output_func_max_lightweight,
    # How to reduce the simulation data (e.g. how often to save)
    :reduction => SaveBatchAsHDF5(datadir, expname, batchsize, nbatches;
        verbose=false
    ),
    :safetycopy => false,
    # How to choose initial condition of the particle
    :prob_func => DposMBvel_2Dxz(
        rng, ne_itp, ne_itp, domain, maxval, tg_itp
    )
)

"""
Keyword arguments that will be passed to DifferentialEquations.solve() for
solving the ensemble problem. 
"""
solve_kwargs = Dict(
    :reltol => 3e-8,
    :abstol => 1e0,
    :trajectories => npart,
    :maxiters => 10_000_000,
    :callback => CallbackSet(out_cb, rel_cb, gca_cb),
    # The number of particles in each simulation batch.
    # Each batch has its datafile.
    :batch_size => batchsize,
)

println("Parameters set...")

# Run the simulation
#==============================================================================#
odeprob = ODEProblem(
    eom,
    zeros(4),
    tspan,
    (
        charge=charge,
        mass=mass,
        fields=fields_itp,
        magneticmoment=0.0,
    )
)
 
prob = EnsembleProblem(odeprob; ensembleprob_kwargs...)
sol = solve(prob, alg, EnsembleSerial(); solve_kwargs...)

# For testing PposPvel and Relativistic callback
u0 = [
        [1.63e7, 1e6, -6.3e6, 1e6],
        [1.63e7, 1e6, -6.3e6, -1e6],
]
mu = [4.1e-11, 4.1e-11]
ensembleprob_kwargs = Dict(
    :prob_func => PposPvel(u0, mu, [tspan for _ in 1:2])
)
prob2 = EnsembleProblem(odeprob; ensembleprob_kwargs...)
solve_kwargs[:trajectories] = 2
solve_kwargs[:batch_size] = 2
sol2 = solve(prob2, alg, EnsembleSerial(); solve_kwargs...)

# Save GCAStates
save_gcastates(fname, fields_itp)

# Test the results
#==============================================================================#
answersol2 = load_object(expdir*"/testresults/answersol2.jld2")
@test all([u.u for u in sol2.u] .≈ answersol2)

df = h5_getall(fname)
answersol1_df = load_object(expdir*"/testresults/answersol1_df.jld2")
syms = names(answersol1_df)
nsyms = length(syms)
testres = BitVector(undef, nsyms)
for i in 1:nsyms
    if syms[i] ∈ ("retmsg", "retcode")
        testres[i] =  all(df[!, syms[i]] .== answersol1_df[!, syms[i]])
    else
        testres[i] =  all(df[!, syms[i]] .≈ answersol1_df[!, syms[i]])
    end
end
@test all(testres)

e0, ef = h5_getenergies(fname)
answersol2_e0, answersol2_ef = load_object(
    expdir*"/testresults/answersol1_energies.jld2"
)
@test all(e0 .≈ answersol2_e0)
@test all(ef .≈ answersol2_ef)

rm(expdir*"/cs_BE.jld2")
rm(expdir*"/cs_tg_normfalse.jld2")
rm(fname)
rm(joinpath(expdir, expname, expname*"_gcastates0.h5"))
rm(joinpath(expdir, expname, expname*"_gcastatesf.h5"))
# Create interpolation objects
# To create electron density interpolators demands EoS-tables so we avoid that.
expdir = Base.source_dir()
brxp = BifrostExperiment("testsnap", expdir)
create_bifrost_itps(
    brxp,
    1,
    BSpline(Cubic()),
    Flat(),
    expdir*"/cs",
    "si";
    variables= ["BE", "tg"],
    normalise=[false, false],
    destagger=[true, false],
    xzinterpolation=true,
)

# EXPEIRMENT PARAMETERS
#===============================================================================#
# Concerning field, duration and bounds
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
charge = -TraceParticles.e              # Charge of particles
mass = TraceParticles.m_e               # Mass of particles
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
# Default rtol for isapprox is sqrt(eps()) which is a bit more than 1e-8
minrtol=1e-8

# These are the columns of the results that should be exactly equal to the
# answers, independent of machine.
exactsyms = (
    "acceptanceratio",
    "charge",
    "mass",
    "retmsg",
    "retcode",
    "x0",
    "y0",
    "z0",
    "t0",
    "weight",
)
# These are the columns that have a slight error in the results across machines,
# but for most particles the relative error is lower than 1e-8
minrtolsyms = (
    "vparal0",
    "magneticmoment",
    "xf",
    "yf",
    "zf",
    "vparalf",
    "tf",
    "maxrl",
)
# These are the columns that have larger errors across machines for certain
# particles. For these particles, the tolerance is adjusted.
rtols_errorprone = Dict(
    "xf" => 1e-6,
    "zf" => 1e-6,
    "tf" => 1e-5,
    "nt" => 1e-2,
    "vparalf" => 1e-2,
    "maxrl" => 1e-4,
    "maxscalesratio" => 1e-3,
    "minlb" => 1e-3,
    "meanrl" => 1e-6,
    "meanlb" => 1e-1,
    "meanscalesratio" => 1e-4,
)

# Retrieve the solutions for which to test the results against.
answersol2 = load_object(expdir*"/testresults/answersol2.jld2")
answersol1_df = load_object(expdir*"/testresults/answersol1_df.jld2")
answersol1_e0, answersol1_ef = load_object(
    expdir*"/testresults/answersol1_energies.jld2"
)

# Get the results
e0, ef = h5_getenergies(fname)
df = h5_getall(fname)

# The particles that have larger errors
errorprone = [75, 81, 1, 14, 49, 70, 46, 93, 84, 44, 53, 80, 27, 9, 13, 41, 30, 39]
# Create a filter for masking out the errorprone particles
totest = BitVector(undef, 100)
totest .= [i ∉ errorprone for i in 1:100]

# Get the column names
syms = names(answersol1_df)
nsyms = length(syms)

# Compare the results of the particles that are not errorprone.
testres = BitVector(undef, nsyms)
testres .= 1
for i in 1:nsyms
    sym = syms[i]
    if sym ∈ exactsyms
        testres[i] =  all(df[totest, sym] .== answersol1_df[totest, sym])
    elseif sym ∈ minrtolsyms
        testres[i] =  all(isapprox.(
            df[totest, sym],
            answersol1_df[totest, sym],
            rtol=minrtol
            )
        )
    end
end

# Compare the results of the particles that are errorprone
testres_errorprone = BitVector(undef, nsyms)
testres_errorprone .= 1
for i in 1:nsyms
    sym = syms[i]
    if sym ∈ keys(rtols_errorprone)
        testres_errorprone[i] =  all(isapprox.(
            df[errorprone, sym],
            answersol1_df[errorprone, sym],
            rtol=rtols_errorprone[sym]
            )
        )
    elseif sym ∈ exactsyms
        testres_errorprone[i] =  all(df[errorprone, sym] .== answersol1_df[errorprone, sym])
    elseif sym ∈ minrtolsyms
        testres_errorprone[i] =  all(isapprox.(
            df[errorprone, sym],
            answersol1_df[errorprone, sym],
            rtol=minrtol
            )
        )

    end
end

# Remove gnerated files before testing.
rm(expdir*"/cs_BE.jld2")
rm(expdir*"/cs_tg_normfalse.jld2")
rm(fname)
rm(joinpath(expdir, expname, expname*"_gcastates0.h5"))
rm(joinpath(expdir, expname, expname*"_gcastatesf.h5"))

# Test the results form the relativistic particles in sol2
@test all([u.u for u in sol2.u] .≈ answersol2)
# Test the columns of the not so errorprone particles
@test all(testres)
# Test the initial energy of all particles
@test all(isapprox.(e0, answersol1_e0; rtol=minrtol))
# The final energy of the not so errorprone particles
@test all(isapprox.(ef[totest], answersol1_ef[totest]; rtol=minrtol))
# Test the columns of the errorprone particles
@test all(testres_errorprone)
# Test the final energy of the errorprone particles
@test all(isapprox.(ef[errorprone], answersol1_ef[errorprone]; rtol=6e-6))

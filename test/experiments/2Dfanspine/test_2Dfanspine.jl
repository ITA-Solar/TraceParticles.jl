using BifrostTools
using DataFrames
using Interpolations
using JLD2
using Logging
using OrdinaryDiffEq
using Random
using SciMLBase
using StaticArrays
using Test
using TraceParticles

expdir = Base.source_dir()

# Create interpolation objects
# To create electron density interpolators demands EoS-tables so we avoid that.
brxp = BifrostExperiment(
    "testsnap",
    expdir
)
create_bifrost_itps(
    brxp,
    1,
    BSpline(Cubic()),
    Flat(),
    expdir * "/cs",
    "si",
    ["BE", "tg", "r"];
    normalise=[false, false, false],
    destagger=[true, false, false],
    xzinterpolation=true,
)

# Experiment parameters
include(joinpath(expdir, "test_2Dfanspine.params.jl"))

# Construct the problem
tpprob = TraceParticlesProblem(params; logger=current_logger())
# Run the simulation
solve(tpprob, params)

#==============================================================================#
(; ensemble_prob, callbackset) = tpprob
(; electromagneticfield) = ensemble_prob.prob.p

# For testing PposPvel and Relativistic callback
u0 = [
    [1.63e7, 1e6, -6.3e6, 1e6],
    [1.63e7, 1e6, -6.3e6, -1e6],
]
mu = [4.1e-11, 4.1e-11]
tspan2 = [(params.t0start, params.tf) for _ in 1:2]
params2 = [
    GCAParams(
        charge=params.charge,
        mass=params.mass,
        magneticmoment=mu[i],
        electromagneticfield=electromagneticfield
    ) for i in 1:2
]
ensembleprob_kwargs = Dict(
    :prob_func => PredefinedICs(u0, tspan2, params2)
)
prob2 = EnsembleProblem(ensemble_prob.prob; ensembleprob_kwargs...)
solve_kwargs = Dict(
    :trajectories => 2,
    :batch_size => 2,
    :abstol => params.abstol,
    :reltol => params.reltol,
    :callback => callbackset,
    :maxiters => params.maxiters,
)
sol2 = solve(prob2, params.alg, EnsembleSerial(); solve_kwargs...)

# Save GCAStates
resultsdir = joinpath(expdir, params.expname)
fname = joinpath(resultsdir, params.expname * ".h5")
save_gcastates(fname, electromagneticfield)

# Test the results
#==============================================================================#
# These are the columns of the results that should be exactly equal to the
# answers, independent of machine.
exactsyms = (
    "nrejections",
    "charge",
    "mass",
    "terminationcode",
    "retcode",
    "x0",
    "y0",
    "z0",
    "t0",
    "weight",
)
# These are the columns that have a slight error in the results across machines,
# but for most particles the relative error is lower than 1e-8
# Default rtol for isapprox is sqrt(eps()) which is a bit more than 1e-8
#minrtol = 1e-20 # for testing
minrtol = 1e-8
rtols = Dict(
    "vparal0" => minrtol,
    "magneticmoment" => minrtol,
    "xf" => minrtol,
    "yf" => minrtol,
    "zf" => minrtol,
    "vparalf" => minrtol,
    "tf" => minrtol,
    "maxrl" => minrtol,
)

# These are the columns that have larger errors across machines for certain
# particles. For these particles, the tolerance is adjusted.
rtols_errorprone = Dict(
    "xf" => 1e-6,
    "zf" => 1e-6,
    "vparalf" => 2e-4,
    "maxrl" => 5e-3,
)

# Retrieve the solutions for which to test the results against.
answersol2 = load_object(expdir * "/answers/answersol2.jld2")
answersol1_df = load_object(expdir * "/answers/answersol1_df.jld2")
answersol1_e0, answersol1_ef = load_object(
    expdir * "/answers/answersol1_energies.jld2"
)
# Get the results
e0, ef = h5_getenergies_gca(fname)
# Check that the results from `save_gcastates` is the same as the results from
# `save_energies`
@test (e0, ef) == h5_getenergies(fname)
df = h5_getall(fname)

# The particles that have larger errors
errorprone = [
    25, 40, 50, 99
]
# Create a filter for masking out the errorprone particles
totest = BitVector(undef, 100)
totest .= [i ∉ errorprone for i in 1:100]

# Get the column names
syms = [exactsyms..., keys(rtols)...]
nsyms = length(syms)

# Compare the results of the particles that are not errorprone.
testres = BitVector(undef, nsyms)
testres .= 1
for i in 1:nsyms
    sym = syms[i]
    if sym ∈ exactsyms
        testres[i] = all(df[totest, sym] .== answersol1_df[totest, sym])
    elseif sym ∈ keys(rtols)
        testres[i] = all(isapprox.(
            df[totest, sym],
            answersol1_df[totest, sym],
            rtol=rtols[sym]
        )
        )
    end
end

# Compare the results of the particles that are errorprone
syms_errorprone = unique([exactsyms..., keys(rtols)..., keys(rtols_errorprone)...])
nsyms_errorprone = length(syms)
testres_errorprone = BitVector(undef, nsyms_errorprone)
testres_errorprone .= 1
for i in 1:nsyms_errorprone
    sym = syms_errorprone[i]
    if sym ∈ keys(rtols_errorprone)
        testres_errorprone[i] = all(isapprox.(
            df[errorprone, sym],
            answersol1_df[errorprone, sym],
            rtol=rtols_errorprone[sym]
        )
        )
    elseif sym ∈ exactsyms
        testres_errorprone[i] = all(df[errorprone, sym] .== answersol1_df[errorprone, sym])
    elseif sym ∈ keys(rtols)
        testres_errorprone[i] = all(isapprox.(
            df[errorprone, sym],
            answersol1_df[errorprone, sym],
            rtol=rtols[sym]
        )
        )

    end
end

# Remove gnerated files before testing.
rm(expdir * "/cs_t1_BE.jld2")
rm(expdir * "/cs_t1_tg_normfalse.jld2")
rm(expdir * "/cs_t1_r_normfalse.jld2")
rm(resultsdir, recursive=true)


# Test the results form the relativistic particles in sol2
@test all([u.u for u in sol2.u] .≈ answersol2)
# Test the columns of the not so errorprone particles
str = ""
for i in eachindex(testres)
    global str
    str *= "$(testres[i]): $(syms[i])\n"
end
#@warn str
str2 = ""
for i in eachindex(testres)
    global str2
    reldiff = abs.((df[totest, syms[i]] ./ answersol1_df[totest, syms[i]]) .- 1)
    mask = findall(x -> x >= minrtol, reldiff)
    str2 *= "$(syms[i]):\n"
    for j in mask
        str2 *= "  $j : $(reldiff[j])\n"
    end
end
#@warn str2
@test all(testres)
# Test the initial energy of all particles
@test all(isapprox.(e0, answersol1_e0; rtol=minrtol))
# The final energy of the not so errorprone particles
@test all(isapprox.(ef[totest], answersol1_ef[totest]; rtol=minrtol))
# Test the columns of the errorprone particles
@test all(testres_errorprone)
# Test the final energy of the errorprone particles
@test all(isapprox.(ef[errorprone], answersol1_ef[errorprone]; rtol=6e-6))

#-------------------------------------------------------------------------------
#           test_rerun.jl
#-------------------------------------------------------------------------------
# Unit tests for re-integrating a selection of particles from a finished
# experiment, keeping their full trajectories instead of the endpoint summary
# the ensemble stored.
#
# A rerun reads an experiment back from disk, so these tests first run a small
# experiment of their own and rerun particles out of it. The field is uniform
# and written to JLD2.
#-------------------------------------------------------------------------------

using HDF5
using Interpolations
using JLD2
using Logging
using OrdinaryDiffEq
using SciMLBase
using Test

using TraceParticles:
    electron_charge,
    electron_mass,
    GCAParams,
    OutOfBounds,
    PredefinedICs,
    Rerun,
    SaveBatchAsHDF5,
    TerminationCode,
    TraceParticlesParameters,
    TraceParticlesProblem,
    guidingcentreapproximation!,
    rerun,
    solve,
    h5_getdataset

if !isdefined(Main, :verbose)
    verbose = 4
end

dummy_field(x, y, z, t) = [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]


"""
    seed_experiment(dir; npart=4, expname="seed", kwargs...)

Run a small guiding-centre experiment in `dir` and return `(params, datafile)`.
Each particle gets its own parallel velocity `i*1e5` and magnetic moment
`i*1e-20`, so a rerun of a selection can be checked to have picked up the right
ones.

The initial conditions are handed over as a `PredefinedICs` problem function
rather than through a file, which keeps the seeded values in plain sight.
`paramsbackup` is off because the backup copies the running script, which here
is the test runner, into the experiment directory.
"""
function seed_experiment(dir; npart=4, expname="seed", kwargs...)
    tf = 1.0
    field(x, y, z, t) = [0, 0, 0], [0, 0, 1e-4]
    u0 = [[0.0, 0.0, 0.0, i*1e5] for i in 1:npart]
    tspans = [(0.0, tf) for _ in 1:npart]
    odeparams = [
        GCAParams(
            charge=electron_charge,
            mass=electron_mass,
            electromagneticfield=field,
            magneticmoment=i*1e-20,
        )
        for i in 1:npart
    ]
    params = TraceParticlesParameters(;
        eom = guidingcentreapproximation!,
        npart = npart,
        tf = tf,
        mass = electron_mass,
        charge = electron_charge,
        alg = Tsit5(),
        maxiters = 100_000,
        reltol = 1e-8,
        abstol = 1e-8,
        expname = expname,
        datadir = dir,
        electromagneticfield = field,
        prob_func = PredefinedICs(u0, tspans, odeparams),
        verbose = false,
        paramsbackup = false,
        log2file=false,
        kwargs...
    )
    solve(TraceParticlesProblem(params; logger=current_logger()), params)
    # The file `SaveBatchAsHDF5` writes, which is what a rerun reads back.
    return params, joinpath(dir, expname, expname*".h5")
end

@testset verbose = verbose ≥ 4 "rerun.jl" begin
    mktempdir() do dir
        npart = 4
        params, datafile = seed_experiment(dir; npart=npart)
        @test isfile(datafile)

        idxs = [4, 2]
        zmax = 3e5
        sol = rerun(
            idxs,
            datafile,
            params;
            alg=Vern9(),
            callback_specs=(OutOfBounds(z=(-zmax, zmax)),)
        )

        # Whole trajectories, not the two-point summary the experiment
        # stored, and one per selected particle.
        @test sol isa SciMLBase.EnsembleSolution
        @info size(sol)
        @test length(sol.u) == length(idxs)
        @test all(s -> s isa SciMLBase.ODESolution, sol.u)
        @test all(s -> length(s.t) > 2, sol.u)
        @test all(SciMLBase.successful_retcode, sol.u)
        @test sol.u[1].alg isa Vern9
        @test sol.u[2].u[end][3] ≈ 2e5*params.tf rtol=1e-6
        @test sol.u[1].prob.p.terminationcode ==
            TerminationCode.OutOfBounds
        @test sol.u[1].u[end][3] ≥ zmax
        @test zmax/4e5 ≤ sol.u[1].t[end] ≤ params.tf

        # Each particle starts from the state and magnetic moment the
        # experiment recorded for it, in the order asked for.
        x0 = h5_getdataset(datafile, "x0")
        y0 = h5_getdataset(datafile, "y0")
        z0 = h5_getdataset(datafile, "z0")
        vparal0 = h5_getdataset(datafile, "vparal0")
        u0 = [[x0[i], y0[i], z0[i], vparal0[i]] for i in idxs]
        mu0 = h5_getdataset(datafile, "magneticmoment")
        for (i, idx) in enumerate(idxs)
            @test sol.u[i].u[1] ≈ u0[i]
            @test sol.u[i].prob.p.magneticmoment == mu0[idx]
            # Seeded as v∥ = idx*1e5 and μ = idx*1e-20.
            @test sol.u[i].u[1][4] ≈ idx*1e5
            @test sol.u[i].prob.p.magneticmoment ≈ idx*1e-20
        end

        # And over the time span the experiment gave it.
        t0 = h5_getdataset(datafile, :t0)
        tf = h5_getdataset(datafile, :tf)
        for (i, idx) in enumerate(idxs)
            @test sol.u[i].t[1] == t0[idx]
            @test sol.u[i].t[end] ≈ tf[idx]
        end

        # The experiment's own parameters are left alone, so it can be
        # rerun again with a different selection.
        @test params.expname == "seed"
        @test params.npart == npart
        @test params.prob_func isa PredefinedICs

        # A misspelled parameter is rejected rather than silently ignored.
        @test_throws ArgumentError rerun(
            idxs, datafile, params; algr=Vern9()
        )
    end
end

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

import TraceParticles as TP
using TraceParticles:
    ElectromagneticFieldInterpolator,
    GCAParams,
    OutOfBounds,
    PredefinedICs,
    Rerun,
    SaveBatchAsHDF5,
    TerminationCode,
    TraceParticlesParameters,
    TraceParticlesProblem,
    guidingcentreapproximation!,
    initialstate_idxs_gca,
    rerun,
    solve,
    tspans_idxs

if !isdefined(Main, :verbose)
    verbose = 4
end

"""
A serialisable uniform field component. JLD2 cannot round-trip an anonymous
function, so the emfield file holds these instead.
"""
struct UniformField{T<:Real}
    value::T
end
(self::UniformField)(x, y, z, t) = self.value
(self::UniformField)(x, y, z) = self.value

"""
Write the JLD2 file holding the 6-element `[bx, by, bz, ex, ey, ez]` vector of
callables that an experiment loads its field from, and return its path.
"""
function write_rerun_emfield(dir; b=(0.0, 0.0, 1e-4), e=(0.0, 0.0, 0.0))
    path = joinpath(dir, "emfield.jld2")
    save_object(path, [UniformField.(b)..., UniformField.(e)...])
    return path
end

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
    field = ElectromagneticFieldInterpolator(
        [UniformField(0.0), UniformField(0.0), UniformField(1e-4),
         UniformField(0.0), UniformField(0.0), UniformField(0.0)]
    )
    u0 = [[0.0, 0.0, 0.0, i*1e5] for i in 1:npart]
    tspans = [(0.0, tf) for _ in 1:npart]
    odeparams = [
        GCAParams(
            charge=-TP.e,
            mass=TP.m_e,
            electromagneticfield=field,
            magneticmoment=i*1e-20,
        )
        for i in 1:npart
    ]
    params = TraceParticlesParameters(;
        eom = guidingcentreapproximation!,
        npart = npart,
        tf = tf,
        mass = TP.m_e,
        charge = -TP.e,
        alg = Tsit5(),
        maxiters = 100_000,
        reltol = 1e-8,
        abstol = 1e-8,
        expname = expname,
        datadir = dir,
        emfield_file = write_rerun_emfield(dir),
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

    @testset verbose = verbose ≥ 5 "rerun" begin
        mktempdir() do dir
            npart = 4
            params, datafile = seed_experiment(dir; npart=npart)
            @test isfile(datafile)

            idxs = [2, 4]
            sol = rerun(idxs, datafile, params; paramsbackup=false)

            # Whole trajectories, not the two-point summary the experiment
            # stored, and one per selected particle.
            @test sol isa SciMLBase.EnsembleSolution
            @test length(sol.u) == length(idxs)
            @test all(s -> s isa SciMLBase.ODESolution, sol.u)
            @test all(s -> length(s.t) > 2, sol.u)
            @test all(SciMLBase.successful_retcode, sol.u)

            # Each particle starts from the state and magnetic moment the
            # experiment recorded for it, in the order asked for.
            u0, mu0 = initialstate_idxs_gca(datafile, idxs)
            for (i, idx) in enumerate(idxs)
                @test sol.u[i].u[1] ≈ u0[i]
                @test sol.u[i].prob.p.magneticmoment == mu0[i]
                # Seeded as v∥ = idx*1e5 and μ = idx*1e-20.
                @test sol.u[i].u[1][4] ≈ idx*1e5
                @test sol.u[i].prob.p.magneticmoment ≈ idx*1e-20
            end

            # And over the time span the experiment gave it.
            for (i, tspan) in enumerate(tspans_idxs(datafile, idxs))
                @test sol.u[i].t[1] == tspan[1]
                @test sol.u[i].t[end] ≈ tspan[2]
            end

            # Uniform B: the guiding centre slides along the field line at a
            # constant parallel velocity and drifts nowhere.
            for (i, idx) in enumerate(idxs)
                @test sol.u[i].u[end][1] ≈ 0.0 atol=1e-6
                @test sol.u[i].u[end][2] ≈ 0.0 atol=1e-6
                @test sol.u[i].u[end][3] ≈ idx*1e5*params.tf rtol=1e-6
                @test sol.u[i].u[end][4] ≈ idx*1e5 rtol=1e-6
            end

            # The selection is taken in the order given, and any
            # `AbstractVector` of indices will do.
            reversed = rerun(
                reverse(idxs), datafile, params; paramsbackup=false
            )
            @test reversed.u[1].u[1] ≈ sol.u[2].u[1]
            @test reversed.u[2].u[1] ≈ sol.u[1].u[1]
            @test length(
                rerun(1:npart, datafile, params; paramsbackup=false).u
            ) == npart

            # The experiment's own parameters are left alone, so it can be
            # rerun again with a different selection.
            @test params.expname == "seed"
            @test params.npart == npart
            @test params.prob_func isa PredefinedICs

            # The particles are read from the file given, not from a name
            # derived from the parameters: a copy of the data elsewhere is
            # rerun just the same.
            elsewhere = joinpath(dir, "copy.h5")
            cp(datafile, elsewhere)
            fromcopy = rerun(idxs, elsewhere, params; paramsbackup=false)
            @test fromcopy.u[1].u[1] ≈ sol.u[1].u[1]
            @test fromcopy.u[1].u[end] ≈ sol.u[1].u[end]

            # The rerun keeps its trajectories in memory rather than writing a
            # data file of its own.
            @test !isfile(joinpath(dir, "seed", "seed_rerun.h5"))
            @test length(filter(
                f -> endswith(f, ".h5"), readdir(dir; join=true)
            )) == 1 # only the copy made above
        end
    end

    @testset verbose = verbose ≥ 5 "rerun parameter changes" begin
        mktempdir() do dir
            params, datafile = seed_experiment(dir)
            idxs = [1, 3]

            # Any parameter can be varied for the rerun. A different
            # integrator solves the same trajectory...
            varied = rerun(
                idxs, datafile, params; alg=Vern9(), paramsbackup=false
            )
            @test varied.u[1].alg isa Vern9
            @test varied.u[1].u[end][3] ≈ 1e5*params.tf rtol=1e-6

            # ...and a termination criterion the experiment did not have stops
            # the particle short, with the code the criterion sets. This is a
            # `DiscreteCallback`, so it fires at the end of the step that left
            # the box, somewhere past the boundary.
            zmax = 1e4
            bounded = rerun(
                idxs, datafile, params;
                callback_specs=(OutOfBounds(z=(-zmax, zmax)),),
                paramsbackup=false,
            )
            @test bounded.u[1].prob.p.terminationcode ==
                TerminationCode.OutOfBounds
            @test bounded.u[1].u[end][3] ≥ zmax
            @test zmax/1e5 ≤ bounded.u[1].t[end] ≤ params.tf

            # A `maxiters` too small to reach the end time stops the particle
            # with the matching return code.
            capped = with_logger(NullLogger()) do
                rerun(idxs, datafile, params; maxiters=3, paramsbackup=false)
            end
            @test capped.u[1].retcode == SciMLBase.ReturnCode.MaxIters
            @test capped.u[1].t[end] < params.tf

            # A misspelled parameter is rejected rather than silently ignored.
            @test_throws ArgumentError rerun(
                idxs, datafile, params; algorithm=Vern9()
            )

            # Solver keywords are not parameters, and are rejected as well.
            @test_throws ArgumentError rerun(
                idxs, datafile, params; saveat=0.1
            )
        end
    end

    @testset verbose = verbose ≥ 5 "rerun experiment directory" begin
        mktempdir() do dir
            params, datafile = seed_experiment(dir)

            # The rerun is an experiment of its own: it gets a directory named
            # after the original, and does not write into the original's.
            before = readdir(dir)
            rerun([1], datafile, params; paramsbackup=false)
            new = setdiff(readdir(dir), before)
            @test length(new) == 1
            @test startswith(only(new), "seed-rerun-")
            @test isdir(joinpath(dir, only(new)))

            # `paramsbackup` copies the running script into that directory. It
            # is on by default, whatever the experiment was run with, and here
            # the running script is the test runner.
            @test !params.paramsbackup
            rerun([1], datafile, params)
            backedup = setdiff(readdir(dir), before, new)
            expdir = joinpath(dir, only(backedup))
            @test isfile(joinpath(expdir, only(backedup)*".jl")) ==
                isfile(PROGRAM_FILE)
        end
    end

end # testset rerun.jl

"""
    solve(prob::TraceParticlesProblem, params::TraceParticlesParameters)
Solve a `TraceParticlesProblem` given a `TraceParticlesParameters` instance.
"""
function SciMLBase.solve(
        prob::TraceParticlesProblem,
        params::TraceParticlesParameters;
        logger = nothing
)
    logger = isnothing(logger) ? prob.logger : logger
    with_logger(logger) do
        (; datadir, expname, abstol, reltol, maxiters, alg, seed) = params
        create_paramsbackup(datadir, expname)
        return solve(
            prob,
            alg,
            abstol,
            reltol,
            maxiters,
            seed,
        )
    end
end

"""
    solve(prob::TraceParticlesProblem; kwargs...)
Solve a `TraceParticlesProblem` and compute the initial and final energy of the
particles using `save_energy`.
"""
function SciMLBase.solve(
    prob::TraceParticlesProblem,
    alg::SciMLBase.AbstractODEAlgorithm,
    abstol::Real,
    reltol::Real,
    maxiters::Number,
    seed::Number,
    kwargs...
)
    # Find the number of processes working in parellel.
    np = nprocs()
    if np > 1
        # If there are more than 1 worker, solve the problem using distributed
        # parallelisation.
        @warn "Running on $np processes. Make sure the the problem function
            and necessary packages are defined on all workers using
            `@everywhere`"
        paral_msg = "on $np processes"
        ensemble_algorithm = EnsembleDistributed()
    elseif Threads.nthreads() > 1
        # If there are multiple threads, solve the problem multi-threaded.
        paral_msg = "on $(Threads.nthreads()) threads"
        ensemble_algorithm = EnsembleThreads()
    else
        # Solve the problem serially.
        paral_msg =  "in serial"
        ensemble_algorithm = EnsembleSerial()
    end
    #==========================================================================#
    # START OF SIMULATION
    # Unpack struct fields
    (; ensemble_prob, callbackset, trajectories, batch_size) = prob

    @info """
------------------------------------------------------------------------
Running program: $PROGRAM_FILE $paral_msg
Start time: $(string(now()))
Host name : $(gethostname())
Running an ensemble of $trajectories particles
TraceParticles.jl package version: $(pkgversion(TraceParticles))
------------------------------------------------------------------------------"""


    time_taken = @timed sim = SciMLBase.solve(
        ensemble_prob,
        alg,
        ensemble_algorithm
        ;
        trajectories=trajectories,
        reltol=reltol,
        abstol=abstol,
        maxiters=maxiters,
        batch_size=batch_size,
        callback=callbackset,
        kwargs...
    );
    @info """Ensemble solved:
    Elapsed time: $(time_taken.time) seconds
    Allocations: $(time_taken.bytes) bytes"""

    if ensemble_prob.reduction == SciMLBase.DEFAULT_REDUCTION
        @warn "No reduction function specified.
        Saving the full ensemble simulation as 'out.jld2'."
        JLD2.@save "out.jld2" sim
    elseif typeof(ensemble_prob.reduction) <: SaveBatchAsHDF5
        try
            filename = get_filename(ensemble_prob.reduction)
            time_taken = @timed save_energy(
                filename,
                ensemble_prob.prob.p.electromagneticfield;
                eom=ensemble_prob.prob.f.f
            )
            @info """Energies computed:
            Elapsed time: $(time_taken.time) seconds
            Allocations: $(time_taken.bytes) bytes"""
        catch e
            @error "Failed to save energies" exception=e
        end
    end
    return sim
end



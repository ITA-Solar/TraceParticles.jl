"""
    solve(prob::TraceParticlesProblem, params::TraceParticlesParameters)
Solve a `TraceParticlesProblem` given a `TraceParticlesParameters` instance.
"""
function SciMLBase.solve(
        prob::TraceParticlesProblem,
        params::TraceParticlesParameters;
        logger = nothing,
        seed = params.seed
)
    logger = isnothing(logger) ? prob.logger : logger
    with_logger(logger) do
        (; datadir, expname, abstol, reltol, maxiters, alg) = params
        if params.paramsbackup
            create_paramsbackup(datadir, expname)
        end
        sim = solve(
            prob,
            alg,
            abstol,
            reltol,
            maxiters,
            seed
        )

        if params.compute_initial_and_final_energy
            if get_reductiontype(prob) <: SaveBatchAsHDF5
                try
                    filename = get_filename(prob)
                    time_taken = @timed save_energy(
                        filename,
                        get_electromagneticfield(prob);
                        eom=get_eom(prob)
                    )
                    logtimed(time_taken, "Energies computed")
                catch e
                    @error "Failed to save energies" exception=e
                end
            else
                @info "Did not compute energies: Uknown reduction"
            end
        end
        return sim
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
        paral_msg = "$np processes"
        ensemble_algorithm = EnsembleDistributed()
    elseif Threads.nthreads() > 1
        # If there are multiple threads, solve the problem multi-threaded.
        paral_msg = "$(Threads.nthreads()) threads"
        ensemble_algorithm = EnsembleThreads()
    else
        # Solve the problem serially.
        paral_msg =  "Serial"
        ensemble_algorithm = EnsembleSerial()
    end
    #==========================================================================#
    # START OF SIMULATION
    # Unpack struct fields
    (; ensemble_prob, callbackset, trajectories, batch_size) = prob

    @info """
------------------------------------------------------------------------
TraceParticles.jl package version: $(pkgversion(TraceParticles))
Program: $PROGRAM_FILE
Experiment name: $(experimentname(prob))
Parallelisation: $paral_msg
Host name : $(gethostname())
Start time: $(string(now()))

Running an ensemble of $(@sprintf "%.1E" trajectories) particles
Random seed: $seed
Equations of motion: $(get_eom(prob))
Max iterations per particle: $maxiters
Problem function: $(typeof(get_prob_func(prob)).name.name)
Output function: $(typeof(get_output_func(prob)).name.name)
Reduction: $(typeof(get_reduction(prob)).name.name)
Number of batches: $(get_nofbatches(prob))
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
        seed=seed,
        kwargs...
    );
    logtimed(time_taken, "Ensemble solved")

    return sim
end

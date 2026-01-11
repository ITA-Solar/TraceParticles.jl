"""
    solve(prob::TraceParticlesProblem, params::TraceParticlesParameters)
Solve a `TraceParticlesProblem` given a `TraceParticlesParameters` instance.
"""
function DifferentialEquations.solve(
        prob::TraceParticlesProblem,
        params::TraceParticlesParameters
)
    (; datadir, expname, abstol, reltol, maxiters, alg) = params
    create_expdir(datadir, expname)
    create_paramsbackup(datadir, expname)
    return solve(prob;
        solver_algorithm=alg,
        abstol=abstol,
        reltol=reltol,
        maxiters=maxiters
    )
end

"""
    solve(prob::TraceParticlesProblem; kwargs...)
Solve a `TraceParticlesProblem` and compute the initial and final energy of the
particles using `save_energy`.

The default solver algorithm is `Tsit5` unless the keyword argument
`solver_algorithm` is set to something different. Other default keyword
arguments are:
- `abstol=1e-6`
- `abstol=1e-3`
- `maxiters=1e5`
"""
function DifferentialEquations.solve(
    prob::TraceParticlesProblem;
    solver_algorithm::SciMLBase.AbstractODEAlgorithm = Tsit5(),
    abstol::Real = 1e-6,
    reltol::Real = 1e-3,
    maxiters::Number = 1e5,
    kwargs...
)
    # Find the number of processes working in parellel.
    np = nprocs()
    println(".................................................................",
        "...............")
    if np > 1
        # If there are more than 1 worker, solve the problem using distributed
        # parallelisation.
        @warn "Running on $np processes. Make sure the the problem function and
            necessary packages are defined on all workers using `@everywhere`"
        ensemble_algorithm = EnsembleDistributed()
    elseif Threads.nthreads() > 1
        # If there are multiple threads, solve the problem multi-threaded.
        println("Running on $(Threads.nthreads()) threads.")
        ensemble_algorithm = EnsembleThreads()
    else
        # Solve the problem serially.
        println("Running in serial.")
        ensemble_algorithm = EnsembleSerial()
    end
    #==========================================================================#
    # START OF SIMULATION
    # Unpack struct fields
    (; ensemble_prob, callbackset, trajectories, batch_size) = prob
    
    println("Running program: $PROGRAM_FILE")
    println("Start time: $(string(now()))")
    println("Host name : $(gethostname())")
    println("Running ensemble of $trajectories particles...")
    println(".................................................................",
        "...............")
    
    @time sim = DifferentialEquations.solve(
        ensemble_prob,
        solver_algorithm,
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

    if ensemble_prob.reduction == SciMLBase.DEFAULT_REDUCTION
        @warn "No reduction function specified.
        Saving the full ensemble simulation as 'out.jld2'."
        JLD2.@save "out.jld2" sim
    elseif typeof(ensemble_prob.reduction) <: SaveBatchAsHDF5
        try
            println("\nPost-processing: Computing energies...")
            filename = get_filename(ensemble_prob.reduction)
            @time save_energy(
                filename,
                ensemble_prob.prob.p.electromagneticfield,
                ensemble_prob.prob.f.f
            )
        catch e
            @error "Failed to save energies:\n $e"
        end
    end
end



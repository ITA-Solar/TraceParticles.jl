function tprun(
    solve_kwargs,
    ensembleprob_kwargs,
    eom,
    npart,
    tspan,
    charge,
    mass,
    fields_itp,
    prob_func_in,
)
    # check_requirements(params_struct)
    # create_workingdir(parameterfile)
    # ensemble_algorithm = get_ensemble_algorithm()
    # ensembleprob_kwargs_internal = add_prob_func_counter(ensemble_algorithm)
    # ensemble_prob = create_ensemble_prob(params_struct)
    # run()

    np = nprocs()
    println(".................................................................",
        "...............")
    if np > 1
        println("Running on $np processes.")
        ensemble_algorithm = EnsembleDistributed()
    elseif Threads.nthreads() > 1
        println("Running on $(Threads.nthreads()) threads.")
        ensemble_algorithm = EnsembleThreads()
    else
        println("Running in serial.")
        ensemble_algorithm = EnsembleSerial()
    end
    solve_args = try
        (alg, ensemble_algorithm)
    catch
        (ensemble_algorithm,)
    end


    function print_particlecounter(count, total)
        @printf "Particle %i (%.2f %%)\r" count count / total * 100
    end

    lk = ReentrantLock()
    particlecounter = 0
    @everywhere function update_particlecounter()
        global particlecounter
        particlecounter += 1
    end

    if ensemble_algorithm == EnsembleDistributed()
        # Distributed
        # Overwrite the the ensemble problem keyword argument `prob_func` with some
        # progress verbose.
        prob_func = (prob, i, repeat) -> begin
            #remote_do(update_particlecounter, 1) # 1 is main process
            #ensembleprob_kwargs.prob_func(prob, i, repeat)
            prob_func_in(prob, i, repeat)
        end
    elseif ensemble_algorithm == EnsembleThreads()
        prob_func = (prob, i, repeat) -> begin
            #global particlecounter
            #global lk
            lock(lk)
            try
                particlecounter += 1
            finally
                unlock(lk)
            end
            prob_func_in(prob, i, repeat)
        end
    else
        # Serial
        # Overwrite the the ensemble problem keyword argument `prob_func` with some
        # progress verbose.
        prob_func = (prob, i, repeat) -> begin
            #print("Particle $i  of $npart \r")
            print_particlecounter(i, npart)
            prob_func_in(prob, i, repeat)
        end
    end



    #===============================================================================#

    #PROBLEM SETUP
    #
    """
    "Sketch" the shape of the ODE-problem. This will be remade by the problem_func
    for each particle in the EnsembleProblem, but the particle independent `eom`,
    `tspan` and `charge`, `mass`, and `fields` may be used.
    """
    prob = ODEProblem(
        eom,
        zeros(4), # Dummy initial condition
        tspan,
        (
            charge=charge,
            mass=mass,
            fields=fields_itp,
            magneticmoment=0.0
        )
    )
    ensemble_prob = EnsembleProblem(
        prob
        ;
        prob_func=prob_func,
        ensembleprob_kwargs...
    )

    #===============================================================================#
    # START OF SIMULATION

    println("Running ensemble of $npart particles...")
    println(".................................................................",
        "...............")

    running = true
    if ensemble_algorithm != EnsembleSerial()
        @async begin
            while running
                #print("Particle $particlecounter of $npart \r")
                print_particlecounter(particlecounter, npart)
                sleep(0.1)
            end
        end
    end

    @time sim = DifferentialEquations.solve(
        ensemble_prob,
        solve_args...
        ;
        solve_kwargs...
    )

    running = false

    if !(:reduction in keys(ensembleprob_kwargs))
        @warn "No reduction function specified.
        Saving the full ensemble simulation as 'out.jld2'."
        JLD2.@save "out.jld2" sim
    else
        try
            println("\nPost-processing: Computing GCAStates...")
            filename = tp.get_filename(ensembleprob_kwargs.reduction)
            @time tp.save_gcastates(filename, fields_itp)
            println("                 Success")
        catch e
            println("                 Failed:\n", e)
        end
    end
    return ensemble_algorithm
end
"""
    rerun_nonthermals(fname, nparticles)
Rerun the `nparticles` most energetic particles from the experiment defined by
`fname`.
"""
function rerun_nonthermals(
    fname::String,
    nparticles::Int
    ;
    kwargs...
)
    u0, mu0, _ = find_nonthermals(fname, nparticles)
    return rerun(
        u0,
        mu0,
        nparticles,
        fname;
        kwargs...
    )
end


"""
    rerun_maxiters(fname, nparticles)
Rerun the `nparticles` particles with most timesteps from the experiment
defined by `fname`.
"""
function rerun_maxiters(
    fname::String,
    nparticles::Int
    ;
    kwargs...
)
    u0, mu0, _ = find_maxiters(fname, nparticles)
    return rerun(
        u0,
        mu0,
        nparticles,
        fname;
        kwargs...
    )
end


"""
    rerun_idxs(fname, idxs)
Rerun the the test particles of index `idxs` from the experiment defined by
`fname`.
"""
function rerun_idxs(
    fname::String,
    idxs::AbstractVector
    ;
    kwargs...
)
    u0, mu0 = h5_getinitialconditions(fname)
    u0 = u0[idxs]
    mu0 = mu0[idxs]
    nparticles = length(idxs)
    return rerun(
        u0,
        mu0,
        nparticles,
        fname;
        kwargs...
    )
end


"""
    rerun(u0, mu0, nparticles, fname; kwargs)
Reruns test particles in the experiment defined by `fname`. The initial
positions and velocities of the particles are arbitrary, defined in `u0`,
and their magnetic moments are defined in `mu0`. Some parameters of the
experiment may be adjusted through keyword arguments. These are

# Keyword arguments
- `reltol`
- `abstol`
- `maxiters`
- `alg`: solver algorithm
- `gcatol`: Tolerance for the guiding centre approximation.
"""
function rerun(
    u0::Vector{<:Vector{<:Real}},
    mu0::Vector{<:Real},
    nparticles::Int,
    fname::String;
    reltol=nothing,
    abstol=nothing,
    maxiters=nothing,
    alg=Tsit5(),
    gcatol=nothing,
)
    include(fname * ".jl")
    tspans = [tspan for _ in 1:nparticles]
    if !isnothing(gcatol)
        solve_kwargs[:callback] = CallbackSet(
            solve_kwargs[:callback].discrete_callbacks[1],
            GCABreakDownCB("lowmem", "2Dxz", tolerance=gcatol)
        )
    end
    s_kwargs = (
        reltol=isnothing(reltol) ? solve_kwargs[:reltol] : reltol,
        abstol=isnothing(abstol) ? solve_kwargs[:abstol] : abstol,
        maxiters=isnothing(maxiters) ? solve_kwargs[:maxiters] : maxiters,
        trajectories=nparticles,
        callback=solve_kwargs[:callback],
    )

    prob = ODEProblem(
        eom,
        u0[1],
        tspans[1],
        (
            charge=charge,
            mass=mass,
            fields=fields_itp,
            magneticmoment=mu0[1],
        )
    )
    ensemble_prob = EnsembleProblem(
        prob
        ;
        prob_func=tp.PposPvel(u0, mu0, tspans),
        # Use default output function. I.e. return the whole solution
        safetycopy=false,
    )

    solve_args = (alg, EnsembleSerial())
    @info "Re-running $fname with $nparticles particles..."
    @time sol = DifferentialEquations.solve(
        ensemble_prob,
        solve_args...
        ;
        s_kwargs...
    )
    return sol
end

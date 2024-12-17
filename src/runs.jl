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
function rerun_nonthermals(
    fname::String,
    nparticles::Int
    ;
    relativegain::Bool=false,
    reltol=nothing,
    abstol=nothing,
    maxiters=10_000,
)
    u0, mu0, _ = find_nonthermals(fname, nparticles; relativegain=relativegain)
    return rerun(
        u0,
        mu0,
        nparticles,
        fname;
        reltol=reltol,
        abstol=abstol,
        maxiters=maxiters,
    )
end


function rerun_maxiters(
    fname::String,
    nparticles::Int
    ;
    reltol=nothing,
    abstol=nothing,
    maxiters=10_000,
)
    u0, mu0, _ = find_maxiters(fname, nparticles)
    return rerun(
        u0,
        mu0,
        nparticles,
        fname;
        reltol=reltol,
        abstol=abstol,
        maxiters=maxiters,
    )
end


function rerun(
    u0::Vector{<:Vector{<:Real}},
    mu0::Vector{<:Real},
    nparticles::Int,
    fname::String;
    reltol=nothing,
    abstol=nothing,
    maxiters=10_000,
)
    include(fname * ".jl")
    tspans = [tspan for _ in 1:nparticles]
    s_kwargs = (
        reltol=isnothing(reltol) ? solve_kwargs.reltol : reltol,
        abstol=isnothing(abstol) ? solve_kwargs.abstol : abstol,
        maxiters=maxiters,
        trajectories=nparticles,
        callback=solve_kwargs.callback
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
    )

    solve_args = (Tsit5(), EnsembleSerial())
    @info "Re-running $fname with $nparticles particles..."
    @time sol = DifferentialEquations.solve(
        ensemble_prob,
        solve_args...
        ;
        s_kwargs...
    )
    return sol
end


function convergencetest(
    u0,
    mu0,
    charge,
    mass,
    fields_itp,
    eom,
    tspan,
    callback;
    alg=Tsit5(),
    reltol=10 .^ LinRange(-14, -3, 12),
    abstol=10 .^ LinRange(-16, -5, 12),
    maxiters=100_000_000
)
    solve_args = (alg, EnsembleSerial())
    npart = length(u0)
    ntols = length(reltol)
    if ntols != length(abstol)
        throw(ArgumentError("reltol and abstol must have the same length"))
    end
    endstates = Matrix{Vector{Float64}}(undef, ntols, npart)
    nt = Matrix{Int64}(undef, ntols, npart)
    tspans = [tspan for _ in 1:npart]
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
        output_func=(sol, i) -> (
            (
                sol.u[end],
                length(sol.t),
            ),
            false
        )
    )
    @time for i in 1:ntols
        solve_kwargs = (
            reltol=reltol[i],
            abstol=abstol[i],
            maxiters=maxiters,
            trajectories=npart,
            callback=callback
        )
        @info "Running with tolerances reltol=$(reltol[i]), abstol=$(abstol[i])"
        @time sol = DifferentialEquations.solve(
            ensemble_prob,
            solve_args...
            ;
            solve_kwargs...
        )
        for j in 1:npart
            endstates[i, j] = sol[j][1]
            nt[i, j] = sol[j][2]
        end
    end
    return endstates, nt, reltol, abstol
end

function find_nonthermals(
    fname::String,
    nparticles::Int;
    relativegain::Bool=false,
)
    batches = h5open(fname * ".h5", "r") do fid
        1:length(fid)
    end
    x0 = h5_getdataset(fname * ".h5", "x0"; batches=batches)
    y0 = h5_getdataset(fname * ".h5", "y0"; batches=batches)
    z0 = h5_getdataset(fname * ".h5", "z0"; batches=batches)
    vparal0 = h5_getdataset(fname * ".h5", "vparal0"; batches=batches)
    mu0 = h5_getdataset(fname * ".h5", "magneticmoment"; batches=batches)
    e0, ef = h5_getenergies(fname; batches=batches)

    sortby = relativegain ? (ef .- e0) ./ e0 : ef
    idxs = sortperm(sortby)[end-nparticles+1:end]
    u0 = [[x0[i], y0[i], z0[i], vparal0[i]] for i in idxs]
    mu0 = mu0[idxs]

    return u0, mu0, idxs
end


function find_maxiters(
    fname::String,
    nparticles::Int;
)
    batches = h5open(fname * ".h5", "r") do fid
        1:length(fid)
    end
    x0 = h5_getdataset(fname * ".h5", "x0"; batches=batches)
    y0 = h5_getdataset(fname * ".h5", "y0"; batches=batches)
    z0 = h5_getdataset(fname * ".h5", "z0"; batches=batches)
    vparal0 = h5_getdataset(fname * ".h5", "vparal0"; batches=batches)
    mu0 = h5_getdataset(fname * ".h5", "magneticmoment"; batches=batches)

    retmsg = h5_getdataset(fname * ".h5", "retmsg"; batches=batches)
    nt = h5_getdataset(fname * ".h5", "nt"; batches=batches)

    sortby = nt
    idxs = sortperm(sortby)[end-nparticles+1:end]
    u0 = [[x0[i], y0[i], z0[i], vparal0[i]] for i in idxs]
    mu0 = mu0[idxs]

    return u0, mu0, idxs
end


function h5_getinitialconditions(
    fname::String;
    mask=nothing
)
    batches = h5open(fname * ".h5", "r") do fid
        1:length(fid)
    end
    x0 = h5_getdataset(fname * ".h5", "x0"; batches=batches)
    y0 = h5_getdataset(fname * ".h5", "y0"; batches=batches)
    z0 = h5_getdataset(fname * ".h5", "z0"; batches=batches)
    vparal0 = h5_getdataset(fname * ".h5", "vparal0"; batches=batches)
    mu0 = h5_getdataset(fname * ".h5", "magneticmoment"; batches=batches)
    u0 = [[x0[i], y0[i], z0[i], vparal0[i]] for i in eachindex(x0)]
    if isnothing(mask)
        return u0, mu0
    else
        return u0[mask], mu0[mask]
    end
end

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

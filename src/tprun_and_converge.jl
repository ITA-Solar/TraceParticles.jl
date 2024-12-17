using HDF5
using Distributed

params_name = ARGS[1]
if length(ARGS) > 1
    np = parse(Int, ARGS[2])
    for i in 1:np
        addprocs(1)
    end
end
paramsfile = joinpath(pwd(), basename(params_name))
include(paramsfile)
try
    mkdir("data")
catch
end

requirements = [
    :expname,
    :eom,
    :charge,
    :mass,
    :fields_itp,
    :tspan,
]
for i in requirements
    try
        eval(i)
    catch
        error("The parameter file must define the variable `$i`")
    end
end

if !isdefined(Main, :expname)
    error("")
end

datadir = try
    datadir
catch
    println("'datadir' not set, using 'data'.")
    "data"
end

expname *= "_convergence"
ensembleprob_kwargs.reduction.expname = expname
expdir = datadir * "/" * expname
paramsbackup = expdir * "/" * expname * ".jl"
try
    mkdir(expdir)
    cp(paramsfile, paramsbackup)
catch
    @warn "The directory $(pwd() * "/" * expdir) already exists.
     Do you want to delete it and re-run? y/n"
    if readline() == "y"
        rm(expdir, recursive=true)
        mkdir(expdir)
        cp(paramsfile, paramsbackup)
    else
        error("Exiting.")
    end
end

#------------------------------------------------------------------------------

# Initial run
run = 1
ensembleprob_kwargs.reduction.postfix = "_$run"
# With real values; convergence at 1e8, 1e-8, 1e-3
solve_kwargs[:maxiters] = 1_000
solve_kwargs[:reltol] = 0.3e-3#[0.3e-4, 0.3e-4, 0.3e-4, 1e-3]
solve_kwargs[:abstol] = 1e4#[1e3, 1e3, 1e3, 1e3]
tolerance_decrement = 0.1

@time ensemble_algorithm = tp.tprun(
    solve_kwargs,
    ensembleprob_kwargs,
    eom,
    npart,
    tspan,
    charge,
    mass,
    fields_itp,
    prob_func,
)


filename_original = tp.get_filename(ensembleprob_kwargs.reduction)

batches, nbatches = tp.h5_nbatches(filename_original)
h5open(filename_original, "cw") do fid
    for b in batches
        h5group = fid["batch_$b"]
        N = length(h5group["x0"])
        write(h5group, "reltol", solve_kwargs[:reltol] * ones(N))
        write(h5group, "abstol", solve_kwargs[:abstol] * ones(N))
    end
end

efprev = tp.h5_getdataset(filename_original[1:end-3] * "_gcastatesf.h5", "energy")
xfprev = tp.h5_getdataset(filename_original, "xf")
yfprev = tp.h5_getdataset(filename_original, "yf")
zfprev = tp.h5_getdataset(filename_original, "zf")

rerun = 1:npart

#reltol = [0.3e-5, 0.3e-5, 0.3e-5, .3e-5]
#abstol = [1e2, 1e2, 1e2, 1e2]

@time while length(rerun) > 0
    global rerun
    global solve_kwargs
    local filename = tp.get_filename(ensembleprob_kwargs.reduction)
    u0, mu0 = tp.h5_getinitialconditions(filename_original[1:end-3], mask=rerun)

    global run += 1
    global alg
    global ensemble_algorithm
    solve_kwargs[:maxiters] /= tolerance_decrement
    #global reltol[4] *= tolerance_decrement
    #global abstol[4] *= tolerance_decrement
    solve_kwargs[:reltol] *= tolerance_decrement
    solve_kwargs[:abstol] *= tolerance_decrement
    solve_kwargs[:trajectories] = length(rerun)
    num_rerunpart = length(rerun)
    if num_rerunpart < solve_kwargs[:batch_size]
        solve_kwargs[:batch_size] = num_rerunpart
    end
    global ensembleprob_kwargs
    ensembleprob_kwargs.reduction.batchnr = 0
    ensembleprob_kwargs.reduction.postfix = "_$run"
    @info """Run $run
        Rerunning  $num_rerunpart particles.
        Maxiters = $(solve_kwargs[:maxiters])
        Reltol   = $(solve_kwargs[:reltol])
        Abstol   = $(solve_kwargs[:abstol])
    """

    tspans = [tspan for i in 1:num_rerunpart]
    ep_kwargs = (
        ensembleprob_kwargs...,
        prob_func=tp.PposPvel(u0, mu0, tspans)
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
        ep_kwargs...
    )
    running = true
    ensemble_algorithm = EnsembleSerial()
    if ensemble_algorithm !== EnsembleSerial()
        @async begin
            while running
                #print("Particle $particlecounter of $npart \r")
                print_particlecounter(particlecounter, npart)
                sleep(0.1)
            end
        end
    end
    @time sol = DifferentialEquations.solve(
        ensemble_prob,
        alg,
        ensemble_algorithm,
        ;
        solve_kwargs...
    )

    filename = tp.get_filename(ensembleprob_kwargs.reduction)
    if !(:reduction in keys(ensembleprob_kwargs))
        @warn "No reduction function specified.
        Saving the full ensemble simulation as 'out.jld2'."
        JLD2.@save "out.jld2" sim
    else
        try
            println("\nPost-processing: Computing GCAStates...")
            @time tp.save_gcastates(filename, fields_itp)
            println("                 Success")
        catch e
            println("                 Failed:\n", e)
        end
    end
    # Check for convergence
    @info filename
    retcode = tp.h5_getdataset(filename, "retcode")
    local ef = tp.h5_getdataset(filename[1:end-3] * "_gcastatesf.h5", "energy")
    local xf = tp.h5_getdataset(filename, "xf")
    local yf = tp.h5_getdataset(filename, "yf")
    local zf = tp.h5_getdataset(filename, "zf")
    x0 = tp.h5_getdataset(filename, "x0")
    y0 = tp.h5_getdataset(filename, "y0")
    z0 = tp.h5_getdataset(filename, "z0")

    global efprev
    global xfprev
    global yfprev
    global zfprev

    energydiff = ef .- efprev
    rel_energydiff = energydiff ./ efprev
    distdiff = Vector{Float64}(undef, num_rerunpart)
    startpos = Vector{Float64}(undef, num_rerunpart)
    rerun_new = Vector{Bool}(undef, num_rerunpart)
    @. startpos = sqrt(x0^2 + y0^2 + z0^2)
    @. distdiff = sqrt((xf - xfprev)^2 + (yf - yfprev)^2 + (zf - zfprev)^2)
    rel_distdiff = distdiff ./ startpos
    @. rerun_new =
        (rel_energydiff > 1e-2) |
        (rel_distdiff > 1e-2) |
        (retcode == "MaxIters") |
        (retcode == "DtLessThanMin")

    # The indices of the particles that has converged
    rerun_converged_idxs = findall(x -> x == false, rerun_new)
    # The indices in the original datafile of the particles that has converged
    update_original_data_idxs = rerun[.!rerun_new]
    # The indices in the original datafile of the particles that did not converge
    rerun = rerun[rerun_new]

    local batches
    local nbatches
    batches, nbatches = tp.h5_nbatches(filename)
    h5open(filename, "cw") do fid
        for b in batches
            h5group = fid["batch_$b"]
            N = length(h5group["x0"])
            write(h5group, "reltol", solve_kwargs[:reltol] * ones(N))
            write(h5group, "abstol", solve_kwargs[:abstol] * ones(N))
            write(h5group, "maxiters", solve_kwargs[:maxiters] * ones(N))
        end
        h5group = create_group(fid, "rerun_idxs")
        write(h5group, "did_converge", rerun_converged_idxs)
        write(h5group, "original_idxs", update_original_data_idxs)
    end

    efprev = ef[rerun_new]
    xfprev = xf[rerun_new]
    yfprev = yf[rerun_new]
    zfprev = zf[rerun_new]
end

#h5_update_with_rerunned_particles!(
#    filename_original,
#    filename_rerun,
#    update_original_data_idxs,
#    rerun_converged_idxs,
#    reltol,
#    abstol,
#    batchsize
#)

function h5_update_with_rerunned_particles!(
    filename_original,
    filename_rerun,
    update_original_data_idxs,
    rerun_converged_idxs,
    reltol,
    abstol,
    batchsize
)
    nupdates = length(update_original_data_idxs)
    if nupdates !== length(rerun_converged_idxs)
        throw(
            ArgumentError(
                "The number of particles in the original datafile that should be updated
                is not equal to the number of particles that has converged."
            )
        )
    end
    fid_og = h5open(filename_original, "cw")
    fid_rerun = h5open(filename_rerun, "r")
    for i in 1:nupdates
        idx_og = update_original_data_idxs[i]
        idx_rerun = rerun_converged_idxs[i]
        batchnr_og = ceil(Int, idx_og / batchsize)
        batchnr_rerun = ceil(Int, idx_rerun / batchsize)
        for j in 1:nfields
        end
    end
    close(fid_og)
    close(fid_rerun)
end

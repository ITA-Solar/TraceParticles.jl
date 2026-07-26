using Dates
using Distributed
using DiffEqCallbacks
using OrdinaryDiffEq
using Printf
using SciMLBase

"""
This script is used as a fast way to run a single test particle simulation
from the command line that adheres to my data management plan.

It is used as follows:
    `julia --proj=path/to/projectenv -- tprun.jl expparams.jl [parallelisation]`

where `expparams.jl` is the file containing the parameters of the simulation and
`parallelisation` is an optional argument that decided where the simulation is
run in serial or using threads or `Distributed.jl`. Currently, only serial
parallelisation is supported.

The script copies the `expparams.jl` file to a backup location in the directory
`data`, which it creates if it doesn't exist. The location of the actual
simulation data depends on the reduction function chosen by the user.
chosen by the user.

## Requirements
As a bare minimum, the `expparams.jl` file must define the following variables:
- `expname::String` : The name of the test particle experiment.
- `eom::Function` : The equation of motion of the test particles.
- `fields_itp::Vector{<:AbstractInterpolation} : The interpolated fields
    that the test particles will be evolved in.
- `charger::Real` : The charge of the test particles.
- `mass::Real` : The mass of the test
- `tspan::Tuple{Real, Real}` : The time span of the simulation.
- `ensembleprob_kwargs::Dict` : Keyword arguments for the `EnsembleProblem`.
- `solve_args::Dict` : Arguments for the `solve`-function.
- `solve_kwargs::Dict` : Keyword arguments for the `solve`-function.

## Positional and keyword arguments to `EnsembleProblem` and `solve`
Keyword arguments to the `EnsembleProblem` constructor and  keyword arguments
to the `solve`-function must be defined in `Dict`s.

The dictionaries are then splashed into place.

E.g. to define a output function, reduction, maxiters, safetycopy,
and solution algorithm, one must, in the user-provided parameter-file,
define the following:
```julia
ensembleprob_kwargs = Dict(
    :prob_func => myprob_func,
    :output_func => myoutput_func,
    :reduction => myreduction
)

solve_kwargs = Dict(
    :trajectories => 10,
    :maxiters => 10000,
    :safetycopy => false,
    )
```

## Parallelisation
To multithread, provide the number of threads through command line.

To run using distributed computing, provide the number of processes through
the command line.

If the script is run with multiple threads and mulitple processes, the script
will use distributed computing.
"""

localARGS = isdefined(Main, :newARGS) ? newARGS : ARGS

commandlineparams = localARGS[1]
paramsfile = joinpath(pwd(), basename(commandlineparams))
if !isdefined(Main, :paramset)
    include(paramsfile)
end

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

if !isdefined(Main, :datadir)
    println("'datadir' not set, using 'data'.")
    datadir = "data"
end

expdir = datadir * "/" * expname
paramsbackup = expdir * "/" * expname * ".jl"
try
    mkdir(expdir)
    cp(paramsfile, paramsbackup)
catch
    println("The directory $(pwd() * "/" * expdir) already exists.
     Do you want to delete it and re-run? y/n")
    if readline() == "y"
        rm(expdir, recursive=true)
        mkdir(expdir)
        cp(paramsfile, paramsbackup)
    else
        error("Exiting.")
    end
end

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
    try
        u0
    catch
        zeros(1) # Dummy initial condition
    end,
    typeof(tspan) <: Tuple{Real,Real} ? tspan : (0, 1),
    if eom == guidingcentreapproximation!
        GCAParams(
            charge=charge,
            mass=mass,
            electromagneticfield=fields_itp,
            magneticmoment=try
                mu0
            catch
                0.0
            end,
        )
    elseif eom == hybridgcafo!
        HybridParams(
            charge=charge,
            mass=mass,
            electromagneticfield=fields_itp,
            magneticmoment=try
                mu0
            catch
                0.0
            end,
        )
    end
)
ensemble_prob = EnsembleProblem(
    prob
    ;
    ensembleprob_kwargs...
)

#===============================================================================#
# START OF SIMULATION

println("Running program: $PROGRAM_FILE")
println("Start time: $(string(now()))")
println("Host name : $(gethostname())")
println("Running ensemble of $npart particles...")
println(".................................................................",
    "...............")

@time sim = solve(
    ensemble_prob,
    solve_args...
    ;
    solve_kwargs...
);

running = false

if !(:reduction in keys(ensembleprob_kwargs))
    @warn "No reduction function specified.
    Saving the full ensemble simulation as 'out.jld2'."
    JLD2.@save "out.jld2" sim
else
    try
        println("\nPost-processing: Computing energies...")
        filename = get_filename(ensembleprob_kwargs[:reduction])
        @time save_energy(
            filename,
            fields_itp,
            ensemble_prob.prob.f.f
        )
        println("                 Success")
    catch e
        println("                 Failed:\n", e)
    end
end

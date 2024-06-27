using DifferentialEquations
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
- `tp_expname::String` : The name of the test particle experiment.
- `eom::Function` : The equation of motion of the test particles.
- `fields_itp::Vector{<:AbstractInterpolation} : The interpolated fields
    that the test particles will be evolved in.
- `charger::Real` : The charge of the test particles.
- `mass::Real` : The mass of the test
- `tspan::Tuple{Real, Real}` : The time span of the simulation.
- `prob_func::Function` or `u0::Vector{<:Real}`: The function that generates
    the initial conditions of the test particles or an initial condition to be
    used for all test particles. `prob_func` must be encapsulated in the
    `NamedTuple` `ensembleprob_kwargs`.
- `trajectories::Int` : The number of test particles to simulate. Must be
    encapsulated in the `NamedTuple` `solve_kwargs`.

## Passing positoinal and keyword arguments to EnsembleProblem and Solve
Keyword arguments to the `EnsembleProblem` constructor and  keyword arguments
to the `solve` function must be defined in separate `NamedTuple`s. The tuples
must be named `ensembleprob_kwargs` and `solve_kwargs`.

The tuples are then splashed into place.

E.g. to define a output function, reduction, maxiters, safetycopy,
and solution algorithm, one must, in the user-provided parameter-file,
define the following:
```julia
ensembleprob_kwargs = (
    prob_func = myprob_func,
    output_func = myoutput_func,
    reduction = myreduction
)

solve_kwargs = (
    trajectories = 10,
    maxiters = 10000,
    safetycopy = false,
    )
```

## Parlallelisation
To multithread, provide the number of threads through command line.

To run using distributed computing, provide the number of processes through
the command line.

If the script is run with multiple threads and mulitple processes, the script
will use distributed computing.
"""
params = ARGS[1]
include(joinpath(pwd(), basename(params)))

try mkdir("data")
catch
end

try expname
catch
    error("'expname' undefined. The parameter file must define the experiment"*
        " name in the parameters file")
end

datadir = try datadir
catch
    pritntln("'datadir' not set, using 'data'.")
    "data"
end

expdir = datadir * "/" * expname
paramsbackup = expdir * "/" * expname * ".jl"
try mkdir(expdir)
    cp(params, paramsbackup)
catch
    println("The directory $(pwd() * "/" * expdir) already exists.
     Do you want to delete it and re-run? y/n")
    if readline() == "y"
        rm(expdir, recursive=true)
        mkdir(expdir)
        cp(params, paramsbackup)
    else
        error("Exiting.")
    end
end

ensemble_algorithm = try np = nprocs()
    if np > 1
        println("Running on $np processes.")
        EnsembleDistributed()
    end
catch
    if Threads.nthreads() > 1
        println("Running on $(Threads.nthreads()) threads.")
        EnsembleThreads()
    else
        println("Running in serial.")
        EnsembleSerial()
    end
end
solve_args = try (alg, ensemble_algorithm)
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
    try u0
    catch
        zeros(1) # Dummy initial condition
    end,
    tspan,
    (
        charge = charge,
        mass = mass,
        fields = fields_itp,
    )
)
ensemble_prob = EnsembleProblem(
        prob
        ;
        ensembleprob_kwargs...
    )

#===============================================================================#
# START OF SIMULATION

println("Running ensemble of $npart particles...")
@time sim = DifferentialEquations.solve(
    ensemble_prob,
    solve_args...
    ;
    solve_kwargs...
    );

if !(:reduction in keys(ensembleprob_kwargs))
    @warn "No reduction function specified.
    Saving the full ensemble simulation as 'out.jld2'."
    JLD2.@save "out.jld2" sim
end

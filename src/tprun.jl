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
    used for all test particles.
"""
params = ARGS[1]
include(joinpath(pwd(), basename(params)))
try mkdir("data")
catch
end
paramsbackup = pwd() * "/data/" * tp_expname * ".jl"
cp(params, paramsbackup)

# Run simulation
try parallelisation = ARGS[2] 
    if parallelisation == "serial"
        include("runserial.jl")
    else 
        error("Only serial parallelisation is supported.")
    end
catch
    include("runserial.jl")
end


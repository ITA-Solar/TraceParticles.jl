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
- `tp_expname::String` : The name of the experiment.
"""
params = ARGS[1]
include(joinpath(pwd(), basename(params)))
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


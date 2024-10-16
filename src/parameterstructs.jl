"""
    parameterstructs.jl
This file defines mutable structs that are useful as containers for the
parameters in a `ODEProblem` in DifferentialEquations.jl. With mutable
structs, the parameters can be altered during the integration, as opposed
to using a `Tuple` or a `NamedTuple`.
"""

"""
    HybridParams
Parameter container for a hybrid GCA/FO particle simulation. The `switch`
determines which EoM to use and the random number generator `rng` is used to
determine the particle phase when switching from GCA to full orbit integration.
"""
mutable struct HybridParams
    charge::Real
    mass::Real
    magneticmoment::Real
    fields
    rng::AbstractRNG
    switch::Int
end

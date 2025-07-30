"""
    parameterstructs.jl
This file defines mutable structs that are useful as containers for the
parameters in a `ODEProblem` in DifferentialEquations.jl. With mutable
structs, the parameters can be altered during the integration, as opposed
to using a `Tuple` or a `NamedTuple`.
"""

mutable struct ParticleParams{T1<:Real,T2}
    charge::T1
    mass::T1
    electromagneticfield::T2

    function ParticleParams(
        ; charge::T1, mass::T1, electromagneticfield::T2,
    ) where {T1<:Real,T2}
        new{T1,T2}(charge, mass, electromagneticfield)
    end
end

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

mutable struct GCAParams{T1<:Real,T2,T3<:Int}
    charge::T1
    mass::T1
    electromagneticfield::T2
    magneticmoment::T1
    weight::T1
    nrejections::T3
    retcode::T3

    function GCAParams(
        ;
        charge::T1,
        mass::T1,
        electromagneticfield::T2,
        magneticmoment::T1,
        weight::T1,
        nrejections::T3,
        retcode::T3,
    ) where {T1<:Real,T2,T3<:Int}
        new{T1,T2,T3}(
            charge,
            mass,
            electromagneticfield,
            magneticmoment,
            weight,
            nrejections,
            retcode,
        )
    end
end

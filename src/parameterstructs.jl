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
    GCAParams
A container for the parameters used in a guiding centre approximation
simulation where the initial state of the particles are drawn using rejection
sampling and they have a statistical weight. The struct also contains a
`terminationcode` that indicates what kind of callback terminated the particle
trajectory.
"""
mutable struct GCAParams{T1<:Real,T2,T3<:Int,T4}
    charge::T1
    mass::T1
    electromagneticfield::T2
    magneticmoment::T1
    weight::T1
    nrejections::T3
    terminationcode::T4

    function GCAParams(
        ;
        charge::T1,
        mass::T1,
        electromagneticfield::T2,
        magneticmoment::T1,
        weight::T1=1.0,
        nrejections::T3=0,
        terminationcode::T4=TerminationCode.NotTerminated,
    ) where {T1<:Real,T2,T3<:Int,T4}
        new{T1,T2,T3,T4}(
            charge,
            mass,
            electromagneticfield,
            magneticmoment,
            weight,
            nrejections,
            terminationcode,
        )
    end
end

"""
    HybridParams
Parameter container for a hybrid GCA/FO particle simulation. The `eomid`
determines which EoM to use and the random number generator `rng` is used to
determine the particle phase when switching from GCA to full orbit integration.
"""
mutable struct HybridParams{
    T1<:Real,
    T2,
    T3<:Int,
    T4,
    T5<:AbstractRNG,
    T6<:Function,
}
    charge::T1
    mass::T1
    electromagneticfield::T2
    magneticmoment::T1
    weight::T1
    nrejections::T3
    terminationcode::T4
    rng::T5
    getphase::T6
    eomid::T3
    nswitches::T3
    initialeomid::T3

    function HybridParams(
        ;
        charge::T1,
        mass::T1,
        electromagneticfield::T2,
        magneticmoment::T1=0.0,
        weight::T1=1.0,
        nrejections::T3=0,
        terminationcode::T4=TerminationCode.NotTerminated,
        rng::T5=Xoshiro(1),
        getphase::T6=(integrator) -> rand(
            integrator.p.rng, 0.0, Float64(2pi)
        ),
        nswitches::T3=0,
        initialeomid::T3=2,
    ) where {T1<:Real,T2,T3<:Int,T4,T5<:AbstractRNG,T6<:Function}
        new{T1,T2,T3,T4,T5,T6}(
            charge,
            mass,
            electromagneticfield,
            magneticmoment,
            weight,
            nrejections,
            terminationcode,
            rng,
            getphase,
            initialeomid,
            nswitches,
            initialeomid,
        )
    end
end

"""
    HybridParamsWithDetection
Parameter container for a hybrid GCA/FO particle simulation. It is similar to
`HybridParams` but with switch times and respective EoM-IDs to use for post
processing. It also has a `maximum` number of switches field.
"""
mutable struct HybridParamsWithDetection{
    T1<:Real,
    T2,
    T3<:Int,
    T4,
    T5<:AbstractRNG,
    T6<:Function,
    T7<:AbstractVector,
    T8<:AbstractVector
}
    charge::T1
    mass::T1
    electromagneticfield::T2
    magneticmoment::T1
    weight::T1
    nrejections::T3
    terminationcode::T4
    rng::T5
    getphase::T6
    eomid::T3
    nswitches::T3
    initialeomid::T3
    initialmagneticmoment::T1
    maxswitches::T3
    timeatswitch::T7
    eomidafterswitch::T8
    magneticmomentafterswitch::T7

    function HybridParamsWithDetection(
        ;
        charge::T1,
        mass::T1,
        electromagneticfield::T2,
        magneticmoment::T1=0.0,
        weight::T1=1.0,
        nrejections::T3=0,
        terminationcode::T4=TerminationCode.NotTerminated,
        rng::T5=Xoshiro(1),
        getphase::T6=(integrator) -> rand(
            integrator.p.rng, 0.0, Float64(pi)
        ),
        nswitches::T3=0,
        initialeomid::T3=2,
        initialmagneticmoment::T1=magneticmoment,
        maxswitches::T3=100,
    ) where {T1<:Real,T2,T3<:Int,T4,T5<:AbstractRNG,T6<:Function}
        timeatswitch = zeros(Float64, maxswitches)
        eomidafterswitch = zeros(Int32, maxswitches)
        magneticmomentafterswitch = zeros(Int32, maxswitches)
        new{T1,T2,T3,T4,T5,T6,Vector{Float64},Vector{Int32}}(
            charge,
            mass,
            electromagneticfield,
            magneticmoment,
            weight,
            nrejections,
            terminationcode,
            rng,
            getphase,
            initialeomid,
            nswitches,
            initialeomid,
            initialmagneticmoment,
            maxswitches,
            timeatswitch,
            eomidafterswitch,
            magneticmomentafterswitch
        )
    end
end

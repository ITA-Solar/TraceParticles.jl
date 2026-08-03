"""
This file defines mutable structs that are useful as containers for the
parameters in a `ODEProblem` in from SciML. With mutable
structs, the parameters can be altered during the integration, as opposed
to using a `Tuple` or a `NamedTuple`.
"""

"""
    FullOrbitParams
A container for the parameters used in a particle simulation solving the
Lorentz force directly, and the initial state of the particles are drawn
using rejection sampling giving them a statistical weight. The struct also
contains a `terminationcode` that indicates what kind of callback terminated
the particle trajectory.
"""
@kwdef mutable struct FullOrbitParams{T1<:Real,T2,T3<:Int,T4,T5}
    charge::T1
    mass::T1
    electromagneticfield::T2
    weight::T5 = 1.0
    nrejections::T3 = 0
    terminationcode::T4 = TerminationCode.NotTerminated
end

"""
    GCAParams
A container for the parameters used in a guiding centre approximation
simulation where the initial state of the particles are drawn using rejection
sampling and they have a statistical weight. The struct also contains a
`terminationcode` that indicates what kind of callback terminated the particle
trajectory.
"""
@kwdef mutable struct GCAParams{T1<:Real,T2,T3<:Int,T4,T5}
    charge::T1
    mass::T1
    electromagneticfield::T2
    magneticmoment::T1
    weight::T5 = 1.0
    nrejections::T3 = 0
    terminationcode::T4 = TerminationCode.NotTerminated
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
    T4<:TerminationCode.T,
    T5<:AbstractRNG,
    T6<:Function,
    T7<:EoMID.T,
}
    charge::T1
    mass::T1
    electromagneticfield::T2
    initialmagneticmoment::T1
    magneticmoment::T1
    weight::T1
    nrejections::T3
    terminationcode::T4
    rng::T5
    getphase::T6
    eomid::T7
    nswitches::T3
    initialeomid::T7

    function HybridParams(
        ;
        charge::T1,
        mass::T1,
        electromagneticfield::T2,
        magneticmoment::T1=NaN,
        weight::T1=1.0,
        nrejections::T3=0,
        initialeomid::T7=EoMID.GuidingCentreApproximation,
        rng::T5=Xoshiro(),
        terminationcode::T4=TerminationCode.NotTerminated,
        getphase::T6=(integrator) -> 2pi*rand(integrator.p.rng),
        nswitches::T3=0,
    ) where {T1<:Real,T2,T3<:Int,T4<:TerminationCode.T,T5<:AbstractRNG,T6<:Function,T7<:EoMID.T}
        new{T1,T2,T3,T4,T5,T6,T7}(
            charge,
            mass,
            electromagneticfield,
            magneticmoment,
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
    T4<:TerminationCode.T,
    T5<:AbstractRNG,
    T6<:Function,
    T7<:AbstractVector,
    T8<:Vector{<:Enum{Int32}},
    T9<:EoMID.T
}
    charge::T1
    mass::T1
    electromagneticfield::T2
    initialmagneticmoment::T1
    magneticmoment::T1
    weight::T1
    nrejections::T3
    terminationcode::T4
    rng::T5
    getphase::T6
    eomid::T9
    nswitches::T3
    initialeomid::T9
    maxswitches::T3
    timeatswitch::T7
    eomidafterswitch::T8
    magneticmomentafterswitch::T7

    function HybridParamsWithDetection(
        ;
        charge::T1,
        mass::T1,
        electromagneticfield::T2,
        magneticmoment::T1=NaN,
        weight::T1=1.0,
        nrejections::T3=0,
        terminationcode::T4=TerminationCode.NotTerminated,
        rng::T5=Xoshiro(1),
        getphase::T6=(integrator) -> 2pi*rand(integrator.p.rng),
        nswitches::T3=0,
        initialeomid::T9=EoMID.FullOrbit,
        maxswitches::T3=100,
    ) where {
            T1<:Real,T2,T3<:Int,T4<:TerminationCode.T,T5<:AbstractRNG,T6<:Function,
            T9<:EoMID.T}
        timeatswitch = zeros(Float64, maxswitches)
        eomidafterswitch = fill(initialeomid, maxswitches)
        magneticmomentafterswitch = zeros(Float64, maxswitches)
        new{T1,T2,T3,T4,T5,T6,Vector{Float64},typeof(eomidafterswitch),T9}(
            charge,
            mass,
            electromagneticfield,
            magneticmoment,
            magneticmoment,
            weight,
            nrejections,
            terminationcode,
            rng,
            getphase,
            initialeomid,
            nswitches,
            initialeomid,
            maxswitches,
            timeatswitch,
            eomidafterswitch,
            magneticmomentafterswitch
        )
    end
end

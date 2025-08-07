"""
    callbacks.jl

This file defines various conditions and affects used as callbacks using
DifferentialEquations.jl.
"""

# -----------------------------------------------------------------------------
# Out of bounds
# -----------------------------------------------------------------------------

"""
    OutOfBoundsCondition
        bounds::NTuple{N, Tuple{Real, Real} where N
        u_idxs::AbstractVector
Returns true if the particle's statevector variable `u_idxs[i]` is outside the
bounds given by `bounds[i]`.

E.g. if the particle has the statevector `(x, y, z, vx, vy, vz)`, then setting
`bounds = ( (0.0, 0.1), (1.0, 2.0) )` and `u_idxs=[1,3]` will make the condition
return true if the particle's `x` is outside `(0.0, 1.0)` or if the particle's
`z` is outside `(1.0, 2.0)`.
"""
struct OutOfBoundsCondition{
    T1<:NTuple{N, Tuple{Real, Real}} where N,
    T2<:AbstractVector,
}
    bounds::T1
    u_idxs::T2
    function OutOfBoundsCondition(bounds, u_idxs)
        if length(u_idxs) != length(bounds)
            throw(
                ArgumentError("The number of bounds and u_idxs must be equal")
            )
        end
        new{typeof(bounds), typeof(u_idxs)}(bounds, u_idxs)
    end
end
function (self::OutOfBoundsCondition)(u, _, _)
    for (i, j) in zip(eachindex(self.bounds), self.u_idxs)
        if u[j] < self.bounds[i][1] || u[j] > self.bounds[i][2]
            return true
        end
    end
    return false
end

# -----------------------------------------------------------------------------
# Relativistic particles
# -----------------------------------------------------------------------------

"""
    RelativisticConditionGCA
Callback condition for checking if the GCA particle is relativistic. Returns
`true` if the particle's kinetic energy is a user defined fraction of the rest
energy. Perpendicular drifts other than E cross B drift are neglected.
"""
mutable struct RelativisticConditionGCA{T<:Real}
    fractionofrestenergy::T

    function RelativisticConditionGCA(; mass::T, fraction) where {T}
        return new{T}(fraction * mass * csqrd)
    end
end
function (self::RelativisticConditionGCA)(u, t, integrator)
    Rx, Ry, Rz, vparal = u[1:4]
    μ = integrator.p.magneticmoment
    mass = integrator.p.mass

    E_vec, B_vec = integrator.p.electromagneticfield(Rx, Ry, Rz, t)
    v_E = exbdrift(B_vec, E_vec)
    B = norm(B_vec)
    vperp = perpendicular_velocity(μ, mass, B)
    energy = kineticenergy(vparal, vperp, v_E, mass)
    return energy > self.fractionofrestenergy
end


"""
    set_limit!(
        condition::Union{
            RelativisticConditionGCA,
            RelativisticConditionGCA_2Dxz
        };
        mass,
        fraction,
)
Sets the energy limit for the the relativistic condition.
"""
function set_limit!(
    condition::RelativisticConditionGCA,
    mass,
    fraction,
)
    condition.energylimit = fraction * mass * csqrd
end




# -----------------------------------------------------------------------------
# Hybrid switching
# -----------------------------------------------------------------------------
"""
    switch2gca_affect!(integrator)
Callback affect which switches to the guiding centre approximation by setting
the parameter `switch` to `1` and transforming `u` to the guiding centre state.
"""
function switch2gca_affect!(integrator)
    B, E = emfieldatpos(integrator.u[1:3], integrator.p.fields)
    get_guidingcentre!(
        integrator.u, # In-place (will be changed)
        integrator.p.magneticmoment,# In-place
        integrator.u[1:3],
        integrator.u[4:6],
        B, E,
        integrator.p.charge,
        integrator.p.mass
    )
    # Set switch to full orbit case
    integrator.p.switch = 1
end


"""
    switch2fo_affect!(integrator)
Callback affect which switches to full orbit integration by setting
the parameter `switch` to `2` and transforming `u` to a full orbit state.
"""
function switch2fo_affect!(integrator)
    B, E = emfieldatpos(integrator.u[1:3], integrator.p.fields)
    phaseangle = rand(integrator.p.rng, 0.0, Float64(pi))
    get_fullorbit!(
        integrator.u,
        B, E,
        integrator.u[1:3],
        integrator.u[4],
        integrator.p.magneticmoment,
        integrator.p.charge,
        integrator.p.mass,
        phaseangle
    )
    # Set switch to full orbit case
    integrator.p.switch = 2
end


"""
    switch_affect!(integrator)
Callback affect that checks the `switch` parameter and switches the EoM, either
from full orbit integration to the guiding centre approximation or the
opposite.
"""
function switch_affect!(integrator)
    println("Switch_affect!")
    println("Switch value before: $(integrator.p.switch)")
    #    println("mu before:           $(integrator.p.magneticmoment)")
    println("u before: $(integrator.u)")
    if integrator.p.switch == 1
        switch2fo_affect!(integrator)
    else
        switch2gca_affect!(integrator)
    end
    println("Switch value after:  $(integrator.p.switch)")
    #    println("mu after:            $(integrator.p.magneticmoment)")
    println("u after:  $(integrator.u)")
end


# -----------------------------------------------------------------------------
# GCA breakdown and validity condition
# -----------------------------------------------------------------------------
"""
    GCAValidCondition
Callback condition for checking GCA assumption. Returns `true` if the ratio
between the particle Larmor radius and the characteristic length of the
magnetic field is lower than a tolerance `switchtol`.
"""
struct GCAValidCondition
    tolerance::Real
end
function (functor::GCAValidCondition)(u, t, integrator)
    ratio = scalesratio(u, t, integrator.p)
    return ratio < functor.tolerance
end


"""
    GCAInvalidCondition
Callback condition for checking GCA assumption. Returns `true` if the ratio
between the particle Larmor radius and the characteristic length of the
magnetic field is higher than a `tolerance`.
"""
struct GCAInvalidCondition
    tolerance::Real
end
function (self::GCAInvalidCondition)(u, t, integrator)
    integrator.p.switch == 1 &&
        scalesratio(u, t, integrator.p, integrator.p.switch) > self.tolerance
end



"""
    MagneticGradientCondition
Callback condition for checking GCA assumption. Returns `true` if the ratio
between the particle Larmor radius and the characteristic length of the
magnetic field is higher than a tolerance `switchtol`.
"""
mutable struct MagneticGradientCondition{T<:Real}
    tolerance::T
end
function (self::MagneticGradientCondition)(u, t, integrator)
    R = SVector(u[1], u[2], u[3]) # The gyrocentre position vector
    ratio = scalesratio(
        R,
        t,
        integrator.p.mass,
        integrator.p.charge,
        integrator.p.magneticmoment,
        integrator.p.electromagneticfield,
    )
    return ratio > self.tolerance
end

"""
    MagneticCurvatureCondition
# Fields
- `tolerance`
A callback condition that returns true if the ratio between the particle's
Larmor radius and curvature radius of the local magnetic field is larger than
the `tolerance`.
"""
struct MagneticCurvatureCondition{T<:Real}
    tolerance::T
end
function (self::MagneticCurvatureCondition)(u, t, integrator)
    R = SVector(u[1], u[2], u[3])
    ratio = magneticcurvatureratio(
        R,
        t,
        integrator.p.mass,
        integrator.p.charge,
        integrator.p.magneticmoment,
        integrator.p.electromagneticfield,
    )
    return ratio > self.tolerance
end

"""
    ParallelElectricFieldCondition
# Fields
- `tolerance`
A callback condition that returns true if the local electric field component
parallel to the magnetic field is larger than the electric field component
perpendicular to the magnetic field by a factor of `tolerance`.
"""
struct ParallelElectricFieldCondition{T<:Real}
    tolerance::T
end
function (self::ParallelElectricFieldCondition)(u, t, integrator)
    Evec, Bvec = integrator.p.electromagneticfield(u[1], u[2], u[3], t)
    b = Bvec / norm(Bvec)
    Eparal = Evec ⋅ b
    Eperp = norm(Evec - Eparal * b)
    return Eparal/Eperp > self.tolerance
end

"""
    set_tol!(
        condition,
        tol
    )
Sets the tolerance for the GCABreakDownCondition.
"""
function set_tol!(
    condition,
    tol
)
    condition.tolerance = tol
end

function hybridswitch_affect!(integrator)
    integrator.p.gca = !integrator.p.gca
end

#-------------------------------------------------------------------------------
# Termination effects
#-------------------------------------------------------------------------------
"""
An EnumX module containing termination codes. The codes are set by callback
affects that terminate the solver. E.g. the callback affect `outofboundsaffect!`
will terminate the integration and set the parameter `terminationcode` to
`TerminationCode.OutOfBounds`.
"""
@enumx TerminationCode begin
    NotTerminated
    OutOfBounds
    MagneticGradient
    MagneticCurvature
    ParallelElectricField
    Relativistic
end

"""
    outofboundsnaffect!(integrator)
Callback affect which terminates the integration and sets the termination
message to `TerminationCode.OutOfBounds`.
"""
function outofboundsaffect!(integrator)
    integrator.p.terminationcode = TerminationCode.OutOfBounds
    terminate!(integrator)
end

"""
    magneticgradientaffect!(integrator)
Callback affect which terminates the integration and sets the return message to
`TerminationCode.MagneticGradient`.
"""
function magneticgradientaffect!(integrator)
    integrator.p.terminationcode = TerminationCode.MagneticGradient
    terminate!(integrator)
end

"""
    magneticcurvatureaffect!(integrator)
Callback affect which terminates the integration and sets the return message to
`TerminationCode.MagneticCurvature`.
"""
function magneticcurvatureaffect!(integrator)
    integrator.p.terminationcode = TerminationCode.MagneticCurvature
    terminate!(integrator)
end

"""
    parallelelectricfieldaffect!(integrator)
Callback affect which terminates the integration and sets the return message to
`TerminationCode.ParallelElectricField`.
"""
function parallelelectricfieldaffect!(integrator)
    integrator.p.terminationcode = TerminationCode.ParallelElectricField
    terminate!(integrator)
end

"""
    relativisticaffect!(integrator)
Callback affect which terminates the integration and sets the return message to
`TerminationCode.Relativistic`.
"""
function relativisticaffect!(integrator)
    integrator.p.terminationcode = TerminationCode.Relativistic
    terminate!(integrator)
end

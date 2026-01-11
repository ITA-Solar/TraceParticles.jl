"""
    callbacks.jl

This file defines various conditions and affects used as callbacks using
DifferentialEquations.jl.
"""

# -----------------------------------------------------------------------------
# Out of bounds
# -----------------------------------------------------------------------------

"""
    struct OutOfBoundsCondition{
        T1<:NTuple{N,Tuple{Real,Real}} where {N},
        T2<:AbstractVector,
    }

## Fields
- `bounds::T1`
- `u_idxs::T2`

Returns true if the particle's statevector variable `u_idxs[i]` is outside the
bounds given by `bounds[i]`.

E.g. if the particle has the statevector `(x, y, z, vx, vy, vz)`, then setting
`bounds = ( (0.0, 0.1), (1.0, 2.0) )` and `u_idxs=[1,3]` will make the condition
return true if the particle's `x` is outside `(0.0, 1.0)` or if the particle's
`z` is outside `(1.0, 2.0)`.
"""
struct OutOfBoundsCondition{
    T1<:NTuple{N,Tuple{Real,Real}} where {N},
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
        new{typeof(bounds),typeof(u_idxs)}(bounds, u_idxs)
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
    Rx, Ry, Rz, vparal = u[1], u[2], u[3], u[4]
    μ = integrator.p.magneticmoment
    mass = integrator.p.mass

    E_vec, B_vec = integrator.p.electromagneticfield(Rx, Ry, Rz, t)
    v_E = exbdrift(B_vec, E_vec)
    B = norm(B_vec)
    vperp = perpendicular_velocity(μ, mass, B)
    energy = kineticenergy(vparal, vperp, v_E, mass)
    return energy > self.fractionofrestenergy
end

# -----------------------------------------------------------------------------
# GCA breakdown and validity conditions
# -----------------------------------------------------------------------------
"""
    HybridSwitchCondition
Callback condition for checking whether it's time to switch equations of motion
when using the hybrid GCA/FO method. The condition returns `true` if the current
equations are full orbit and the GCA is determined to be valid, or if the
current equations are GCA but the GCA is determined to be invalid.

The GCA is determined to be invalid in the following cases:
- The ratio between the particle's Larmor radius and the characteristic length
  of the magnetic field strength is larger than `tol_gradient`.
- The ratio between the particle's Larmor radius and the curvature radius of
  the magnetic field is larger than `tol_curvature`.
- The ratio between the parallel and perpendicular electric field components
  is larger than `tol_efield`.

The GCA is determined to be valid if the points above are false even with an
additional factor `switchback_tol` applied to the tolerances. The switchback
tolerance is used to avoid rapid oscillations between the GCA and full orbit.
"""
struct HybridSwitchCondition{T<:Real}
    tol_gradient::T
    tol_curvature::T
    tol_efield::T
    switchback_tol::T # Set to 0.9 by Adam Stanier
end
function (self::HybridSwitchCondition)(u, t, integrator)
    x, y, z = u[1], u[2], u[3]
    pos = SVector(x, y, z)
    q, m = integrator.p.charge, integrator.p.mass
    Eratio = eratio(pos, t, integrator.p.electromagneticfield)
    if integrator.p.eomid == 1
        # Current scheme is guiding centre approximation
        kratio = magneticcurvatureratio(
            pos, t, m, q, integrator.p.magneticmoment,
            integrator.p.electromagneticfield
        )
        gratio = scalesratio(
            pos, t, m, q, integrator.p.magneticmoment,
            integrator.p.electromagneticfield
        )

        c1 = abs(gratio) > self.tol_gradient
        c2 = abs(kratio) > self.tol_curvature
        c3 = abs(Eratio) > self.tol_efield
        if c1 || c2 || c3
            return true
        else
            return false
        end
    elseif integrator.p.eomid == 2
        vel = SVector(u[4], u[5], u[6]) # The full orbit velocity vector
        kratio = magneticcurvatureratio(
            pos, vel, t, m, q,
            integrator.p.electromagneticfield
        )
        gratio = scalesratio(
            pos, vel, t, m, q,
            integrator.p.electromagneticfield
        )
        c1 = abs(gratio) < self.tol_gradient * self.switchback_tol
        c2 = abs(kratio) < self.tol_curvature * self.switchback_tol
        c3 = abs(Eratio) < self.tol_efield * self.switchback_tol
        if c1 && c2 && c3
            return true
        else
            return false
        end
    end
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
    Eperp = eperp(Evec, b, Eparal)
    return Eparal / Eperp > self.tolerance
end

#-------------------------------------------------------------------------------
# Callback affects
#-------------------------------------------------------------------------------
"""
    hybridswitchaffect!(integrator)
Callback affect that checks the `switch` parameter and switches the EoM, either
from full orbit integration to the guiding centre approximation or the
opposite.
"""
function hybridswitchaffect!(integrator)
    if integrator.p.eomid == 1
        switch2fo_affect!(integrator)
    elseif integrator.p.eomid == 2
        switch2gca_affect!(integrator)
    else
        @warn """Trying to switch EoM, but switch parameter is neither 1 nor 2.
Current value: $(integrator.p.eomid)."""
    end
    integrator.p.nswitches += 1
end

"""
    hybridswitchaffect_withdetection!(integrator)
Callback affect that checks the `switch` parameter and switches the EoM, either
from full orbit integration to the guiding centre approximation or the
opposite. Also records the time of the switch and the new EoM ID in the
`switchtimes` and `neweomid` arrays, respectively.
"""
function hybridswitchaffect_withdetection!(integrator)
    integrator.p.nswitches += 1
    integrator.p.timeatswitch[integrator.p.nswitches] = integrator.t
    if integrator.p.eomid == 1
        integrator.p.eomidafterswitch[integrator.p.nswitches] = 2
        switch2fo_affect!(integrator)
        integrator.p.magneticmomentafterswitch[integrator.p.nswitches] = NaN
    elseif integrator.p.eomid == 2
        integrator.p.eomidafterswitch[integrator.p.nswitches] = 1
        switch2gca_affect!(integrator)
        integrator.p.magneticmomentafterswitch[integrator.p.nswitches] =
            integrator.p.magneticmoment
    else
        @warn """Trying to switch EoM, but switch paramer is neither 1 nor 2.
Current value: $(integrator.p.eomid)."""
    end
    if integrator.p.nswitches == integrator.p.maxswitches
        @warn "Maximum number of switches reached: $(integrator.p.maxswitches)."
        integrator.p.terminationcode = TerminationCode.MaxSwitches
        terminate!(integrator)
    end
end

"""
    switch2gca_affect!(integrator)
Callback affect which switches to the guiding centre approximation by setting
the parameter `switch` to `1` and transforming `u` to the guiding centre state.
"""
function switch2gca_affect!(integrator)
    integrator.p.eomid = 1
    integrator.p.magneticmoment = get_guidingcentre!(
        integrator.u,
        integrator.t,
        integrator.p.electromagneticfield,
        integrator.p.charge,
        integrator.p.mass
    )
end


"""
    switch2fo_affect!(integrator)
Callback affect which switches to full orbit integration by setting
the parameter `switch` to `2` and transforming `u` to a full orbit state.
"""
function switch2fo_affect!(integrator)
    integrator.p.eomid = 2
    phaseangle = integrator.p.getphase(integrator)
    get_fullorbit!(
        integrator.u,
        integrator.t,
        integrator.p.electromagneticfield,
        integrator.p.magneticmoment,
        integrator.p.charge,
        integrator.p.mass,
        phaseangle
    )
end

#-------------------------------------------------------------------------------
# Callback affects that terminate
#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
# Other functions
#-------------------------------------------------------------------------------
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

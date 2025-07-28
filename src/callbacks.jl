"""
    callbacks.jl

This file defines various conditions and affects used as callbacks using
DifferentialEquations.jl.
"""

# -----------------------------------------------------------------------------
# Out of bounds
# -----------------------------------------------------------------------------

"""
    OutOfDomainCondition_3D
        xbounds::Tuple{<:Real, <:Real}
        ybounds::Tuple{<:Real, <:Real}
        zbounds::Tuple{<:Real, <:Real}
Returns true if the particle is outside the bounds.
"""
mutable struct OutOfDomainCondition_3D
    xbounds::Tuple{<:Real,<:Real}
    ybounds::Tuple{<:Real,<:Real}
    zbounds::Tuple{<:Real,<:Real}
end
function (self::OutOfDomainCondition_3D)(u, _, _)
    return u[1] <= self.xbounds[1] || u[1] >= self.xbounds[2] ||
           u[2] <= self.ybounds[1] || u[2] >= self.ybounds[2] ||
           u[3] <= self.zbounds[1] || u[3] >= self.zbounds[2]
end


"""
    OutOfDomainCondition_2Dxz
        xbounds::Tuple{<:Real, <:Real}
        zbounds::Tuple{<:Real, <:Real}
Returns true if the particle is outside the 2D bounds in x and z.
"""
mutable struct OutOfDomainCondition_2Dxz
    xbounds::Tuple{<:Real,<:Real}
    zbounds::Tuple{<:Real,<:Real}
end
function (self::OutOfDomainCondition_2Dxz)(u, _, _)
    return u[1] <= self.xbounds[1] || u[1] >= self.xbounds[2] ||
           u[3] <= self.zbounds[1] || u[3] >= self.zbounds[2]
end


"""
    outside2dnullpointzoom(u, t, integrator)
Callback condition for particles exiting the zoomed-in domain in the 2D
nullpoint Bifrost simulation used for the results in the Japan Hinode/IRIS
poster of septermber 2023.

The effect of this callback should be termination, which means that the
callback should be created like this:
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition,affect!)
"""
function outside2dnullpointzoom(u, _, _) # (u, t, integrator)
    return u[1] <= 14.52e6 || u[1] >= 18.83e6 ||
           u[3] <= -7.79e6 || u[3] >= -3.88e6
end


"""
    outofdomainaffect!(integrator)
Callback affect which terminates the integration and sets the return message to
"OutofDomain".
"""
function outofdomainaffect!(integrator)
    integrator.p.userdata.retmsg = "OutofDomain"
    terminate!(integrator)
end


# -----------------------------------------------------------------------------
# Relativistic particles
# -----------------------------------------------------------------------------

"""
    RelativisticConditionGCA
Callback condition for checking if the GCA particle is relativistic. Returns
`true` if the particle's kinetic energy is a user defined fraction of the rest
energy. Perpendicular drifts other than E cross B drift are neglected.

Default fraction is 0.002, which is 1022 eV for an electron.
"""
mutable struct RelativisticConditionGCA
    fractionofrestenergy::Real

    function RelativisticConditionGCA(; mass, fraction)
        return new(fraction * mass * csqrd)
    end
end
function (self::RelativisticConditionGCA)(u, _, integrator)
    Rx, Ry, Rz, vparal = u[1:4]
    μ = integrator.p.magneticmoment
    mass = integrator.p.mass

    E_vec, B_vec = integrator.p.electromagneticfield(Rx, Ry, Rz)
    v_E = exbdrift(B_vec, E_vec)
    B = norm(B_vec)
    vperp = perpendicular_velocity(μ, mass, B)
    energy = kineticenergy(vparal, vperp, v_E, mass)
    return energy > self.fractionofrestenergy
end


"""
    RelativisticConditionGCA_2Dxz
Callback condition for checking if the GCA particle is relativistic. Returns
`true` if the particle's kinetic energy is a user defined fraction of the rest
energy. Perpendicular drifts other than E cross B drift are neglected.

Default fraction is 0.002, which is 1022 eV for an electron.
"""
mutable struct RelativisticConditionGCA_2Dxz
    energylimit::Real

    function RelativisticConditionGCA_2Dxz(; mass, fraction)
        return new(fraction * mass * csqrd)
    end
end
function (self::RelativisticConditionGCA_2Dxz)(u, _, integrator)
    Rx, _, Rz, vparal = u[1:4]
    μ = integrator.p.magneticmoment
    mass = integrator.p.mass

    E_vec, B_vec = integrator.p.electromagneticfield(Rx, Rz)
    v_E = exbdrift(B_vec, E_vec)
    B = norm(B_vec)
    vperp = perpendicular_velocity(μ, mass, B)
    energy = kineticenergy(vparal, vperp, v_E, mass)
    return energy > self.energylimit
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
    condition::Union{RelativisticConditionGCA,RelativisticConditionGCA_2Dxz};
    mass,
    fraction,
)
    condition.energylimit = fraction * mass * csqrd
end


"""
    relativisticaffect!(integrator)
Callback affect which terminates the integration and sets the return message to
"Relativistic".
"""
function relativisticaffect!(integrator)
    integrator.p.userdata.retmsg = "Relativistic"
    terminate!(integrator)
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
    GCABreakDownCondition
Callback condition for checking GCA assumption. Returns `true` if the ratio
between the particle Larmor radius and the characteristic length of the
magnetic field is higher than a tolerance `switchtol`.
"""
mutable struct GCABreakDownCondition
    tolerance::Real
end
function (self::GCABreakDownCondition)(u, t, integrator)
    R = SVector(u[1], u[2], u[3]) # The gyrocentre position vector
    ratio = scalesratio(
        R,
        integrator.p.mass,
        integrator.p.charge,
        integrator.p.magneticmoment,
        integrator.p.electromagneticfield,
    )
    return ratio > self.tolerance
end


"""
    GCABreakDownCondition_2Dxz
Callback condition for checking GCA assumption. Returns `true` if the ratio
between the particle Larmor radius and the characteristic length of the
magnetic field is higher than a tolerance `switchtol`.
"""
mutable struct GCABreakDownCondition_2Dxz
    tolerance::Real
end
function (self::GCABreakDownCondition_2Dxz)(u, t, integrator)
    R = SVector(u[1], u[2], u[3]) # The gyrocentre position vector in 2D
    ratio = scalesratio(
        R,
        integrator.p.mass,
        integrator.p.charge,
        integrator.p.magneticmoment,
        integrator.p.electromagneticfield,
    )
    return ratio > self.tolerance
end


"""
    set_tol!(
        condition::Union{GCABreakDownCondition, GCABreakDownCondition_2Dxz},
        tol
    )
Sets the tolerance for the GCABreakDownCondition.
"""
function set_tol!(
    condition::Union{GCABreakDownCondition,GCABreakDownCondition_2Dxz},
    tol
)
    condition.tolerance = tol
end


"""
    gcabreakdownaffect!(integrator)
Callback affect which terminates the integration and sets the return message to
"GCABreakDown".
"""
function gcabreakdownaffect!(integrator)
    integrator.p.userdata.retmsg = "GCABreakDown"
    terminate!(integrator)
end

function hybridswitch_affect!(integrator)
    integrator.p.gca = !integrator.p.gca
end


#------------------------------------------------------------------------------
# MAX callback for LARMOR RADIUS
# Condition
function maxcondition_larmorradius_2Dxz(u, _, integrator)
    R = [u[1], u[3]]
    r_L = larmorradius(R, integrator.p)
    return r_L > integrator.p.larmorradius.max
end
# Affect
function maxaffect_larmorradius_2Dxz!(integrator)
    R = [integrator.u[1], integrator.u[3]]
    r_L = larmorradius(R, integrator.p)
    integrator.p.larmorradius.max = r_L
    integrator.p.larmorradius.max_u = copy(integrator.u)
    integrator.p.larmorradius.max_t = integrator.t
end

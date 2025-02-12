"""
    callbacks.jl

This file defines various conditions and affects used as callbacks using
DifferentialEquations.jl.
"""

# -----------------------------------------------------------------------------
# Out of bounds
# -----------------------------------------------------------------------------

"""
    OutOfDomainCondition_2Dxz
        xbounds::Tuple{<:Real, <:Real}
        zbounds::Tuple{<:Real, <:Real}
Returns true if the particle is outside the 2D bounds in x and z.
"""
struct OutOfDomainCondition_2Dxz
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
function outofdomainaffect!(integrator)
    integrator.p.userdata.retmsg = "OutofDomain"
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
    GCABreakDownCondition
Callback condition for checking GCA assumption. Returns `true` if the ratio
between the particle Larmor radius and the characteristic length of the
magnetic field is higher than a tolerance `switchtol`.
"""
struct GCABreakDownCondition
    tolerance::Real
end
function (functor::GCABreakDownCondition)(u, t, integrator)
    ratio = scalesratio(u, t, integrator.p)
    return ratio > functor.tolerance
end


struct GCABreakDownCondition_2Dxz
    tolerance::Real
end
function (functor::GCABreakDownCondition_2Dxz)(u, t, integrator)
    ratio = scalesratio_2Dxz(u, t, integrator.p)
    return ratio > functor.tolerance
end

struct GCABreakDownCondition_highmem_2Dxz
    tolerance::Real
end
function (functor::GCABreakDownCondition_highmem_2Dxz)(u, t, integrator)
    ratio = scalesratio_highmem_2Dxz(u, t, integrator.p)
    return ratio > functor.tolerance
end


function GCABreakDownCB(
    method::String,
    dimensionality::String
    ;
    tolerance::Real=1.0
)
    if method * dimensionality == "lowmem2Dxz"
        condition = GCABreakDownCondition_2Dxz(tolerance)
    elseif method * dimensionality == "highmem2Dxz"
        condition = GCABreakDownCondition_highmem_2Dxz(tolerance)
    else
        error("Unknown method.")
    end
    return DiscreteCallback(condition, gcabreakdownaffect!)
end

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

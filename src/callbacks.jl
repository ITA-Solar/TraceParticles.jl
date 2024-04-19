
abstract type AbstractTPCallback end

struct DiscreteTPCallback{F1, F2} <: AbstractTPCallback
    condition::F1
    affect!  ::F2
    function DiscreteTPCallback(condition::F1, affect!::F2) where {F1, F2}
        new{F1, F2}(condition, affect!)
    end
end
function DiscreteTPCallback()
    condition(x) = false
    affect! = nothing
    return DiscreteTPCallback(condition, affect!)
end

# -----------------------------------------------------------------------------
# Out of bounds
# -----------------------------------------------------------------------------
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
function outside2dnullpointzoom(u,_,_) # (u, t, integrator)
    return u[1] <= 14.52e6 || u[1] >= 18.83e6 ||
            u[3] <= -7.79e6 || u[3] >= -3.88e6
end

function killparticle!(integrator)
    terminate!(integrator)
end


# -----------------------------------------------------------------------------
# Hybrid switching
# -----------------------------------------------------------------------------
function gcabreakdown_highmem2Dxz(u, t, integrator)
    ratio = scalesratio_highmem_2Dxz(u, t, integrator.p)
    return ratio > integrator.p.switchtol
end
function gcabreakdown_2Dxz(u, t, integrator)
    ratio = scalesratio_2Dxz(u, t, integrator.p)
    return ratio > integrator.p.switchtol
end
function scalesratio(
    R::Vector{<:Real},
    itpvec::Vector{<:AbstractInterpolation},
    ∇B::Vector{<:Real},
    params::NamedTuple
    )
    μ = params.magneticmoment
    q = params.charge
    m = params.mass
    B = norm([itp(R...) for itp in itpvec])
    L_B = characteristicfieldlength(B, ∇B)
    vperp = perpendicular_velocity(μ, m, B)
    r_L = larmorradius(m, vperp, q, B)
    return r_L/L_B
end
function scalesratio_2Dxz(u, t, params)
    R = [u[1], u[3]]
    itpvec = params.fields[1:3]
    ∇B = ∇(R, itpvec)
    scalesratio(R, itpvec, ∇B, params)
end
function scalesratio_highmem_2Dxz(u, t, params)
    R = [u[1], u[3]]
    itpvec = params.fields
    ∇B = [itpvec[i](R) for i in 7:9]
    scalesratio(R, itpvec, ∇B, params)
end


function hybridswitch_affect!(integrator)
    integrator.p.gca = !integrator.p.gca
end


#____/\_____/\_________________________________________________________________
# Save maximum values
mutable struct MaxValue{T1, T2, T3}
    max::T1
    max_u::T2
    max_t::T3
end
mutable struct MinValue{T1, T2, T3}
    min::T1
    min_u::T2
    min_t::T3
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
# Helper function
function larmorradius(
    R::Vector{<:Real},
    params::NamedTuple
    )
    itpvec = params.fields[1:3]
    μ = params.magneticmoment
    q = params.charge
    m = params.mass
    r_L = larmorradius(R, μ, q, m, itpvec)
end

#------------------------------------------------------------------------------
# MIN callback for characteristic field length of the magnetic field.
function mincondition_fieldlength_2Dxz(u, _, integrator)
    itpvec = integrator.p.fields[1:3]
    L_B = characteristicfieldlength([u[1], u[3]], itpvec)
    return L_B < integrator.p.fieldlength.min
end
function mincondition_fieldlength_highmemory_2Dxz(u, _, integrator)
    itpvec = integrator.p.fields
    R = [u[1],u[3]]
    L_B = characteristicfieldlength(R, itpvec[1:3], itpvec[7:9])
    return L_B < integrator.p.fieldlength.min
end

# Affects
function minaffect_fieldlength_2Dxz!(integrator)
    itpvec = integrator.p[:fields][1:3]
    L_B = characteristicfieldlength([integrator.u[1], integrator.u[3]], itpvec)
    integrator.p.fieldlength.min = L_B
    integrator.p.fieldlength.min_u = copy(integrator.u)
    integrator.p.fieldlength.min_t = integrator.t
end
function minaffect_fieldlength_highmemory_2Dxz!(integrator)
    itpvec = integrator.p.fields
    R = [integrator.u[1],integrator.u[3]]
    L_B = characteristicfieldlength(R, itpvec[1:3], itpvec[7:9])
    integrator.p.fieldlength.min = L_B
    integrator.p.fieldlength.min_u = copy(integrator.u)
    integrator.p.fieldlength.min_t = integrator.t
end

#------------------------------------------------------------------------------
# MAX callback for SCALESRATIO
function maxcondition_scalesratio_2Dxz(u, t, integrator)
    ratio = scalesratio_2Dxz(u, t, integrator.p)
    return ratio > integrator.p.scalesratio.max
end
function maxcondition_scalesratio_highmemory_2Dxz(u, t, integrator)
    ratio = scalesratio_highmem_2Dxz(u, t, integrator.p)
    return ratio > integrator.p.scalesratio.max
end

function maxaffect_scalesratio_2Dxz!(integrator)
    ratio = scalesratio_2Dxz(integrator.u, integrator.t, integrator.p)
    integrator.p.scalesratio.max = ratio
    integrator.p.scalesratio.max_u = copy(integrator.u)
    integrator.p.scalesratio.max_t = integrator.t
end
function maxaffect_scalesratio_highmemory_2Dxz!(integrator)
    ratio = scalesratio_highmem_2Dxz(integrator.u, integrator.t, integrator.p)
    integrator.p.scalesratio.max = ratio
    integrator.p.scalesratio.max_u = copy(integrator.u)
    integrator.p.scalesratio.max_t = integrator.t
end

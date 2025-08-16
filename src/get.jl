#
# get.jl
#
# Created 11.02.25
#
# For fetching data and calculating observables from ODE-solutions and 
# ensembles stored as DataFrames.
#
#____/\_____/\_________________________________________________________________

"""
    getinitialstate(df::DataFrame)
Return the initial state  of particle simulation stored as a DataFrame.
"""
function getinitialstate(df::DataFrame)
    return df[:, [:x0, :y0, :z0, :vparal0]]
end


"""
    getfinalstate(df::DataFrame)
Return the final state of particle simulation stored as a DataFrame.
"""
function getfinalstate(df::DataFrame)
    return df[:, [:xf, :yf, :zf, :vparalf]]
end


"""
    getfinaltimes(df::DataFrame)
Return the start time of the particle simulation stored as a DataFrame.
"""
function getinitialtimes(df::DataFrame)
    return df.t0
end


"""
    getfinaltimes(df::DataFrame)
Return the end time of the particle simulation stored as a DataFrame.
"""
function getfinaltimes(df::DataFrame)
    return df.tf
end


"""
    getmagneticmoments(df::DataFrame)
Return the magnetic moments of the particle simulation stored as a DataFrame.
"""
function getmagneticmoments(df::DataFrame)
    return df.magneticmoment
end


"""
Get the observable `sym` from an `ODESolution` or a vector of `ODESolution`s at
a single time or over a range of times.

# Methods
    get_observable(vectorofsolutions, symbol; t, kwargs...)
    get_observable(solution, symbol, t; kwargs...)
    get_observable(solution, symbol; tspan, times, nsteps, kwargs...)
"""
function get_observable(
    sol::Union{EnsembleSolution,VectorOfArray},
    sym::Symbol;
    t=nothing,
    kwargs...
)
    if sym == :t0
        return [last(u.t) for u in sol]
    elseif sym == :tf
        return [last(u.t) for u in sol]
    elseif sym == :x0
        return [u[1, 1] for u in sol]
    elseif sym == :xf
        return [u[1, end] for u in sol]
    elseif sym == :x0
        return [u[3, 1] for u in sol]
    elseif sym == :zf
        return [u[3, end] for u in sol]
    elseif sym == :e0
        return [get_energy(u, first(u.t); kwargs...) for u in sol]
    elseif sym == :ef
        return [get_energy(u, last(u.t); kwargs...) for u in sol]
    elseif sym == :nt
        return [length(u.t) for u in sol]
    elseif sym == :terminationcode
        return [u.prob.p.userdata.terminationcode for u in sol]
    elseif isnothing(t)
        throw(ArgumentError(
            "Observable $sym is not implemented without a specified time." *
            "Try setting the keyword argument `t`, or maybe you want to pass" *
            " a single `ODESolution` instead of several?`"
        ))
    end
    # Variables at specific times
    if sym == :x
        return [u(t)[1] for u in sol]
    elseif sym == :y
        return [u(t)[2] for u in sol]
    elseif sym == :z
        return [u(t)[3] for u in sol]
    elseif sym == :vparal
        return [u(t)[4] for u in sol]
    else
        error("Observable $sym is not implemented.")
    end
end


function get_observable(
    sol::ODESolution,
    sym::Symbol;
    tspan=(first(sol.t), last(sol.t)),
    nsteps=1000,
    times=range(tspan[1], tspan[2], length=nsteps),
    magnitude=false,
    component=nothing,
    EoM="GCA",
)

    if sym == :x
        obs = [sol(t)[1] for t in times]
    elseif sym == :y
        obs = [sol(t)[2] for t in times]
    elseif sym == :z
        obs = [sol(t)[3] for t in times]
    elseif sym == :vparal
        obs = [sol(t)[4] for t in times]
    elseif sym == :vx
        obs = [sol(t)[4] for t in times]
    elseif sym == :vy
        obs = [sol(t)[5] for t in times]
    elseif sym == :vz
        obs = [sol(t)[6] for t in times]
    elseif sym == :t
        obs = times
    elseif sym == :mu || sym == :magneticmoment
        obs = [
            TraceParticles.magneticmoment(
                sol(t, idxs=1:3),
                sol(t, idxs=4:6),
                t,
                sol.prob.p.mass,
                sol.prob.p.electromagneticfield
            ) for t in times
        ]
    elseif sym == :vperp
        if EoM == "GCA"
            obs = [
                TraceParticles.perpendicular_velocity(
                    sol(t, idxs=1:3),
                    t,
                    sol.prob.p.magneticmoment,
                    sol.prob.p.mass,
                    sol.prob.p.electromagneticfield
                ) for t in times
            ]
        elseif EoM == "FO"
            obs = [
                TraceParticles.perpendicular_velocity(
                    sol(t, idxs=1:3),
                    sol(t, idxs=4:6),
                    t,
                    sol.prob.p.electromagneticfield,
                ) for t in times
            ]
        end
    elseif sym == :r_L
        if EoM == "GCA"
            obs = [
                larmorradius(
                    sol(t, idxs=1:3),
                    t,
                    sol.prob.p.magneticmoment,
                    sol.prob.p.charge,
                    sol.prob.p.mass,
                    sol.prob.p.electromagneticfield
                ) for t in times
            ]
        elseif EoM == "FO"
            obs = [
                larmorradius(
                    sol(t, idxs=1:3),
                    sol(t, idxs=4:6),
                    t,
                    sol.prob.p.mass,
                    sol.prob.p.charge,
                    sol.prob.p.electromagneticfield
                ) for t in times
            ]
        end
    elseif sym == :L_B
        obs = [
            characteristicfieldlength(
                sol(t, idxs=1:3), t, sol.prob.p.electromagneticfield
            ) for t in times
        ]
    elseif sym == :scalesratio
        if EoM == "GCA"
            obs = [
                scalesratio(
                    sol(t, idxs=1:3),
                    t,
                    sol.prob.p.mass,
                    sol.prob.p.charge,
                    sol.prob.p.magneticmoment,
                    sol.prob.p.electromagneticfield,
                ) for t in times
            ]
        elseif EoM == "FO"
            obs = [
                scalesratio(
                    sol(t, idxs=1:3),
                    sol(t, idxs=4:6),
                    t,
                    sol.prob.p.mass,
                    sol.prob.p.charge,
                    sol.prob.p.electromagneticfield,
                ) for t in times
            ]
        end
    elseif sym == :exbdrift
        obs = [
            exbdrift(
                sol(t, idxs=1:3), t, sol.prob.p.electromagneticfield
            ) for t in times
        ]
    elseif sym == :magneticcurvatureratio || sym == :mcr
        if EoM == "GCA"
            obs = [
                magneticcurvatureratio(
                    sol(t, idxs=1:3),
                    t,
                    sol.prob.p.mass,
                    sol.prob.p.charge,
                    sol.prob.p.magneticmoment,
                    sol.prob.p.electromagneticfield
                ) for t in times
            ]
        elseif EoM == "FO"
            obs = [
                magneticcurvatureratio(
                    sol(t, idxs=1:3),
                    sol(t, idxs=4:6),
                    t,
                    sol.prob.p.mass,
                    sol.prob.p.charge,
                    sol.prob.p.electromagneticfield
                ) for t in times
            ]
        end
    elseif sym == :gyrofrequency
        obs = [norm(sol.prob.p.electromagneticfield(sol(t)[1:3]..., t)[2]) *
               sol.prob.p.charge / sol.prob.p.mass for t in times]
    elseif sym == :gyroperiod
        factor = 2π * sol.prob.p.mass /
                 sol.prob.p.charge
        obs = factor / [
            norm(sol.prob.p.electromagneticfield(sol(t)[1:3]..., t)[2])
            for t in times
        ]
    elseif sym == :pitchangle
        if EoM == "GCA"
            obs = [
                pitchangle(
                    sol.prob.p.electromagneticfield(sol(t, idxs=1:3)..., t)[2],
                    sol(t, idxs=4),
                    sol.prob.p.mass,
                    sol.prob.p.magneticmoment,
                ) for t in times
            ]
        elseif EoM == "FO"
            B = [
                sol.prob.p.electromagneticfield(sol(t, idxs=1:3)..., t)[2]
                for t in times
            ]
            vel = [sol(t, idxs=4:6) for t in times]
            obs = pitchangle.(vel, B)
        end
    elseif sym == :energy
        obs = [get_energy(sol, t; EoM=EoM) for t in times]
    elseif sym == :fermi
        obs = [get_fermi(sol, t) for t in times]
    elseif sym == :betatron
        obs = [get_betatron(sol, t) for t in times]
    elseif sym == :polarisationacc
        obs = [get_polarisationacc(sol, t) for t in times]
    elseif sym == :parallelpower
        obs = [get_parallelpower(sol, t) for t in times]
    elseif sym == :driftenergy
        obs = [get_driftenergy(sol, t) for t in times]
    elseif sym == :parallelenergy
        obs = [get_parallelenergy(sol, t) for t in times]
    elseif sym == :perpenergy
        obs = [get_perpenergy(sol, t) for t in times]
    elseif sym == :eparal
        obs = [get_eparal(sol, t) for t in times]
    elseif sym == :eperp
        obs = [get_eperp(sol, t) for t in times]
    elseif sym == :efieldratio
        eparal = [get_eparal(sol, t) for t in times]
        eperp = [norm(get_eperp(sol, t)) for t in times]
        obs = eparal ./ eperp
    elseif sym == :efield
        obs = [
            sol.prob.p.electromagneticfield(sol(t)[1:3]..., t)[1]
            for t in times
        ]
    elseif sym == :bfield
        obs = [
            sol.prob.p.electromagneticfield(sol(t)[1:3]..., t)[2]
            for t in times
        ]
    elseif sym == :b
        bfield = [
            sol.prob.p.electromagneticfield(sol(t)[1:3]..., t)[2]
            for t in times
        ]
        obs = [b / norm(b) for b in bfield]
    elseif sym in fieldnames(GCAState)
        charge = sol.prob.p.charge
        mass = sol.prob.p.mass
        magneticmoment = sol.prob.p.magneticmoment
        emfield = sol.prob.p.electromagneticfield
        obs = [
            getfield(
                GCAState(
                    sol(t)[1:4],
                    t,
                    charge,
                    mass,
                    magneticmoment,
                    emfield
                ),
                sym)
            for t in times
        ]
    else
        error("Observable $sym is not implemented.")
    end
    if !isnothing(component)
        obs = [o[component] for o in obs]
    end
    if magnitude
        obs = norm.(obs)
    end
    return obs
end


function get_observable(
    sol::ODESolution,
    observable::Symbol,
    t::Real
    ;
    kwargs...
)
    return get_observable(sol, observable; times=[t], nsteps=1, kwargs...)[1]
end

"""
    get_fermi(
        sol::ODESolution,
        t::Real,
    )
Calculate the Fermi acceleration at time `t` of `sol`, in eV. That is the
energy change due to 1st order Fermi acceleration at time `t`.
"""
function get_fermi(
    sol::ODESolution,
    t::Real;
)
    q = sol.prob.p.charge
    gcastate = GCAState(
        sol(t, idxs=1:4),
        t,
        q,
        sol.prob.p.mass,
        sol.prob.p.magneticmoment,
        sol.prob.p.electromagneticfield;
    )
    v_C = gcastate.curvaturedrift
    E = gcastate.efield
    return v_C ⋅ E * q * J2eV
end


"""
    get_betatron(
        sol::ODESolution,
        t::Real,
    )
Calculate the betatron acceleration at time `t` of `sol`, in eV. That is the
energy change due to betatron acceleration at time `t`.
"""
function get_betatron(
    sol::ODESolution,
    t::Real;
)
    q = sol.prob.p.charge
    gcastate = GCAState(
        sol(t, idxs=1:4),
        t,
        q,
        sol.prob.p.mass,
        sol.prob.p.magneticmoment,
        sol.prob.p.electromagneticfield;
    )
    v_gradB = gcastate.∇Bdrift
    E = gcastate.efield
    return v_gradB ⋅ E * q * J2eV
end


"""
    get_parallelpower(
        sol::ODESolution,
        t::Real;
    )
Calculate the energy change per time due to parallel electric fields_itp.
"""
function get_parallelpower(
    sol::ODESolution,
    t::Real;
)
    q = sol.prob.p.charge
    m = sol.prob.p.mass
    gcastate = GCAState(
        sol(t, idxs=1:4),
        q,
        t,
        sol.prob.p.mass,
        sol.prob.p.magneticmoment,
        sol.prob.p.electromagneticfield;
    )
    vparal = gcastate.state[4]
    paralacc = gcastate.paralacc
    return vparal * paralacc * m * J2eV
end


"""
    get_polarisationacc(
        sol::ODESolution,
        t::Real,
    )
Calculate the polarisation acceleration at time `t` of `sol`, in eV.
That is the energy change due to polarisation acceleration at time `t`.
"""
function get_polarisationacc(
    sol::ODESolution,
    t::Real;
)
    q = sol.prob.p.charge
    gcastate = GCAState(
        sol(t, idxs=1:4),
        t,
        q,
        sol.prob.p.mass,
        sol.prob.p.magneticmoment,
        sol.prob.p.electromagneticfield;
    )
    v_P = gcastate.polarisationdrift
    E = gcastate.efield
    return v_P ⋅ E * q * J2eV
end


"""
    get_energy(
        sol::ODESolution,
        t::Real
        ;
        units="eV",
        switch=sol.prob.p.switch,
    )
Calculate the kinetic energy of the `ODESolution` `sol` at time `t`.
"""
function get_energy(
    sol::ODESolution,
    t::Real
    ;
    units="eV",
    EoM="GCA",
)
    if EoM == "FO" # Full orbit
        energy = kineticenergy(
            sol(t, idxs=4:6),
            sol.prob.p.mass,
        )
    elseif EoM == "GCA" # Guiding centre approximation
        energy = kineticenergy(GCAState(
            sol(t, idxs=1:4),
            t,
            sol.prob.p.charge,
            sol.prob.p.mass,
            sol.prob.p.magneticmoment,
            sol.prob.p.electromagneticfield,
        ))
    end
    if units == "eV"
        return energy * J2eV
    else
        return energy
    end
end

"""
    get_driftenergy(
        sol::ODESolution,
        t::Real,
    )
Calculate the kinetic energy of the `ODESolution` `sol` at time `t` from
perpendicular drifts.
"""
function get_driftenergy(
    sol::ODESolution,
    t::Real;
)
    q = sol.prob.p.charge
    m = sol.prob.p.mass
    gcastate = GCAState(
        sol(t)[1:4],
        t,
        q,
        m,
        sol.prob.p.magneticmoment,
        sol.prob.p.electromagneticfield;
    )
    v_E = gcastate.exbdrift
    v_gradB = gcastate.∇Bdrift
    v_C = gcastate.curvaturedrift
    v_P = gcastate.polarisationdrift
    return norm(v_E .+ v_gradB .+ v_C .+ v_P)^2 * 0.5 * m * J2eV
end


"""
    get_parallelenergy(
        sol::ODESolution,
        t::Real;
    )
Calculate the parallel kinetic energy of the `ODESolution` `sol` at time `t`.

"""
function get_parallelenergy(
    sol::ODESolution,
    t::Real;
)
    return 0.5 * sol(t)[4]^2 * sol.prob.p.mass * J2eV
end


"""
    get_perpenergy(
        sol::ODESolution,
        t::Real;
    )
Calculate the perpendicular kinetic energy of the `ODESolution` `sol` at time `t`.
"""
function get_perpenergy(
    sol::ODESolution,
    t::Real;
)
    _, B = sol.prob.p.electromagneticfield(sol(t)[1:3]...)
    return norm(B) * sol.prob.p.magneticmoment * J2eV
end

"""
    get_eparal(
        sol::ODESolution,
        t::Real;
    )
Calculate the parallel velocity of the `ODESolution` `sol` at time `t`.
"""
function get_eparal(
    sol::ODESolution,
    t::Real;
)
    efield, bfield = sol.prob.p.electromagneticfield(sol(t, idxs=1:3)..., t)
    B = norm(bfield)
    b = bfield / B
    return efield ⋅ b
end

function get_eperp(
    sol::ODESolution,
    t::Real;
)
    efield, bfield = sol.prob.p.electromagneticfield(sol(t, idxs=1:3)..., t)
    B = norm(bfield)
    b = bfield / B
    return efield - (efield ⋅ b) * b
end

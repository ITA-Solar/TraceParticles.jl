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
Get the observable `sym` from an `ODESolution` or a vector of `ODESolution`s.

# Methods
    get_observable(vectorofsolutions, symbol; t=nothing, kwargs...)
    get_observable(
        solution, symbol;
        times=nothing, 
        tspan=(first(sol.t), last(sol.t)), 
        nsteps=1000, 
        kwargs...
)
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
    elseif sym == :retmsg
        return [u.prob.p.userdata.retmsg for u in sol]
    elseif isnothing(t)
        throw(ArgumentError(
            "Observable $sym is not implemented without a specified time." *
            "Try setting the keyword argument `t`."
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
    kwargs...
)
    if sym == :x
        return [sol(t)[1] for t in times]
    elseif sym == :y
        return [sol(t)[2] for t in times]
    elseif sym == :z
        return [sol(t)[3] for t in times]
    elseif sym == :vparal
        return [sol(t)[4] for t in times]
    elseif sym == :exbdrift
        return [get_exbdrift(sol, t) for t in times]
    elseif sym == :fermi
        return [get_fermi(sol, t) for t in times]
    elseif sym == :betatron
        return [get_betatron(sol, t) for t in times]
    elseif sym == :polarisationacc
        return [get_polarisationacc(sol, t) for t in times]
    elseif sym == :energy
        return [get_energy(sol, t; kwargs...) for t in times]
    elseif sym == :driftenergy
        return [get_driftenergy(sol, t) for t in times]
    elseif sym == :parallelenergy
        return [get_parallelenergy(sol, t) for t in times]
    elseif sym == :perpenergy
        return [get_perpenergy(sol, t) for t in times]
    elseif sym in fieldnames(GCAState)
        charge = sol.prob.p.charge
        mass = sol.prob.p.mass
        magneticmoment = sol.prob.p.magneticmoment
        fields = sol.prob.p.fields
        obs = [
            getfield(
                GCAState(
                    sol(t)[1:4],
                    charge,
                    mass,
                    magneticmoment,
                    fields;
                    time=t,
                    dimensionality="2Dxz",
                ),
                sym)
            for t in times
        ]
        if !isnothing(component)
            obs = [o[component] for o in obs]
        end
        if magnitude
            return norm.(obs)
        else
            return obs
        end
    elseif sym == :scalesratio
        r_L = [
            getfield(
                GCAState(
                    sol(t)[1:4],
                    sol.prob.p.charge,
                    sol.prob.p.mass,
                    sol.prob.p.magneticmoment,
                    sol.prob.p.fields;
                    time=t,
                    dimensionality="2Dxz",
                ),
                :r_L)
            for t in times
        ]
        L_B = [
            getfield(
                GCAState(
                    sol(t)[1:4],
                    sol.prob.p.charge,
                    sol.prob.p.mass,
                    sol.prob.p.magneticmoment,
                    sol.prob.p.fields;
                    time=t,
                    dimensionality="2Dxz",
                ),
                :L_B)
            for t in times
        ]
        return r_L ./ L_B
    end
end
function get_observable(
    sol::ODESolution,
    observable::Function,
    ;
    tspan=(first(sol.t), last(sol.t)),
    nsteps=1000,
    times=range(tspan[1], tspan[2], length=nsteps),
    kwargs...
)
    obs = [
        observable(sol, t; kwargs...)
        for t in times
    ]
    return obs
end


"""
    get_stateidx(
        sol::ODESolution,
        idx::Int;
        tspan=(first(sol.t), last(sol.t)),
        nsteps=1000,
        times=range(tspan[1], tspan[2], length=nsteps),
    )
Return the `idx`-th component of the state vector of `sol` at times `times`.
"""
function get_stateidx(
    sol::ODESolution,
    idx::Int;
    tspan=(first(sol.t), last(sol.t)),
    nsteps=1000,
    times=range(tspan[1], tspan[2], length=nsteps),
)
    return [sol(t)[idx] for t in times]
end

"""
    get_exbdrift(
        sol::ODESolution,
        t::Real,
    )
Calculate the magnitude of the ExB-drift at time `t` of `sol`.
"""
function get_exbdrift(
    sol::ODESolution,
    t::Real,
)
    x, z = sol(t)[1], sol(t)[3]
    B, E = emfieldatpos(
        [x, z],
        sol.prob.p.fields,
    )
    norm(exbdrift(B, E))
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
        sol(t)[1:4],
        q,
        sol.prob.p.mass,
        sol.prob.p.magneticmoment,
        sol.prob.p.fields;
        time=t,
        dimensionality="2Dxz",
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
        sol(t)[1:4],
        q,
        sol.prob.p.mass,
        sol.prob.p.magneticmoment,
        sol.prob.p.fields;
        time=t,
        dimensionality="2Dxz",
    )
    v_gradB = gcastate.∇Bdrift
    E = gcastate.efield
    return v_gradB ⋅ E * q * J2eV
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
        sol(t)[1:4],
        q,
        sol.prob.p.mass,
        sol.prob.p.magneticmoment,
        sol.prob.p.fields;
        time=t,
        dimensionality="2Dxz",
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
    switch=sol.prob.p.switch,
)
    if switch == 2 # Full orbit
        energy = kineticenergy(
            sol(t)[4:6],
            sol.prob.p.mass,
        )
    elseif switch == 1 # GCA
        energy = kineticenergy(
            sol(t)[1:4],
            sol.prob.p.charge,
            sol.prob.p.mass,
            sol.prob.p.magneticmoment,
            sol.prob.p.fields,
            t;
            dimensionality="2Dxz",
        )
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
        q,
        m,
        sol.prob.p.magneticmoment,
        sol.prob.p.fields;
        time=t,
        dimensionality="2Dxz",
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
    x, z = sol(t)[1], sol(t)[3]
    B, _ = emfieldatpos(
        [x, z],
        sol.prob.p.fields,
    )
    return norm(B) * sol.prob.p.magneticmoment * J2eV
end

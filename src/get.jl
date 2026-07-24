#
# get.jl
#
# Created 11.02.25
#
# For fetching data and calculating observables from ODE-solutions and 
# ensembles stored as DataFrames.
#
#____/\_____/\_________________________________________________________________
#
const eomidmap = Dict(
    1 => "GCA",
    2 => "FO",
)

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
    sol::Union{EnsembleSolution,Vector},
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
    elseif sym == :pitchanglef
        return [get_pitchangle(u, last(u.t); kwargs...) for u in sol]
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
    df::DataFrame, # particledata
    sym::Symbol;
    emfield=nothing
)
    npart = nrow(df)
    result = Vector{Float64}(undef, npart)
    EoMsf = [eomidmap[eomid] for eomid in df.eomid]
    EoMs0 = [eomidmap[eomid] for eomid in df.initialeomid]
    if sym == :displacement
        result .= norm.([[
        df.xf[i] - df.x0[i],
        df.yf[i] - df.y0[i],
        df.zf[i] - df.z0[i],
        ] for i in 1:npart])
    elseif sym == :pitchanglef
        if isnothing(emfield)
            error("`emfield` must be provided to calculate pitch angle.")
        end
        for i in 1:npart
            _, Bvec = emfield(df.xf[i], df.yf[i], df.zf[i], df.tf[i])
            if EoMsf[i] == "GCA"
                vparal = df.vxf[i]
                mass = df.mass[i]
                magneticmoment = df.magneticmoment[i]
                result[i] = pitchangle(Bvec, vparal, mass, magneticmoment)
            elseif EoMsf[i] == "FO"
                vel = [df.vxf[i], df.vyf[i], df.vzf[i]]
                result[i] = pitchangle(vel, Bvec)
            end
        end
    elseif sym == :pitchangle0
        if isnothing(emfield)
            error("`emfield` must be provided to calculate pitch angle.")
        end
        for i in 1:npart
            _, Bvec = emfield(df.x0[i], df.y0[i], df.z0[i], df.t0[i])
            if EoMs0[i] == "GCA"
                vparal = df.vx0[i]
                mass = df.mass[i]
                magneticmoment = df.initialmagneticmoment[i]
                result[i] = pitchangle(Bvec, vparal, mass, magneticmoment)
            elseif EoMs0[i] == "FO"
                vel = [df.vx0[i], df.vy0[i], df.vz0[i]]
                result[i] = pitchangle(vel, Bvec)
            end
        end
    else
        throw(ArgumentError("Observable $sym not implemented"))
    end
    return result
end


function get_observable(
    sol::ODESolution,
    sym::Symbol;
    tspan=(first(sol.t), last(sol.t)),
    nsteps=1000,
    times=range(tspan[1], tspan[2], length=nsteps),
    magnitude=false,
    component=nothing,
    EoM=nothing
)


    observables = Dict{Symbol,Function}([
        :vperp => get_vperp;
        [:larmorradius, :r_L] .=> get_larmorradius;
        [:scalesratio, :mgr] .=> get_scalesratio;
        [:magneticcurvatureratio, :mcr] .=> get_mcratio;
        :pitchangle => get_pitchangle;
        :energy => get_energy;
        [:magneticmoment, :mu] .=>
            (sol, t; kwargs...) -> get_magneticmoment(sol, t)
        :gyrofrequency => (sol, t, kwargs...) -> get_gyrofrequency(sol, t);
        :gyroperiod => (sol, t, kwargs...) -> get_gyroperiod(sol, t)
    ])
    gcaobservable = Dict{Symbol,Function}(
        :fermi => get_fermi,
        :betatron => get_betatron,
        :polarisationacc => get_polarisationacc,
        :parallelpower => get_parallelpower,
        :driftenergy => get_driftenergy,
        :parallelenergy => get_parallelenergy,
        :perpenergy => get_perpenergy,
    )
    fieldquantities = Dict{Symbol,Function}([
        :eparal => get_eparal;
        :eperp => get_eperp;
        :efield => (sol, t) -> begin
            sol.prob.p.electromagneticfield(sol(t)[1:3]..., t)[1]
        end
        :bfield => (sol, t) -> begin
            sol.prob.p.electromagneticfield(sol(t)[1:3]..., t)[2]
        end
        :b => (sol, t) -> begin
            bfield = sol.prob.p.electromagneticfield(sol(t)[1:3]..., t)[2]
            bfield / norm(bfield)
        end
        [:efieldratio, :epr] .=> (sol, t) -> begin
            eparal = get_eparal(sol, t)
            eperp = norm(get_eperp(sol, t))
            eparal / eperp
        end
        :exbdrift => (sol, t) -> begin
            exbdrift(
                sol(t, idxs=1:3), t, sol.prob.p.electromagneticfield
            )
        end
        [:characteristicfieldlength, :L_B] .=> (sol, t) -> begin
            characteristicfieldlength(
                sol(t, idxs=1:3), t, sol.prob.p.electromagneticfield
            )
        end
    ])

    # First, check if user is requesting a primary variable
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
    elseif sym == :eomid
        obs = [get_eomid(sol, t) for t in times]
        # Next, check if user is requesting a standard observable
    elseif haskey(observables, sym)
        func = observables[sym]
        # These observables requires knowledge of the EoM and possibly
        # magnetic moment.
        # Check if the solution is using a hybrid GCA/FO scheme with storage
        # of switch times and the EoMs used after each switch. Construct an
        # array of strings that identify the EoM used at each time step. If
        # the solution is not using a hybrid scheme, we just use the EoM
        # specified by the `EoM` keyword argument.
        if (
            hasproperty(sol.prob.p, :eomidafterswitch) &&
            hasproperty(sol.prob.p, :timeatswitch) &&
            hasproperty(sol.prob.p, :magneticmomentafterswitch)
        )
            eomids = [get_eomid(sol, t) for t in times]
            EoMs = [eomidmap[eomid] for eomid in eomids]
            magneticmoments = [
                get_magneticmoment(sol, t; EoM=e) for (t, e) in zip(times, EoMs)
            ]
        elseif EoM == "GCA"
            EoMs = fill("GCA", length(times))
            magneticmoments = fill(sol.prob.p.magneticmoment, length(times))
        elseif EoM == "FO"
            EoMs = fill("FO", length(times))
            magneticmoments = [
                magneticmoment(
                    sol(t, idxs=1:3),
                    sol(t, idxs=4:6),
                    t,
                    sol.prob.p.mass,
                    sol.prob.p.electromagneticfield
                ) for t in times
            ]
        else
            error("EoM $EoM is uknown to `get_observable`.")
        end
        zipper = zip(times, EoMs, magneticmoments)
        obs = [
            func(sol, t, EoM=e, magneticmoment=mu)
            for (t, e, mu) in zipper
        ]
        # Next, check if user is requesting a GCA-only observable
    elseif haskey(gcaobservable, sym)
        @warn "Observable $sym assumes a GCA solution."
        func = gcaobservable[sym]
        obs = [
            func(sol, t) for t in times
        ]
        # Next, check if user is requesting a field quantity
    elseif haskey(fieldquantities, sym)
        func = fieldquantities[sym]
        obs = [func(sol, t) for t in times]
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
    magneticmoment=nothing,
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
            magneticmoment,
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

function findswitchidx(sol::ODESolution, t::Real)
    mask = sol.prob.p.timeatswitch .!= 0
    timeatswitch = sol.prob.p.timeatswitch[mask]
    if length(timeatswitch) != sol.prob.p.nswitches
        error(
            "Switch times do not match the number of switches. Is there" *
            " a switch at time t = 0?"
        )
    end
    i = findfirst(
        x -> x >= t,
        vcat(timeatswitch, last(sol.t)),
    )
    return mask, i
end

function get_eomid(
    sol::ODESolution,
    t::Real;
)
    mask, i = findswitchidx(sol, t)
    eomidafterswitch = sol.prob.p.eomidafterswitch[mask]
    if i == 1
        return sol.prob.p.initialeomid
    else
        return eomidafterswitch[i-1]
    end
end

function get_vperp(
    sol::ODESolution,
    t::Real;
    EoM="GCA",
    magneticmoment=nothing
)
    if EoM == "GCA"
        return perpendicular_velocity(
            sol(t, idxs=1:3),
            t,
            magneticmoment,
            sol.prob.p.mass,
            sol.prob.p.electromagneticfield
        )
    elseif EoM == "FO"
        return perpendicular_velocity(
            sol(t, idxs=1:3),
            sol(t, idxs=4:6),
            t,
            sol.prob.p.electromagneticfield,
        )
    end
end

function get_larmorradius(
    sol::ODESolution,
    t::Real;
    EoM="GCA",
    magneticmoment=nothing
)
    if EoM == "GCA"
        return larmorradius(
            sol(t, idxs=1:3),
            t,
            magneticmoment,
            sol.prob.p.charge,
            sol.prob.p.mass,
            sol.prob.p.electromagneticfield,
        )
    elseif EoM == "FO"
        return larmorradius(
            sol(t, idxs=1:3),
            sol(t, idxs=4:6),
            t,
            sol.prob.p.mass,
            sol.prob.p.charge,
            sol.prob.p.electromagneticfield,
        )
    end
end

function get_scalesratio(
    sol::ODESolution,
    t::Real;
    EoM="GCA",
    magneticmoment=nothing,
)
    if EoM == "GCA"
        return scalesratio(
            sol(t, idxs=1:3),
            t,
            sol.prob.p.mass,
            sol.prob.p.charge,
            magneticmoment,
            sol.prob.p.electromagneticfield,
        )
    elseif EoM == "FO"
        return scalesratio(
            sol(t, idxs=1:3),
            sol(t, idxs=4:6),
            t,
            sol.prob.p.mass,
            sol.prob.p.charge,
            sol.prob.p.electromagneticfield,
        )
    end
end

function get_mcratio(
    sol::ODESolution,
    t::Real;
    EoM="GCA",
    magneticmoment=nothing
)
    if EoM == "GCA"
        return magneticcurvatureratio(
            sol(t, idxs=1:3),
            t,
            sol.prob.p.mass,
            sol.prob.p.charge,
            magneticmoment,
            sol.prob.p.electromagneticfield,
        )
    elseif EoM == "FO"
        return magneticcurvatureratio(
            sol(t, idxs=1:3),
            sol(t, idxs=4:6),
            t,
            sol.prob.p.mass,
            sol.prob.p.charge,
            sol.prob.p.electromagneticfield,
        )
    end
end

function get_pitchangle(
    sol::ODESolution,
    t::Real;
    EoM="GCA",
    magneticmoment=nothing
)
    if EoM == "GCA"
        return pitchangle(
            sol.prob.p.electromagneticfield(sol(t, idxs=1:3)..., t)[2],
            sol(t, idxs=4),
            sol.prob.p.mass,
            magneticmoment,
        )
    elseif EoM == "FO"
        B = sol.prob.p.electromagneticfield(sol(t, idxs=1:3)..., t)[2]
        vel = sol(t, idxs=4:6)
        return pitchangle(vel, B)
    end
end

function get_gyrofrequency(
    sol::ODESolution,
    t::Real;
)
    return norm(sol.prob.p.electromagneticfield(sol(t)[1:3]..., t)[2]) *
           sol.prob.p.charge / sol.prob.p.mass
end

function get_gyroperiod(
    sol::ODESolution,
    t::Real;
)
    return 2π * sol.prob.p.mass /
           (norm(sol.prob.p.electromagneticfield(sol(t)[1:3]..., t)[2]) *
            sol.prob.p.charge)
end

function get_magneticmoment(
    sol::ODESolution,
    t::Real;
    EoM="GCA",
)
    if EoM == "FO"
        return magneticmoment(
            sol(t, idxs=1:3),
            sol(t, idxs=4:6),
            t,
            sol.prob.p.mass,
            sol.prob.p.electromagneticfield
        )
    elseif EoM == "GCA"
        mask, i = findswitchidx(sol, t)
        magneticmomentafterswitch = sol.prob.p.magneticmomentafterswitch[mask]
        if i == 1
            return sol.prob.p.initialmagneticmoment
        else
            return magneticmomentafterswitch[i-1]
        end
    end
end

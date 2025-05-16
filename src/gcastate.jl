"""
    GCAState(
        state            ::Vector{<:Real},
        charge           ::Real,
        mass             ::Real,
        magneticmoment   ::Real,
        time             ::Real,
        kineticenergy    ::Real,
        perpvelocity     ::Real,
        larmorradius     ::Real,
        gyrofrequency    ::Real,
        cosofpitchangle  ::Real,
        exbdrift         ::Vector{<:Real},
        ∇Bdrift          ::Vector{<:Real},
        curvaturedrift   ::Vector{<:Real},
        polarisationdrift::Vector{<:Real},
        mirroracc        ::Real,
        paralacc         ::Real,
        magneticfield    ::Vector{<:Real},
        electricfield    ::Vector{<:Real},
        magneticfieldlength::Real
        )
A struct to hold the state (guiding centre position and its parallel velocity)
of a particle in the guiding centre approximation, all its drifts, the external
electromagnetic field, relevant parameters, and auxiliary quantities.

# Constructors
    GCAState(
        state ::Vector{<:Real},
        q     ::Real,
        m     ::Real,
        μ     ::Real,
        electromagneticfield::Any;
        time=nothing,
        )
Constructs a `GCAState` object from;
* `state` - The guiding centre state; position and parallel velocity
* `vparal` - The parallel velocity of the guiding centre
* `q` - The charge of the particle
* `m` - The mass of the particle
* `μ` - The magnetic moment of the particle
* `electromagneticfield` - A function or functor that returns the
  electric and magnetic field vectors at an arbitrary position in space.
**Keyword** arguments:
* `time` - The time of the state

    GCAState(
        state ::Vector{<:Real},
        q     ::Real, 
        m     ::Real,
        μ     ::Real,
        bfield::Vector{<:Real},
        efield::Vector{<:Real}
        ;
        time=nothing
        )
Constructs a partially filled `GCAState` object, neglecting all drifts except
the ExB-drift.

    GCAState(
        solution::DataFrame,
        fields::Vector{<:AbstractInterpolation}
        ;
        kwargs...
        )
From DataFrame test-particle solution and the field interpolation vector used in
the simulation, construct two ensembles of `GCAState`s. One for the initial
states of the particles and one for the final states.
"""
struct GCAState
    state::Vector{<:Real}
    charge::Real
    mass::Real
    μ::Real
    time::Real
    energy::Real
    vperp::Real
    r_L::Real
    ω_c::Real
    β::Real
    exbdrift::Vector{<:Real}
    ∇Bdrift::Vector{<:Real}
    curvaturedrift::Vector{<:Real}
    polarisationdrift::Vector{<:Real}
    mirroracc::Real
    paralacc::Real
    dbdtacc::Real
    bfield::Vector{<:Real}
    efield::Vector{<:Real}
    dbdt::Vector{<:Real}
    L_B::Real

    function GCAState(
        state::Vector{<:Real},
        q::Real,
        m::Real,
        μ::Real,
        electromagneticfield::Any;
        time=0.0,
    )
        R = state[1:3]
        vparal = state[4]
        q⁻¹ = 1 / q # Inverse of q - to replace division with multiplication
        # Use the gyrocentre position interpolate the vectors
        # scalars from the interpolation objects.
        E_vec, B_vec = electromagneticfield(R...)
        B = norm(B_vec)   # The magnetic field strength
        B⁻¹ = 1 / B         # Inverse of B - to replace divition with multipl.
        b̂ = B_vec * B⁻¹     # An unit vector pointing in the direction of the
        ExBdrift = (E_vec × b̂) / B # The E cross B-drift

        jacobian_matrix = ForwardDiff.jacobian(R) do x
            E_vec_fd, B_vec_fd = electromagneticfield(x...)
            B_fd = norm(B_vec_fd)
            b̂_fd = B_vec_fd / B_fd
            return [b̂_fd; E_vec_fd × b̂_fd / B_fd; B_fd]
        end
        ∇b̂ = jacobian_matrix[1:3, :]
        ∇ExB = jacobian_matrix[4:6, :]
        ∇B = jacobian_matrix[7, :]

        # Total time derivatives. Assumes ∂/∂t = 0,
        dbdt = vparal * (∇b̂ * b̂) + ∇b̂ * ExBdrift
        dExBdt = vparal * (∇ExB * b̂) + ∇ExB * ExBdrift

        # Calculate other drifts
        ∇Bdrift = gradbdrift(b̂, ∇B, μ, B⁻¹, q⁻¹)
        Rdrift = curvaturedrift(b̂, dbdt, vparal, B⁻¹, q⁻¹, m)
        Pdrift = polarisationdrift(b̂, dExBdt, B⁻¹, q⁻¹, m)
        # Compute parallel acceleration
        mirroracc = magneticmirror_acceleration(b̂, ∇B, μ, m)
        paralacc = parallel_acceleration(b̂, E_vec, q, m)
        dbdtacc = fermi_acceleration(ExBdrift, ∇Bdrift, dbdt)
        # Compute the perpendicular velocity
        # Other auxiliary quantities
        vperp = perpendicular_velocity(μ, m, B)
        Ek = kineticenergy(
            vparal, vperp, ExBdrift, ∇Bdrift, Rdrift, Pdrift, m
        )
        L_B = characteristicfieldlength(B, ∇B)
        r_L = larmorradius(m, vperp, q, B)
        ω_c = gyrofrequency(m, q, B)
        β = cosineof_pitchangle(B_vec, vparal, m, μ)

        return new(
            state,
            q, m, μ, time,
            Ek, vperp, r_L, ω_c, β,
            ExBdrift, ∇Bdrift, Rdrift, Pdrift,
            mirroracc, paralacc, dbdtacc,
            B_vec, E_vec,
            dbdt,
            L_B,
        )
    end

    function GCAState(
        state::Vector{<:Real},
        q::Real,
        m::Real,
        μ::Real,
        bfield::Vector{<:Real},
        efield::Vector{<:Real}
        ;
        time=-1
    )
        vparallel = state[4]
        B = norm(bfield)
        v_E = exbdrift(bfield, efield)
        vperp = perpendicular_velocity(μ, m, B)
        Ek = kineticenergy(vparallel, vperp, v_E, m)
        new(
            state, q, m, μ, time, Ek, vperp, -1, -1, -1, v_E,
            -1ones(3), -1ones(3), -1ones(3),
            -1, -1, -1,
            bfield, efield,
            -1ones(3),
            -1
        )
    end

    function GCAState(
        solution::DataFrame,
        electromagneticfield::Any
        ;
        lowmem=true,
        kwargs...
    )
        return init_and_final_ensembleofgcastates(
            solution, electromagneticfield; lowmem=lowmem, kwargs...
        )
    end
end


################################################################################
"""
Get functions
"""

"""
    get_drifts(gcastate::GCAState)
Return all drifts components in a vector.
"""
function get_drifts(gcastate::GCAState)
    return [
        gcastate.exbdrift,
        gcastate.∇Bdrift,
        gcastate.curvaturedrift,
        gcastate.polarisationdrift
    ]
end


#-------------------------------------------------------------------------------
# To retrieve various quantities from a vector GCAStates
function get_x(gcastates::Vector{GCAState})
    [gcastate.state[1] for gcastate in gcastates]
end
function get_y(gcastates::Vector{GCAState})
    [gcastate.state[2] for gcastate in gcastates]
end
function get_z(gcastates::Vector{GCAState})
    [gcastate.state[3] for gcastate in gcastates]
end
function get_vparal(gcastates::Vector{GCAState})
    [gcastate.state[4] for gcastate in gcastates]
end
function get_mu(gcastates::Vector{GCAState})
    [gcastate.μ for gcastate in gcastates]
end
function get_time(gcastates::Vector{GCAState})
    [gcastate.time for gcastate in gcastates]
end
function get_energy(gcastates::Vector{GCAState})
    [gcastate.energy for gcastate in gcastates]
end
function get_vperp(gcastates::Vector{GCAState})
    [gcastate.vperp for gcastate in gcastates]
end
function get_larmorradius(gcastates::Vector{GCAState})
    [gcastate.r_L for gcastate in gcastates]
end
function get_beta(gcastates::Vector{GCAState})
    [gcastate.β for gcastate in gcastates]
end
function get_gyrofrequency(gcastates::Vector{GCAState})
    [gcastate.ω_c for gcastate in gcastates]
end
function get_bfield(gcastates::Vector{GCAState})
    [gcastate.bfield for gcastate in gcastates]
end
function get_B(gcastates::Vector{GCAState})
    [norm(gcastate.bfield) for gcastate in gcastates]
end
function get_E(gcastates::Vector{GCAState})
    [norm(gcastate.efield) for gcastate in gcastates]
end
function get_efield(gcastates::Vector{GCAState})
    [gcastate.efield for gcastate in gcastates]
end
function get_fieldlength(gcastates::Vector{GCAState})
    [gcastate.L_B for gcastate in gcastates]
end
function get_scalesratio(gcastates::Vector{GCAState})
    [gcastate.r_L / gcastate.L_B for gcastate in gcastates]
end
function get_exb(gcastates::Vector{GCAState})
    [norm(gcastate.exbdrift) for gcastate in gcastates]
end
function get_gradBdrift(gcastates::Vector{GCAState})
    [norm(gcastate.∇Bdrift) for gcastate in gcastates]
end
function get_Rdrift(gcastates::Vector{GCAState})
    [norm(gcastate.curvaturedrift) for gcastate in gcastates]
end
function get_Pdrift(gcastates::Vector{GCAState})
    [norm(gcastate.polarisationdrift) for gcastate in gcastates]
end
function get_mirroracc(gcastates::Vector{GCAState})
    [gcastate.mirroracc for gcastate in gcastates]
end
function get_paralacc(gcastates::Vector{GCAState})
    [gcastate.paralacc for gcastate in gcastates]
end
function get_mass(gcastates::Vector{GCAState})
    gcastates[1].mass
end
function get_drifts(gcastates::Vector{GCAState})
    [
        get_exb(gcastates),
        get_gradBdrift(gcastates),
        get_Rdrift(gcastates),
        get_Pdrift(gcastates),
    ]
end


"""
    get_driftnorm(gcastates::Vector{GCAState})
Return the total drift velocity at each state in a vector of GCAStates.
"""
function get_driftnorm(gcastates::Vector{GCAState})
    [norm(
        gcastate.exbdrift .+
        gcastate.∇Bdrift .+
        gcastate.curvaturedrift .+
        gcastate.polarisationdrift
    ) for gcastate in gcastates]
end


#-------------------------------------------------------------------------------
# To retrieve various quantities from a vector of a vector of GCAStates.
function get_initialenergies(states::Vector{Vector{GCAState}})
    [s[1].energy for s in states]
end
function get_finalenergies(states::Vector{Vector{GCAState}})
    [s[end].energy for s in states]
end
function get_initialx(states::Vector{Vector{GCAState}})
    [s[1].state[1] for s in states]
end
function get_finalx(states::Vector{Vector{GCAState}})
    [s[end].state[1] for s in states]
end
function get_initialz(states::Vector{Vector{GCAState}})
    [s[1].state[3] for s in states]
end
function get_finalz(states::Vector{Vector{GCAState}})
    [s[end].state[3] for s in states]
end


################################################################################
"""
    DataFrames.DataFrame(gcastates::Vector{GCAState})
Constructs a DataFrame from a vector of `GCAState` objects. The DataFrame
contains the fields of the `GCAState` objects as columns. Vector fields are
collapsed into one column per component.
"""
function DataFrames.DataFrame(gcastates::Vector{GCAState}; lowmem=false)
    df = DataFrame(
        Dict(
            field => [getfield(state, field) for state in gcastates]
            for field in fieldnames(GCAState)
        )
    )
    newcols = Dict(
        :bx => (:bfield, 1),
        :by => (:bfield, 2),
        :bz => (:bfield, 3),
        :ex => (:efield, 1),
        :ey => (:efield, 2),
        :ez => (:efield, 3),
        :exbdrift_x => (:exbdrift, 1),
        :exbdrift_y => (:exbdrift, 2),
        :exbdrift_z => (:exbdrift, 3),
        :∇Bdrift_x => (:∇Bdrift, 1),
        :∇Bdrift_y => (:∇Bdrift, 2),
        :∇Bdrift_z => (:∇Bdrift, 3),
        :curvaturedrift_x => (:curvaturedrift, 1),
        :curvaturedrift_y => (:curvaturedrift, 2),
        :curvaturedrift_z => (:curvaturedrift, 3),
        :polarisationdrift_x => (:polarisationdrift, 1),
        :polarisationdrift_y => (:polarisationdrift, 2),
        :polarisationdrift_z => (:polarisationdrift, 3),
        :dbdt_x => (:dbdt, 1),
        :dbdt_y => (:dbdt, 2),
        :dbdt_z => (:dbdt, 3),
    )
    if !lowmem
        newcols[:Rx] = (:state, 1)
        newcols[:Ry] = (:state, 2)
        newcols[:Rz] = (:state, 3)
        newcols[:vparal] = (:state, 4)
    end

    # Create new columns for the vector field components
    for (field, (source, idx)) in newcols
        column = stack(df[:, source])
        data = column[idx, :]
        df[!, field] = data
    end
    # Remove the vector fields
    for source in values(newcols)
        try
            select!(df, Not(source[1]))
        catch
        end
    end
    # Change order of columns
    if !lowmem
        select!(df,
            [
                :Rx, :Ry, :Rz, :vparal, :charge, :mass, :μ, :time, :energy, :vperp,
                :r_L, :ω_c, :β, :exbdrift_x, :exbdrift_y,
                :exbdrift_z, :∇Bdrift_x, :∇Bdrift_y, :∇Bdrift_z, :curvaturedrift_x,
                :curvaturedrift_y, :curvaturedrift_z, :polarisationdrift_x,
                :polarisationdrift_y, :polarisationdrift_z, :mirroracc, :paralacc,
                :dbdtacc, :bx, :by, :bz, :ex, :ey, :ez, :dbdt_x, :dbdt_y,
                :dbdt_z, :L_B
            ]
        )
    else
        select!(df,
            [
                :energy, :vperp,
                :r_L, :ω_c, :β, :exbdrift_x, :exbdrift_y,
                :exbdrift_z, :∇Bdrift_x, :∇Bdrift_y, :∇Bdrift_z, :curvaturedrift_x,
                :curvaturedrift_y, :curvaturedrift_z, :polarisationdrift_x,
                :polarisationdrift_y, :polarisationdrift_z, :mirroracc, :paralacc,
                :dbdtacc, :bx, :by, :bz, :ex, :ey, :ez, :dbdt_x, :dbdt_y,
                :dbdt_z, :L_B
            ]
        )
    end
    return df
end


"""
    timeseriesofgcastates(
        solution::ODESolution,
        times::AbstractVector=solution.t
        ;
        components="all",
        )
Constructs a vector of `GCAState` objects from an `ODESolution` object.
Exploits the inherent solution interpolation in a `ODESolution` object to
construct `GCAState`s at arbitrary times.

# Keyword arguments
- `components`: The components of the guiding centre drift to include in the
  in the `GCAState` objects. Default: "all"
"""
function timeseriesofgcastates(
    solution::ODESolution,
    times::AbstractVector=solution.t
    ;
    components="all",
)
    timeseries = Vector{GCAState}(undef, length(times))
    Threads.@threads for i in eachindex(times)
        interpolated_solution = solution(times[i], idxs=1:4)
        if components == "all"
            timeseries[i] = GCAState(
                interpolated_solution,
                solution.prob.p.charge,
                solution.prob.p.mass,
                solution.prob.p.magneticmoment,
                solution.prob.p.fields,
                ;
                time=times[i],
            )
        else
            itpvec = solution.prob.fields
            R = interpolated_solution[1:3]
            bfield = [itpvec[i](R...) for i in 1:3]
            efield = [itpvec[i](R...) for i in 4:6]
            timeseries[i] = GCAState(
                interpolated_solution,
                solution.prob.p[1],
                solution.prob.p[2],
                solution.prob.p[3],
                bfield,
                efield
                ;
                time=times[i]
            )
        end
    end
    return timeseries
end


"""
    init_and_final_ensembleofgcastates(
        solution::Vector{Tuple{Vararg{Any}}},
        electromagneticfield::Any,
        ;
        components="all",
        charge=tp.e,
        mass=tp.m_e,
        )
Caclulate ensembles of `GCAState` objects from the initial and final states of
an ensemble solution, using the electromangetic field function/functor
`electromagneticfield`.

The solution is assumed to be a vector of a tuple of any type
(e.g. an output function). The function needs to be able to extract the initial
and final states of the solution using the `getinitialstates` and
`getfinalstates` functions. The function also needs to be able to extact the
magnetic moment out of the solution using the function `getmagneticmoments`.

Set components to "all" to calculate include all drift components in the
calculation. Otherwise, the function will only use the ExB-drift component.
"""
function init_and_final_ensembleofgcastates(
    solution::Any,
    electromagneticfield::Any
    ;
    kwargs...
)
    initialstates = getinitialstate(solution)
    finalstates = getfinalstate(solution)
    initialtimes = getinitialtimes(solution)
    finaltimes = getfinaltimes(solution)
    magneticmoments = getmagneticmoments(solution)
    initialgcastates = ensembleofgcastates(
        initialstates,
        magneticmoments,
        initialtimes,
        electromagneticfield
        ;
        kwargs...
    )
    finalgcastates = ensembleofgcastates(
        finalstates,
        magneticmoments,
        finaltimes,
        electromagneticfield
        ;
        kwargs...
    )
    return initialgcastates, finalgcastates
end


"""
    ensembleofgcastates(
        solution::DataFrame,
        magneticmoments::Vector{<:Real},
        times::Vector{<:Real},
        args...
        ;
        kwargs...
        )
    ensembleofgcastates(
        solution::Vector{<:Vector{<:Real}},
        magneticmoments::Vector{<:Real},
        times::Vector{<:Real},
        electromagneticfield::Any
        ;
        components="exb",
        charge=tp.e,
        mass=tp.m_e,
        )

Construct a ensemble of `GCAState`s from a `solution` of a simulation using the
GCA. The `solution` contains multiple GCA state-vectors. `magneticmoments` and
`times` contain their respective magnetic moments and times.

Require an callable object that returns the electromagnetic field at a given
position in 3D space.

________________________________________________________________________________

    ensembleofgcastates(
        solution::Vector{Tuple{Any, Any, Any, Vector, Vector, Vector, Vector, Any, Any}}
        ;
        components="exb",
        charge=tp.e,
        mass=tp.m_e,
        )
Harcoded for evaluating a particular deprecated output function which stored
a vector of (u0, uf, mu, B0, E0, Bf, Ef, nt, time)

"""
function ensembleofgcastates(
    solution::DataFrame,
    magneticmoments::Vector{<:Real},
    times::Vector{<:Real},
    args...
    ;
    lowmem=true,
    kwargs...
)
    nrows = size(solution, 1)
    gcastates = ensembleofgcastates(
        [Vector(solution[i, :]) for i in 1:nrows],
        magneticmoments,
        times,
        args...
        ;
        kwargs...
    )
    DataFrame(gcastates; lowmem=lowmem)
end
function ensembleofgcastates(
    solution::Vector{<:Vector{<:Real}},
    magneticmoments::Vector{<:Real},
    times::Vector{<:Real},
    electromagneticfield::Any
    ;
    components="all",
    charge=e,
    mass=m_e,
)
    nofparticles = length(solution)
    particlestates = Vector{GCAState}(undef, nofparticles)
    for i in 1:nofparticles
        if components == "all"
            particlestates[i] = GCAState(
                solution[i],
                charge,
                mass,
                magneticmoments[i],
                electromagneticfield
                ;
                time=times[i],
            )
        elseif components == "exb"
            R = solution[i][1:3]
            efield, bfield = electromagneticfield(R...)
            particlestates[i] = GCAState(
                solution[i],
                charge,
                mass,
                magneticmoments[i],
                bfield,
                efield
                ;
                time=times[i]
            )
        end
    end
    return particlestates
end

#____/\_____/\_________________________________________________________________
#
# Functions uses GCAState types
#

"""
    kineticenergy(
        state::Vector{<:Real}, # Guiding centre state
        q::Real, # Particle charge
        m::Real, # Particle mass
        μ::Real, # Particle magnetic moment
        electromagneticfield::Any,
)
Kinetic energy of a charged particle based on its guiding centre, parallel
velocity, magnetic moment, and the electromagnetic field in which it is
embedded.
"""
function kineticenergy(
    state::Vector{<:Real}, # Guiding centre state
    q::Real, # Particle charge
    m::Real, # Particle mass
    μ::Real, # Particle magnetic moment
    electromagneticfield::Any,
    time::Real
)
    kineticenergy(GCAState(
        state,
        q,
        m,
        μ,
        electromagneticfield;
        time=time,
    ))
end
"""
    kineticenergy(
        gcastate::GCAState
    )
Kinetic energy of a charged particle based on its guiding centre state
represented by a `GCAState` object.
"""
function kineticenergy(gcastate::GCAState)
    kineticenergy(
        gcastate.state[4],
        gcastate.vperp,
        gcastate.exbdrift,
        gcastate.∇Bdrift,
        gcastate.curvaturedrift,
        gcastate.polarisationdrift,
        gcastate.mass
    )
end

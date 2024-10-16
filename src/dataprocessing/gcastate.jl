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
        itpvec::Vector{<:AbstractInterpolation}; 
        time=nothing,
        dimensionality="2Dxz"
        )
Constructs a GCAState object from;
* `state` - The guiding centre state; position and parallel velocity
* `vparal` - The parallel velocity of the guiding centre
* `q` - The charge of the particle
* `m` - The mass of the particle
* `μ` - The magnetic moment of the particle
* `itpvec` - A vector of interpolation objects for the electromagnetic fields
Optional kqyword arguments;
* `time` - The time of the state
* `dimensionality` - The dimensionality of the problem (default: "2Dxz")

_______________________________________________________________________________

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
Constructs a partially filled GCAState object, neglecting all drifts except
the ExB-drift.

_______________________________________________________________________________

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
    L_B::Real

    function GCAState(
        state ::Vector{<:Real},
        q     ::Real,
        m     ::Real,
        μ     ::Real, 
        itpvec::Vector{<:AbstractInterpolation}; 
        time=0.0,
        dimensionality="2Dxz"
        )
        if dimensionality == "2Dxz"
            R = [state[1], state[3]]
        else
            R = state[1:3]
        end
        vparal = state[4]
        q⁻¹ = 1/q # Inverse of q - to replace division with multiplication
        # Use the gyrocentre position interpolate the vectors
        # scalars from the interpolation objects.
        B_vec = [itpvec[i](R...) for i in 1:3]
        E_vec = [itpvec[i](R...) for i in 4:6]
        B = norm(B_vec)   # The magnetic field strength
        B⁻¹ = 1/B         # Inverse of B - to replace divition with multipl.
        b̂ = B_vec*B⁻¹     # An unit vector pointing in the direction of the
        ExBdrift = (E_vec × b̂)/B # The E cross B-drift
    
        if length(itpvec) == 27
            ∇B = [itpvec[i](R...) for i in 7:9]
            ∇b̂ = reshape([itpvec[i](R...) for i in 10:18], 3, 3)
            ∇ExB = reshape([itpvec[i](R...) for i in 19:27], 3, 3)
        elseif length(itpvec) == 6
            jacobian_matrix = ForwardDiff.jacobian(R) do x
                vec = [itp(x...) for itp in itpvec]
                B_vec_fd = vec[1:3]
                E_vec_fd = vec[4:6]
                B_fd = norm(B_vec_fd)
                b̂_fd = B_vec_fd/B_fd
                return [b̂_fd; E_vec_fd × b̂_fd/B_fd; B_fd]
            end
            # Add zeros-column representing derivatives along the y-axis
            if dimensionality == "2Dxz"
                jacobian_matrix = [
                    jacobian_matrix[:,1];;
                    zeros(eltype(R), 7);;
                    jacobian_matrix[:,2]
                    ]
            end
            ∇b̂ = jacobian_matrix[1:3,:]
            ∇ExB = jacobian_matrix[4:6,:]
            ∇B = jacobian_matrix[7,:]
        end
    
        # Total time derivatives. Assumes ∂/∂t = 0,
        dbdt = vparal * (∇b̂ * b̂) + ∇b̂*ExBdrift
        dExBdt = vparal * (∇ExB * b̂) + ∇ExB*ExBdrift
    
        # Calculate other drifts
        ∇Bdrift = gradbdrift(b̂, ∇B, μ, B⁻¹, q⁻¹)
        Rdrift = curvaturedrift(b̂, dbdt, vparal, B⁻¹, q⁻¹, m)
        Pdrift = polarisationdrift(b̂, dExBdt, B⁻¹, q⁻¹, m)
        # Compute parallel acceleration
        mirroracc = magneticmirror_acceleration(b̂, ∇B, μ, m)
        paralacc = parallel_acceleration(b̂, E_vec, q, m)
        dbdtacc = drift_dbdt_acceleration(ExBdrift, ∇Bdrift, dbdt)
        # Compute the perpendicular velocity
        # Other auxiliary quantities
        vperp = √(2B*μ/m) # The perpendicular velocity
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
            L_B,
            )
    end

   function GCAState(
        state    ::Vector{<:Real},
        q        ::Real,
        m        ::Real,
        μ        ::Real,
        bfield   ::Vector{<:Real},
        efield   ::Vector{<:Real}
        ;
        time=-1
        )
        vparallel = state[4]
        v_E, B, _ = exbdrift(bfield, efield)
        vperp = perpendicular_velocity(μ, m, B)
        Ek = kineticenergy(vparallel, vperp, v_E, m)
        new(
            state, q, m, μ, time, Ek, vperp, -1,-1, -1, v_E,
            -1ones(3), -1ones(3), -1ones(3), -1, -1, -1, bfield, efield, -1
            )
    end

    function GCAState(
        solution::DataFrame,
        fields::Vector{<:AbstractInterpolation}
        ;
        lowmem=true,
        kwargs...
        )
        return init_and_final_ensembleofgcastates(
            solution, fields; lowmem=lowmem, kwargs...
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
    [gcastate.r_L/gcastate.L_B for gcastate in gcastates]
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


"""
    get_forces(gcastates::Vector{GCAState})
Calculate and return the forces acting on a GCA-particle. The parallel forces
are straightforward, since we already have the parallel acceleration terms in
the `GCAState`. The perpendicular drift forces and the "gyration force" needs
to be evaluated through differentiation. Since the latter forces are evaluated
at a time between the `GCAState`s, the parallel forces are interpolated in
time.
"""
function get_forces(gcastates::Vector{GCAState})
    # Perpendicular acceleration
    # Use first order upwind finite difference scheme to evaluate acceleration
    drifts = get_drifts(gcastates)
    times = get_time(gcastates)
    Δdrifts = [diff(drift) for drift in drifts]
    Δt = diff(times)
    perpendicular_acceleration = [ Δv ./ Δt for Δv in Δdrifts]

    # "Gyration force"
    Δvperp = diff(get_vperp(gcastates))
    append!(perpendicular_acceleration, [Δvperp ./ Δt])

    # Energy change
    ΔE = diff(get_energy(gcastates))
    ΔEΔt = ΔE ./ Δt
    # Take the pairwise average of the parallel acceleration to get the
    # acceleration at the same time.
    pairwiseaverage(x) = ((x .+ circshift(x, -1))/2)[1:end-1]
    mirroracc = pairwiseaverage(get_mirroracc(gcastates))
    parallelacc = pairwiseaverage(get_paralacc(gcastates))
    acceleration = [[parallelacc, mirroracc]; perpendicular_acceleration]
    force = get_mass(gcastates)*acceleration

    # Return a NamedTuple of the time and force components.
    return (
        time = pairwiseaverage(times),
        parallel = force[1],
        mirror = force[2],
        exb = force[3],
        gradb = force[4],
        curvature = force[5],
        polarisation = force[6],
        vperp = force[7],
        energy = ΔEΔt
    )
end

#-------------------------------------------------------------------------------
# To retrieve various quantities from a vector of a vector of GCAStates.
function get_initialenergies(states::Vector{Vector{GCAState}})
    e0 = [s[1].energy for s in states]
end
function get_finalenergies(states::Vector{Vector{GCAState}})
    ef = [s[end].energy for s in states]
end
function get_initialx(states::Vector{Vector{GCAState}})
    e0 = [s[1].state[1] for s in states]
end
function get_finalx(states::Vector{Vector{GCAState}})
    ef = [s[end].state[1] for s in states]
end
function get_initialz(states::Vector{Vector{GCAState}})
    e0 = [s[1].state[3] for s in states]
end
function get_finalz(states::Vector{Vector{GCAState}})
    ef = [s[end].state[3] for s in states]
end


################################################################################
"""
    DataFrame(gcastates::Vector{GCAState})
Constructs a DataFrame from a vector of `GCAState` objects. The DataFrame
contains the fields of the `GCAState` objects as columns. Vector fields are
collapsed into one column per component.
"""
function DataFrame(gcastates::Vector{GCAState}; lowmem=false)
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
        try select!(df, Not(source[1]))
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
            :dbdtacc, :bx, :by, :bz, :ex, :ey, :ez, :L_B
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
            :dbdtacc, :bx, :by, :bz, :ex, :ey, :ez, :L_B
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
        dimensionality="2Dxz"
        )
Constructs a vector of `GCAState` objects from an `ODESolution` object.
Exploits the inherent solution interpolation in a `ODESolution` object to
construct `GCAState`s at arbitrary times.

# Keyword arguments
- `components`: The components of the guiding centre drift to include in the
  in the `GCAState` objects. Default: "all"
- `dimensionality`: The dimensionality of the background fields. I.e. the
dimension of the interpolation objects in the solution parameters.
Default: "2Dxz".
"""
function timeseriesofgcastates(
    solution::ODESolution,
    times::AbstractVector=solution.t
    ;
    components="all",
    dimensionality="2Dxz"
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
                dimensionality=dimensionality
                )
        else
            itpvec = solution.prob.fields
            if dimensionality == "2Dxz"
                Rx = interpolated_solution[1]
                Rz = interpolated_solution[3]
                bfield = [itpvec[i](Rx, Rz) for i in 1:3]
                efield = [itpvec[i](Rx, Rz) for i in 4:6]
            else
                error("Dimensionality not implemented")
            end
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
        itpvec::Vector{<:AbstractInterpolation},
        ;
        components="all",
        charge=tp.e,
        mass=tp.m_e,
        )
Caclulate ensembles of `GCAState` objects from the initial and final states of
an ensemble solution, using the electromangetic field interpolation objects
`itpvec`.

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
    itpvec::Vector{<:AbstractInterpolation},
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
        itpvec
        ;
        kwargs...
        )
    finalgcastates = ensembleofgcastates(
        finalstates,
        magneticmoments,
        finaltimes,
        itpvec
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
        itpvec::Vector{<:AbstractInterpolation},
        ;
        components="exb",
        charge=tp.e,
        mass=tp.m_e,
        dimensionality="2Dxz"
        )

Construct a ensemble of `GCAState`s from a `solution` of a simulation using the
GCA. The `solution` contains multiple GCA state-vectors. `magneticmoments` and
`times` contain their respective magnetic moments and times.

Require the electromagnetic field interpolation object vector as a additional
argument.

________________________________________________________________________________

    ensembleofgcastates(
        solution::Vector{Tuple{Any, Any, Any, Vector, Vector, Vector, Vector, Any, Any}}
        ;
        components="exb",
        charge=tp.e,
        mass=tp.m_e,
        dimensionality="2Dxz"
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
        [Vector(solution[i,:]) for i in 1:nrows],
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
    itpvec::Vector{<:AbstractInterpolation},
    ;
    components="exb",
    charge=tp.e,
    mass=tp.m_e,
    dimensionality="2Dxz",
    lowmem=false,
    )
    nofparticles = length(solution)
    particlestates = Vector{GCAState}(undef, nofparticles)
    Threads.@threads for i in 1:nofparticles
        if components == "all"
            particlestates[i] = GCAState(
                solution[i],
                charge,
                mass,
                magneticmoments[i],
                itpvec
                ;
                time=times[i],
                dimensionality=dimensionality,

                )
        elseif components == "exb"
            if dimensionality == "2Dxz"
                Rx = solution[i][1]
                Rz = solution[i][3]
                bfield = [itpvec[i](Rx, Rz) for i in 1:3]
                efield = [itpvec[i](Rx, Rz) for i in 4:6]
            else
                error("Dimensionality not implemented")
            end
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

function ensembleofgcastates(
    solution::Vector{Tuple{Any, Any, Any, Vector, Vector, Vector, Vector, Any, Any}}
    ;
    components="exb",
    charge=tp.e,
    mass=tp.m_e,
    dimensionality="2Dxz"
    )
    nofparticles = length(solution)
    initialstates = Vector{GCAState}(undef, nofparticles)
    finalstates = Vector{GCAState}(undef, nofparticles)
    for i in 1:nofparticles
        if components == "exb"
            initialstates[i] = GCAState(
                solution[i],
                charge,
                mass,
                solution[i][3], # magnetic moment
                solution[i][4], # magnetic field
                solution[i][5], # electric field
                ;
                time=solution[i][9],
                )
            finalstates[i] = GCAState(
                solution[i],
                charge,
                mass,
                solution[i][3], # magnetic moment
                solution[i][6], # magnetic field
                solution[i][7], # electric field
                ;
                time=solution[i][9],
                )
        else 
            error("Field interpolation objects unknown")
        end
    end
    return initialstates, finalstates
end

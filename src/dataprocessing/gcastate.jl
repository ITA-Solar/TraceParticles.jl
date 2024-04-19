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
        exbdrift         ::Vector{<:Real},
        ∇Bdrift          ::Vector{<:Real},
        curvaturedrift   ::Vector{<:Real},
        polarisationdrift::Vector{<:Real},
        mirroracc        ::Real,
        parallelacc      ::Real,
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
Constructs a partially filled GCAState object, neglectin all drifts except
the ExB-drift.
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
    exbdrift::Vector{<:Real}
    ∇Bdrift::Vector{<:Real}
    curvaturedrift::Vector{<:Real}
    polarisationdrift::Vector{<:Real}
    mirroracc::Real
    paralacc::Real
    bfield::Vector{<:Real}
    efield::Vector{<:Real}
    L_B::Real

    function GCAState(
        state ::Vector{<:Real},
        q     ::Real,
        m     ::Real,
        μ     ::Real, 
        itpvec::Vector{<:AbstractInterpolation}; 
        time=nothing,
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
            # Calculate the gradient of the magnetic field strength
            ∇B = ForwardDiff.gradient(R) do x
                sqrt(itpvec[1](x...)^2 + itpvec[2](x...)^2 + itpvec[3](x...)^2)
            end
            ∇B = [∇B[1], 0f0, ∇B[2]]
            #
            # Calculate the gradient of the magnetic field direction
            jacobian_matrix = stack(
                [Interpolations.gradient(itp, R...) for itp in itpvec],
                dims=1
                )
            # Add zeros-column representing derivatives along the y-axis
            if dimensionality == "2Dxz"
                jacobian_matrix = [
                    jacobian_matrix[:,1];;
                    zeros(eltype(R), 6);;
                    jacobian_matrix[:,2]
                    ]
            end
            ∇B_vec = jacobian_matrix[1:3,:]
            ∇b̂ = (∇B_vec - b̂ * ∇B')*B⁻¹
            #
            # Calculate the Jacobian matrix of the ExB-drift
            ∇E_vec = jacobian_matrix[4:6,:]
            skewE = skewsymmetric_matrix(E_vec)
            skewb = skewsymmetric_matrix(b̂)
            ∇ExB = (-skewb*∇E_vec + skewE*∇b̂ - ExBdrift * ∇B')*B⁻¹
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
        # Compute the perpendicular velocity
        # Other auxiliary quantities
        vperp = √(2B*μ/m) # The perpendicular velocity
        Ek = kineticenergy(
            vparal, vperp, ExBdrift, ∇Bdrift, Rdrift, Pdrift, m
            )
        L_B = characteristicfieldlength(B, ∇B)
        r_L = larmorradius(m, vperp, q, B)
        ω_c = gyrofrequency(m, q, B)
    
        return new(
            state,
            q, m, μ, time,
            Ek, vperp, r_L, ω_c,
            ExBdrift, ∇Bdrift, Rdrift, Pdrift,
            mirroracc, paralacc,
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
            state, q, m, μ, time, Ek, vperp, -1,-1, v_E,
            [-1], [-1], [-1], -1, -1, bfield, efield, -1
            )
    end
end


function get_x(gcastates::Vector{GCAState})
    [gcastate.state[1] for gcastate in gcastates]
end
function get_y(gcastates::Vector{GCAState})
    [gcastate.state[2] for gcastate in gcastates]
end
function get_z(gcastates::Vector{GCAState})
    [gcastate.state[3] for gcastate in gcastates]
end
function get_vparallel(gcastates::Vector{GCAState})
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
function get_gyrofrequency(gcastates::Vector{GCAState})
    [gcastate.ω_c for gcastate in gcastates]
end
function get_bfield(gcastates::Vector{GCAState})
    [gcastate.bfield for gcastate in gcastates]
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
function get_drifts(gcastate::GCAState)
    return [
        gcastate.exbdrift,
        gcastate.∇Bdrift,
        gcastate.curvaturedrift,
        gcastate.polarisationdrift
        ]
end


"""
    kineticenergy(
        gcastate::GCAState,
        )
Calculates the kinetic energy of a charged particle given the state of the
particle in the guiding centre approximation.
"""
function kineticenergy(
    gcastate::GCAState,
    )
    if any(isnothing, get_drifts(gcastate))
        kineticenergy(
            gcastate.state[4],
            gcastate.vperp,
            gcastate.exbdrift,
            gcastate.mass
            )
    else
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
    for i in eachindex(times)
        interpolated_solution = solution(times[i], idxs=1:4)
        if components == "all"
            timeseries[i] = GCAState(
                interpolated_solution,
                solution.prob.p[1],
                solution.prob.p[2],
                solution.prob.p[3],
                solution.prob.p[4]
                ;
                time=times[i],
                dimensionality=dimensionality
                )
        else
            itpvec = solution.prob.p[4]
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
"""
function ensembleofgcastates(
    solution::Vector{Vector{<:Real}},
    magneticmoments::Vector{<:Real},
    time::Vector{<:Real},
    itpvec::Vector{<:AbstractInterpolation},
    ;
    components="exb",
    charge=tp.e,
    mass=tp.m_e,
    dimensionality="2Dxz"
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
                itpvec
                ;
                time=time[i],
                dimensionality=dimensionality
                )
        elseif components == "exb"
            bfield = [itpvec[i](solution[i][1]...) for i in 1:3]
            efield = [itpvec[i](solution[i][1]...) for i in 4:6]
            particlestates[i] = GCAState(
                solution[i],
                charge,
                mass,
                magneticmoments[i],
                bfield,
                efield
                ;
                time=time[i]
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
    solution::Vector{Tuple{Vararg{Any}}},
    itpvec::Vector{<:AbstractInterpolation},
    ;
    components="",
    charge=tp.e,
    mass=tp.m_e,
    )
    initialstates = getinitialstates(solution)
    finalstates = getfinalstates(solution)
    initialtimes = getinitialtimes(solution)
    finaltimes = getfinaltimes(solution)
    magneticmoments = getmagneticmoments(solution)
    initialgcastates = ensembleofgcastates(
        initialstates,
        magneticmoments,
        initialtimes,
        itpvec
        ;
        components=components, charge=charge, mass=mass
        )
    finalgcastates = ensembleofgcastates(
        finalstates, 
        magneticmoments,
        finaltimes,
        itpvec
        ;
        components=components, charge=charge, mass=mass
        )
    return initialgcastates, finalgcastates
end

# Created 15.01.24 by eilfso
#
#
#-------------------------------------------------------------------------------
# Created 15.01.24
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                physics.jl
#
#-------------------------------------------------------------------------------
#  Contains physics of charge particles in plasma.
#-------------------------------------------------------------------------------

"""
    larmorradius(mass, vperp, charge, magneticfieldstrength)
Return the Larmor radius of a charged particle.
"""
function larmorradius(mass, vperp::Real, charge, magneticfieldstrength::Real)
    return abs(mass * vperp / (charge * magneticfieldstrength))
end
function larmorradius(
    R::AbstractVector,
    μ::Real,
    q::Real,
    m::Real,
    itpvec::Vector{<:Any}
)
    B = norm([itp(R...) for itp in itpvec])
    vperp = perpendicular_velocity(μ, m, B)
    return larmorradius(m, vperp, q, B)
end
function larmorradius(
    R::AbstractVector,
    params::NamedTuple
)
    B = norm(params.electromagneticfield(R...)[2])
    μ = params.magneticmoment
    q = params.charge
    m = params.mass
    vperp = perpendicular_velocity(μ, m, B)
    return larmorradius(m, vperp, q, B)
end

"""
    gyrofrequency(mass, charge, magneticfieldstrength)
Return the gyrofrequency of a charged particle.
"""
function gyrofrequency(mass, charge, magneticfieldstrength)
    return abs(charge * magneticfieldstrength / mass)
end


"""
Return the perpendicular velocity of a charged particle in a electromagnetic
field.

# Methods
    perpendicular_velocity(velocity, magneticfield, electricfield)
    perpendicular_velocity(magnetic_moment, mass, magneticfieldstrength)
    perpendicular_velocity(velocity, b_vec, vparal)
"""
function perpendicular_velocity(
    magnetic_moment::Real,
    mass::Real,
    magneticfieldstrength::Real
)::Real
    sqrt(2magnetic_moment * magneticfieldstrength / mass)
end
function perpendicular_velocity(
    vel::AbstractVector,
    magneticfield::AbstractVector,
    electricfield::AbstractVector,
)
    B = norm(magneticfield) # Magnetic field strength
    b_vec = magneticfield / B   # Magnetic field direction (unit vector)
    ExBdrift = exbdrift(magneticfield, electricfield) # get E cross B drift,
    vel_in_E_frame = vel - ExBdrift
    vparal = vel ⋅ b_vec
    vperp = perpendicular_velocity(vel_in_E_frame, b_vec, vparal)
    return vperp
end
function perpendicular_velocity(
    gyrationvelocity::AbstractVector,
    b_vec::AbstractVector,
    vparal::Real,
)
    return gyrationvelocity - vparal * b_vec
end


"""
    magneticmoment(perpendicular_velocity, mass, magneticfieldstrength)
Return the magnetic moment of a chargred particle in a magnetic field.
"""
function magneticmoment(perpendicular_velocity, mass, magneticfieldstrength)
    0.5mass * perpendicular_velocity^2 / magneticfieldstrength
end


"""
Return the characteristic length of a field.

# Methods
    characteristicfieldlength(fieldstrength, fieldstrengthgradient)
    characteristicfieldlength(position, itpvec)
"""
function characteristicfieldlength(
    fieldstrength::Real,
    fieldstrengthgradient::AbstractVector
)
    return fieldstrength / norm(fieldstrengthgradient)
end
function characteristicfieldlength(
    position::AbstractVector,
    itpvec::Vector{<:Any}
)
    fieldstrength = norm([itpvec[i](position...) for i in 1:3])
    grad = ∇(position, itpvec)
    return characteristicfieldlength(fieldstrength, grad)
end
function characteristicfieldlength(
    position::AbstractVector,
    field::Function,
)
    fieldstrength = norm(field(position...))
    grad = ∇(position, field)
    return characteristicfieldlength(fieldstrength, grad)
end


"""
Calculates the ratio between the Larmor radius and the characteristic field
length of the magnetic field.

# Methods
    scalesratio(R, itpvec, ∇B, params)
    scalesratio(pos, vel, itpvec, ∇B, params)
    scalesratio(B, ∇B, vperp, q, m)
    scalesratio(statevector, time, params, switch)
    scalesratio(R, _, params)

# Arguments
- `R`, the guiding centre location.
- `pos`, the 3 component particle position.
- `vel`, the 3 component particle velocity.
- `itpvec`, the magnetic field, as a vector of interpolation objects for each
    component
- `∇B`, the gradient of the magnetic field, as a vector of components.
- `params`, `NamedTuple` containing the parameters charge, mass and magnetic
    moment.
- `statevector`, the state vector of the particle.
- `time`, the time of the state vector.
- `switch`, an integer value that determines the type of calculation to be
    performed. 1 for guiding centre approximation, 2 for full orbit.
"""
function scalesratio(
    R::AbstractVector,
    itpvec::Vector{<:Any},
    ∇B::AbstractVector,
    params::Any
)
    μ = params.magneticmoment
    m = params.mass
    B = norm([itp(R...) for itp in itpvec[1:3]])
    vperp = perpendicular_velocity(μ, m, B)
    return scalesratio(B, ∇B, vperp, params.charge, params.mass)
end
function scalesratio(
    pos::AbstractVector,
    vel::AbstractVector,
    itpvec::Vector{<:AbstractInterpolation},
    ∇B::AbstractVector,
    params::Any
)
    Bvec = [itp(pos...) for itp in itpvec[1:3]]
    Evec = [itp(pos...) for itp in itpvec[4:6]]
    B = norm(Bvec)
    vperp = norm(perpendicular_velocity(vel, Bvec, Evec))
    return scalesratio(B, ∇B, vperp, params.charge, params.mass)
end
function scalesratio(
    B::Real,
    ∇B::AbstractVector,
    vperp::Real,
    q::Real,
    m::Real
)
    L_B = characteristicfieldlength(B, ∇B)
    r_L = larmorradius(m, vperp, q, B)
    return r_L / L_B
end
function scalesratio(
    statevector::AbstractVector,
    time::Real,
    params::NamedTuple,
    switch::Int
)
    pos = statevector[1:3]
    ∇B = ∇(pos, params.fields)
    if switch == 1
        scalesratio(pos, params.fields, ∇B, params)
    elseif switch == 2
        vel = statevector[4:6]
        scalesratio(pos, vel, params.fields, ∇B, params)
    else
        error("Unknown switch value.")
    end
end
function scalesratio(R, _, params)
    emfields = params.electromagneticfield
    B = norm(emfields(R...)[2])
    vperp = perpendicular_velocity(params.magneticmoment, params.mass, B)
    ∇B = ∇(R, (R...) -> emfields(R...)[2])
    scalesratio(B, ∇B, vperp, params.charge, params.mass)
end


"""
Return the kinetic energy of a particle.

# Methods
    kineticenergy(velocity::Any, mass::Any)
    kineticenergy(velocity::Vector, mass::Any)
"""
function kineticenergy(velocity::Real, mass::Real)
    0.5 * mass * velocity^2
end
function kineticenergy(velocity::Vector, mass::Real)
    kineticenergy(norm(velocity), mass)
end


"""
    kineticenergy(
        parallel_velocity::Real,
        vperp            ::Real,
        exbdrift         ::AbstractVector,
        ∇Bdrift          ::AbstractVector,
        curvaturedrift   ::AbstractVector,
        polarisationdrift::AbstractVector,
        mass             ::Real,
        )
Return the kinetic energy of charged particle given the particle mass,
its drifts, and `parallel_velocity` and perpendicular (`vperp`)
velocity of its *guiding centre*.
"""
function kineticenergy(
    parallel_velocity::Real,
    vperp::Real,
    exbdrift::AbstractVector,
    ∇Bdrift::AbstractVector,
    curvaturedrift::AbstractVector,
    polarisationdrift::AbstractVector,
    mass::Real,
)
    vdrift = norm(exbdrift + ∇Bdrift + curvaturedrift + polarisationdrift)
    return 0.5 * mass * (parallel_velocity^2 + vperp^2 + vdrift^2)
end


"""
    kineticenergy(
         parallel_velocity::Real,
         vperp            ::Real,
         exbdrift         ::AbstractVector,
         mass             ::Real,
    )
Return the kinetic energy of charged particle given the particle mass,
its E cross B drift, and `parallel_velocity` and perpendicular (`vperp`)
velocity of its *guiding centre*. Drifts other than the E cross B drift
are neglected.
"""
function kineticenergy(
    parallel_velocity::Real,
    vperp::Real,
    exbdrift::AbstractVector,
    mass::Real,
)
    velsquared = parallel_velocity^2 + vperp^2 + norm(exbdrift)^2
    return 0.5 * mass * velsquared
end



"""
Calculate the E cross B drift in a given electromagnetic field. 

# Methods
    exbdrift(
        magneticfield::AbstractVector,
        electricfield::AbstractVector
    )
    exbdrift(
        magneticfield_direction::AbstractVector,
        electricfield::AbstractVector,
        inverted_magneticfieldstrength::Real,
    )
"""
function exbdrift(
    b::AbstractVector,
    E_vec::AbstractVector,
    B_inv::Real,
)
    return (E_vec × b) * B_inv
end
function exbdrift(
    magneticfield::AbstractVector,
    electricfield::AbstractVector
)
    B = norm(magneticfield) # Magnetic field strength
    B_inv = 1 / B
    b̂ = magneticfield * B_inv     # Magnetic field direction (unit vector)
    return exbdrift(b̂, electricfield, B_inv)
end


"""
    gradbdrift(
        b    ::AbstractVector,
        ∇B   ::AbstractVector,
        μ    ::Real,
        B_inv::Real,
        q_inv::Real,
        )
Calculate the ∇B-drift given the inverse of the magnetic field strength `B_inv`,
the magnetic field direction `b` (a unit vector), the gradient of the magnetic
field strength `∇B`, inverse of particle charge `q_inv` and magnetic moment `μ`.
"""
function gradbdrift(
    b::AbstractVector, # The direction of the magnetic field
    ∇B::AbstractVector, # The gradient of the magnetic field strength
    μ::Real, # the magnetic moment of the particle
    B_inv::Real, # The inverse of the magnetic field strength
    q_inv::Real  # The inverse of the charge of the particle
)
    return q_inv * B_inv * μ * (b × ∇B)
end


"""
    curvaturedrift(
        b     ::AbstractVector,  # The direction of the magnetic fielda
        dbdt  ::AbstractMatrix,# Material derivative of the magn. field direction
        vparal::Real, # Particle velocity component parallel to the magnetic field
        B_inv ::Real, # The inverse of the magnetic field strength
        q_inv ::Real, # The inverse of the charge of the particle
        mass  ::Real  # The particle mass
        )
Calculate the curvature drift of a charged particle in an electromagnetic field.
"""
function curvaturedrift(
    b::AbstractVector,  # The direction of the magnetic fielda
    dbdt::AbstractVector, # Material derivative of the magn. field direction
    vparal::Real, # Particle velocity component parallel to the magnetic field
    B_inv::Real, # The inverse of the magnetic field strength
    q_inv::Real, # The inverse of the charge of the particle
    mass::Real  # The particle mass
)
    return q_inv * B_inv * mass * b × (vparal * dbdt)
end


"""
    polarisationdrift(
        b     ::AbstractVector,  # The direction of the magnetic fielda
        dExBdt::AbstractMatrix,# Material derivative of the E cross B drift
        B_inv ::Real, # The inverse of the magnetic field strength
        q_inv ::Real, # The inverse of the charge of the particle
        mass  ::Real  # The particle mass
        )
Calculate the polarisation drift of a charged particle in an electromagnetic
field.
"""
function polarisationdrift(
    b::AbstractVector,  # The direction of the magnetic fielda
    dExBdt::AbstractVector,# Material derivative of the E cross B drift
    B_inv::Real, # The inverse of the magnetic field strength
    q_inv::Real, # The inverse of the charge of the particle
    mass::Real  # The particle mass
)
    return q_inv * B_inv * mass * b × dExBdt
end


"""
    magneticmirror_acceleration(
        b ::AbstractVector,
        ∇B::AbstractVector,
        μ ::Real,
        m ::Real,
        )
Calculate the magnetic mirror force acting on a guiding centre particle in a
magnetic field, given by the magnetic field direction `b`, the gradient of the
magnetic field strength `∇B`, the magnetic moment `μ` and the particle mass `m`.
"""
function magneticmirror_acceleration(
    b::AbstractVector,  # The direction of the magnetic fielda
    ∇B::AbstractVector,  # The gradient of the magnetic field strength
    μ::Real, # The magnetic moment of the particle
    m::Real, # The particle mass
)
    return -μ * b ⋅ ∇B / m
end


"""
    parallel_acceleration(
        b    ::AbstractVector,
        E_vec::AbstractVector,
        q    ::Real,
        m    ::Real,
        )
Calculate the acceleration of a charged particle due to electric fields
`E_vec` parallal to the magnetic field direction `b`.
"""
function parallel_acceleration(
    b::AbstractVector,  # The direction of the magnetic fielda
    E_vec::AbstractVector,# The electric field
    q::Real, # The charge of the particle
    m::Real, # The particle mass
)
    return q * (E_vec ⋅ b) / m
end


"""
    fermi_acceleration(
        ExBdrift::AbstractVector,
        ∇Bdrift::AbstractVector,
        dbdt   ::AbstractVector,
        )
    fermi_acceleration(
        ExBdrift::AbstractVector,
        dbdt   ::AbstractVector,
        )
Calculate the acceleration of a charged particle due to an E cross B drift
along the material derivative of the magnetic field direction. This
acceleration term is associated with a curvature drift along the electric
field. It is also termed 1st order Fermi acceleration by Birn et al. 2017.
"""
function fermi_acceleration(
    ExBdrift::AbstractVector,
    ∇Bdrift::AbstractVector,
    dbdt::AbstractVector,
)
    return (ExBdrift + ∇Bdrift) ⋅ dbdt
end
function fermi_acceleration(
    ExBdrift::AbstractVector,
    dbdt::AbstractVector,
)
    return ExBdrift ⋅ dbdt
end


"""
    fieldgradients(
        R::AbstractVector,
        electromagneticfield::Any
    )
Calculate gradients of the electromagnetic field given by the function/functor
`electromagneticfield`. The gradient that are calculated are
- ∇b, ∇ExB and ∇B

The magnetic and electric field vectors are created by the call
```
E, B = electromagneticfield(R...)
```
"""
function fieldgradients(
    R::SVector{3},
    electromagneticfield::Any
)
    # Interpolate the electromagnetic field to the guiding centre position.
    E_vec, B_vec = electromagneticfield(R...)

    # Calculate the field gradients.
    J = ForwardDiff.jacobian(R) do x
        E_vec_fd, B_vec_fd = electromagneticfield(x...)
        B_fd = norm(B_vec_fd)
        b_fd = B_vec_fd / B_fd
        out = vcat(b_fd, E_vec_fd × b_fd / B_fd)
        return [out; B_fd]
    end

    ∇b = SMatrix{3,3}(
        J[1], J[2], J[3], J[8], J[9], J[10], J[15], J[16], J[17]
    )
    ∇ExB = SMatrix{3,3}(
        J[4], J[5], J[6], J[11], J[12], J[13], J[18], J[19], J[20]
    )
    ∇B = SVector{3}(J[7], J[14], J[21])

    return ∇b, ∇ExB, ∇B, B_vec, E_vec
end


"""
    drifts(R, vparal, q, m, μ, electromagneticfield)
Calculates and returns the perpendicular drifts in the guiding centre
approximation
"""
function drifts(
    R::AbstractVector,
    vparal::Real,
    q::Real,
    m::Real,
    μ::Real,
    electromagneticfield::Any,
)
    ∇b, ∇ExB, ∇B, B_vec, E_vec = fieldgradients(
        R, electromagneticfield
    )

    B = norm(B_vec)
    B_inv = 1 / B
    b = B_vec * B_inv
    q_inv = 1 / q
    ExBdrift = exbdrift(B_vec, E_vec)

    vparal_vec = vparal * b
    total_velocity = vparal_vec + ExBdrift
    dbdt = ∇b * total_velocity
    dExBdt = ∇ExB * total_velocity

    ∇Bdrift = gradbdrift(b, ∇B, μ, B_inv, q_inv)
    Rdrift = curvaturedrift(b, dbdt, vparal, B_inv, q_inv, m)
    Pdrift = polarisationdrift(b, dExBdt, B_inv, q_inv, m)

    return ExBdrift, ∇Bdrift, Rdrift, Pdrift
end


"""
    gca_drift_and_acceleration(
        ∇b         ::AbstractMatrix,
        ∇ExBdrift  ::AbstractMatrix,
        ∇B         ::AbstractVector,
        B_vec      ::AbstractVector,
        E_vec      ::AbstractVector,
        vparal     ::Real,
        q          ::Real,
        m          ::Real,
        μ          ::Real,
        )
Calculate the drifts and accelerations of a charged particle in a 
electromagnetic field. The function returns the total time derivative
of the guiding centre position and parallel velocity.
"""
function gca_drift_and_acceleration(
    ∇b::AbstractMatrix,
    ∇ExBdrift::AbstractMatrix,
    ∇B::AbstractVector,
    B_vec::AbstractVector,
    E_vec::AbstractVector,
    vparal::Real,
    q::Real,
    m::Real,
    μ::Real,
)
    q_inv = 1 / q       # Inverse of the charge -- to save computation time.
    # Maybe not necessary in Julia?
    B = norm(B_vec)   # The magnetic field strength
    B_inv = 1 / B       # Inverse of B - to replace divition with multiplication
    b = B_vec * B_inv   # An unit vector pointing in the direction of the
    #  magnetic field

    ExBdrift = (E_vec × b) * B_inv # The E cross B-drift
    # Electric field component parallel to the magnetic field
    Eparal = E_vec ⋅ b
    # Calculate ∇B-drift
    ∇Bdrift = gradbdrift(b, ∇B, μ, B_inv, q_inv)
    #∇Bdrift = q_inv*B_inv*μ*(b × ∇B)

    # Total time derivatives. Assumes ∂/∂t = 0, and neglects other
    # drifts than ExB
    vparal_vec = vparal * b
    total_velocity = vparal_vec + ExBdrift
    dbdt = ∇b * total_velocity
    dExBdt = ∇ExBdrift * total_velocity
    #dbdt = vparal * (∇b * b) + ∇b*ExBdrift
    #dExBdt = vparal * (∇ExB * b) + ∇ExB*ExBdrift

    # Compute the perpendicular velocity
    dRperpdt =
        ExBdrift + ∇Bdrift + q_inv * B_inv * m * b × (vparal * dbdt + dExBdt)

    # Compute the acceleration along the magnetic field lines
    # With correction proposed by Birn et al., 2004:
    dvparaldt = (q * Eparal - μ * b ⋅ ∇B) / m + (ExBdrift + ∇Bdrift) ⋅ dbdt

    # Compute the velocity
    #dRdt = vparal*b + Rperp
    dRdt = vparal_vec + dRperpdt
    return [dRdt; dvparaldt]
end


"""
    get_guidingcentre!(
        u            ::Vector{<:Real},
        mu           ::Real,
        pos          ::Vector{<:Real},
        vel          ::Vector{<:Real},
        magneticfield::Vector{<:Real},
        electricfield::Vector{<:Real},
        charge       ::Real,
        mass         ::Real
        )
Calculates and returns the guiding-centre position, parallel velocity, and
 magnetic moment of a charged particle in a electromagnetic field.
"""
function get_guidingcentre!(
    u::Vector{<:Real},
    pos::Vector{<:Real},
    vel::Vector{<:Real},
    magneticfield::SVector{3},
    electricfield::SVector{3},
    charge::Real,
    mass::Real,
)
    B, b_vec, B_inv = strength_direction_invers(magneticfield)
    vel_in_E_frame = vel - exbdrift(b_vec, electricfield, B_inv)
    # Calculate the guiding centre posistion
    R = pos + mass / (charge * B) * (vel_in_E_frame × b_vec)
    # Calculate the velocity parallell to the magnetic field -- vparal
    # and the magnetic moment -- mu
    vparal = vel ⋅ b_vec
    vperp = perpendicular_velocity(vel_in_E_frame, b_vec, vparal)
    mu = magneticmoment(norm(vperp), mass, B)
    # In place return
    u[1:3] = R
    u[4] = vparal
    return mu
end


"""
    get_fullorbit(
        magneticfield::Vector{<:Real},
        electricfield::Vector{<:Real},
        R            ::Vector{<:Real},
        vparal       ::Real,
        μ            ::Real,
        mass         ::Real,
        charge       ::Real,
        phaseangle   ::Real,
        )
Get the position and velocity from the guiding centre position `R`, parallel
velocity `vparal`, and magnetic moment `μ` of a charged particle with `mass`
and `charge` in an electromagnetic field.

Requires an arbitrary `phaseangle`
of the gyromotion because we are going from 5 parameters to 6. The phase angle
is with respect to vector perpendicular to the magnetic field and guiding
centre position vector.
"""
function get_fullorbit(
    magneticfield::SVector{3},
    electricfield::SVector{3},
    R::SVector{3}, # Guiding centre position
    vparal::Real,           # Velocity parallel to the magnetic field
    μ::Real,           # Magnetic moment of particle
    charge::Real,
    mass::Real,
    phaseangle::Real,           # Arbitrary phase angle of gyration.
)
    B = norm(magneticfield) # Magnetic field strength
    b = magneticfield / B   # Magnetic field direction (unit vector)
    v_exb = exbdrift(magneticfield, electricfield) # get E cross B drift,
    # magnetic field strength and magnetic field direction
    vperp = perpendicular_velocity(μ, mass, B)
    r_L = larmorradius(mass, vperp, charge, B)
    e₁ = (R × b) / norm(R)
    e₂ = e₁ × b
    sθ, cθ = sincos(phaseangle)
    position = R + r_L * (cθ * e₁ + sθ * e₂)
    velocity = vparal * b + v_exb + vperp * (cθ * e₂ - sθ * e₁)
    return [position; velocity]
end


"""
    get_fullorbit!(
        u            ::Vector{<:Real},
        magneticfield::Vector{<:Real},
        electricfield::Vector{<:Real},
        R            ::Vector{<:Real},
        vparal       ::Real,
        μ            ::Real,
        mass         ::Real,
        charge       ::Real,
        phaseangle   ::Real,
        )
Get the position and velocity from the guiding centre position `R`, parallel
velocity `vparal`, and magnetic moment `μ` of a charged particle with `mass`
and `charge` in an electromagnetic field.

Requires an arbitrary `phaseangle`
of the gyromotion because we are going from 5 parameters to 6. The phase angle
is with respect to vector perpendicular to the magnetic field and guiding
centre position vector.
"""
function get_fullorbit!(
    u::Vector{<:Real},
    magneticfield::Vector{<:Real},
    electricfield::Vector{<:Real},
    R::Vector{<:Real}, # Guiding centre position
    vparal::Real,           # Velocity parallel to the magnetic field
    μ::Real,           # Magnetic moment of particle
    charge::Real,
    mass::Real,
    phaseangle::Real,           # Arbitrary phase angle of gyration.
)
    B = norm(magneticfield) # Magnetic field strength
    b = magneticfield / B   # Magnetic field direction (unit vector)
    v_exb = exbdrift(magneticfield, electricfield) # get E cross B drift,
    # magnetic field strength and magnetic field direction
    vperp = perpendicular_velocity(μ, mass, B)
    r_L = larmorradius(mass, vperp, charge, B)
    e₁ = (R × b) / norm(R)
    e₂ = e₁ × b
    #e₂ = b × e₁
    sθ, cθ = sincos(phaseangle)
    position = R + r_L * (cθ * e₁ + sθ * e₂)
    #velocity = vparal*b + v_exb + vperp*(cθ*e₂ - sθ*e₁) !OLD!
    velocity = vparal * b + v_exb + sign(charge) * vperp * (cθ * e₂ + sθ * e₁)
    # In-place return
    u[:] = [position; velocity]
end

"""
    cosineof_pitchangle(
        magneticfield    ::AbstractVector,
        parallel_velocity::Real,
        mass             ::Real,
        magneticmoment   ::Real,
        )
Calculates the cosine of the pitch angle of a charge particle in a magnetic
field, given the particle's parallel velocity, mass and magnetic moment.
"""
function cosineof_pitchangle(
    magneticfield::AbstractVector,
    parallel_velocity::Real,
    mass::Real,
    magneticmoment::Real,
)
    B = norm(magneticfield)
    return sqrt(1 / (2B * magneticmoment / (mass * parallel_velocity^2) + 1))
end


"""
    lorentzfactor(speed::Real)
The relativistic Lorentz factor in SI units.
"""
function lorentzfactor(speed::Real)
    return 1 / √(1 - speed^2 * csqrdinv)
end


"""
    kineticspeed(kineticenergy::Real, mass::Real)
Calculate the speed of a particle from its `kineticenergy` and `mass`.
"""
function kineticspeed(kineticenergy::Real, mass::Real)
    return √(2kineticenergy / mass)
end


"""
    mfp_close_cc_collisions_thermaltarget(v_p, q_p, q_t, m_p, m_t, T_t, n_t)
Calculate the mean free path of a particle colliding with cold target
in the close Coulomb collision approximation.
The parameters are:
- `v_p`: the particle speed
- `q_p`: the particle charge
- `q_t`: the target charge
- `m_p`: the particle mass
- `m_t`: the target mass
- `T_t`: the target temperature
- `n_t`: the target number density
The result is the mean free path in meters.
"""
function mfp_close_cc_thermaltarget(v, q, q_t, m, m_t, T_t, n_t)
    mu = m * m_t / (m + m_t) # reduced mass
    α = 16π * TraceParticles.ϵ_0^2 * m * mu / (n_t * q^2 * q_t^2)
    return α * v^4
end

"""
    mfp_distant_cc_thermaltarget(v, q, q_t, m, m_t, T_t, n_t; coulomb_logarithm=20)
Calculate the mean free path of a particle colliding with a cold target
in the distant Coulomb collision approximation.
"""
function mfp_distant_cc_thermaltarget(v, q, q_t, m, m_t, T_t, n_t; coulomb_logarithm=20)
    return mfp_close_cc_thermaltarget(v, q, q_t, m, m_t, T_t, n_t) / coulomb_logarithm
end

"""
    collisionaltime(meanfreepath, speed)
Calculate the collisional time of a particle given its mean free path.
"""
function collisionaltime(meanfreepath, speed)
    return meanfreepath / speed
end

"""
    spitzercollisionaltime(q, m, T, n; coulomb_logarithm=20)
Calculate the collisional time in SI-units of a particle in a plasma due to
Coulomb collisions.
Source: Somov 2012 "Plasma Astrophysics, Part I, section 8.3.1, page 158.
"""
function spitzercollisionaltime(q, m, T, n; coulomb_logarithm=20)
    # From Somov 2012 "Plasma Astrophysics, Part I, section 8.3.1, page 158
    numerator = 16π * TraceParticles.ϵ_0^2 * m^2 * (3 * TraceParticles.k_B * T / m)^(3 / 2)
    denominator = (n * q^4 * 8 * coulomb_logarithm) * 0.714
    return numerator / denominator
end

"""
    spitzercollisionaltime_cgs(q, m, t, n; coulomb_logarithm=20)
calculate the collisional time in cgs-units of a particle in a plasma due to
coulomb collisions.
source: Somov 2012 "plasma astrophysics, part i, section 8.3.1, page 158.
"""
function spitzercollisionaltime_cgs(q, m, T, n; coulomb_logarithm=20)
    # From Somov 2012 "Plasma Astrophysics, Part I, section 8.3.1, page 158
    # in CGS units
    numerator = m^2 * (3 * TraceParticles.k_B_cgs * T / m)^(3 / 2)
    denominator = (π * n * e^4 * 8 * coulomb_logarithm) * 0.714
    return numerator / denominator
end

"""
    spitzercollisionaltime_chen(q, m, t, n; coulomb_logarithm=20)
calculate the collisional time in SI-units of a particle in a plasma due to
coulomb collisions.
source: Chen 2016 "Introduction to Plasma Physics and Controlled Fusion,
"""
function spitzercollisionaltime_chen(q, m, T, n; coulomb_logarithm=20)
    # From Somov 2012 "Plasma Astrophysics, Part I, section 8.3.1, page 158
    numerator = 16π * TraceParticles.ϵ_0^2 * m^2 * (TraceParticles.k_B * T / m)^(3 / 2)
    denominator = n * q^4 * coulomb_logarithm
    return numerator / denominator
end

"""
    localcollisiontime(particledata, densityitp, temperatureitp)
Calculate the Coulomb collision time for an ensemble of particles, where
the local density and temperature at the particle's initial position is used.
"""
function localcollisiontime(particledata, densityitp, temperatureitp)
    x0 = h5_getdataset(particledata, "x0")
    y0 = h5_getdataset(particledata, "y0")
    z0 = h5_getdataset(particledata, "z0")
    charge = h5_getdataset(particledata, "charge")
    mass = h5_getdataset(particledata, "mass")
    n_itp = load_object(densityitp)
    T_itp = load_object(temperatureitp)
    n = [n_itp(x0[i], y0[i], z0[i]) / TraceParticles.m_p for i in eachindex(x0)]
    T = [T_itp(x0[i], y0[i], z0[i]) for i in eachindex(x0)]
    # Assuming the like-particle collisions are the most frequent
    return spitzercollisionaltime.(charge, mass, T, n)
end

"""
    dreicerfield_cgs(n, T; coulomb_logarithm=20)
Calculate the Dreicer field in CGS units.
"""
function dreicerfield_cgs(n, T; coulomb_logarithm=20)
    nominator = 4π * TraceParticles.e_cgs^3 * coulomb_logarithm * n
    denominator = TraceParticles.k_B_cgs * T
    return nominator / denominator
end

"""
    criticalvelocity_cgs(n, T, E, m; coulomb_logarithm=20)
Calculate the critical velocity of particles with mass `m` in CGS units.
The critical velocity is the velocity at which the Dreicer field
is equal to the electric field `E`, or the velocity at which the frictional
drag from Coulomb collisions (of thermal like-particle targets) is equal to
the acceleration from the parallel electric field.
"""
function criticalvelocity_cgs(n, T, E, m; coulomb_logarithm=20)
    E_D = dreicerfield_cgs(n, T; coulomb_logarithm=coulomb_logarithm)
    v_th = sqrt(TraceParticles.k_B_cgs * T / m)
    return v_th * sqrt(E_D / E)
end

"""
    localcriticalvelocity_cgs(
        particledata,
        densityitp,
        temperatureitp,
        electricfielditp;
        coulomb_logarithm=20
    )
Calculate the critical velocity of an ensemble of particles in CGS units,
using the local density, temperature and parallel electric field at the
particle's initial position.
"""
function localcriticalvelocity_cgs(
    particledata,
    densityitp,
    temperatureitp,
    electricfielditp;
    coulomb_logarithm=20
)
    x0 = h5_getdataset(particledata, "x0")
    y0 = h5_getdataset(particledata, "y0")
    z0 = h5_getdataset(particledata, "z0")
    mass = h5_getdataset(particledata, "mass")
    n_itp = load_object(densityitp)
    T_itp = load_object(temperatureitp)
    E_itp = load_object(electricfielditp)
    n = [n_itp(x0[i], y0[i], z0[i]) * 1e-3 / TraceParticles.m_p_cgs for i in eachindex(x0)]
    T = [T_itp(x0[i], y0[i], z0[i]) for i in eachindex(x0)]
    E = [abs(E_itp(1e2x0[i], 1e2y0[i], 1e2z0[i])) for i in eachindex(x0)]
    return criticalvelocity_cgs.(
        n, T, E, mass;
        coulomb_logarithm=coulomb_logarithm
    )
end

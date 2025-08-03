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
Return the Larmor radius of a charged particle.
# Methods
    larmorradius(mass, vperp, charge, magneticfieldstrength)
    larmorradius(vel, mass, charge, electricfield, magneticfield)
    # Methods with interpolation objects
    larmorradius(R, μ, q, m, emfields) # for guiding centre approximation
    larmorradius(position, velocity, mass, charge, emfields) # direct solution
"""
function larmorradius(
    mass::Real,
    vperp::Real,
    charge::Real,
    magneticfieldstrength::Real
)::Real
    return abs(mass * vperp / (charge * magneticfieldstrength))
end
function larmorradius(
    vel::AbstractVector,
    mass::Real,
    charge::Real,
    electricfield::AbstractVector,
    magneticfield::AbstractVector,
)
    vperp = norm(perpendicular_velocity(
        vel, magneticfield, electricfield
    ))
    B = norm(magneticfield)
    return larmorradius(mass, vperp, charge, B)
end
function larmorradius(
    R::AbstractVector,
    μ::Real,
    q::Real,
    m::Real,
    emfields
)
    B = norm(emfields(R[1], R[2], R[3])[2])
    vperp = perpendicular_velocity(μ, m, B)
    return larmorradius(m, vperp, q, B)
end
function larmorradius(
    position::AbstractVector,
    velocity::AbstractVector,
    mass::Real,
    charge::Real,
    emfields,
)
    E, B = emfields(position[1], position[2], position[3])
    return larmorradius(velocity, mass, charge, E, B)
end

"""
    gyrofrequency(mass, charge, magneticfieldstrength)
Return the gyrofrequency of a charged particle.
"""
function gyrofrequency(
    mass::Real,
    charge::Real,
    magneticfieldstrength::Real,
)::Real
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
    gyrationvelocity::AbstractVector,
    b_vec::AbstractVector,
    vparal::Real,
)
    return gyrationvelocity - vparal * b_vec
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

"""
Return the magnetic moment of a chargred particle in a magnetic field.
# Methods
    magneticmoment(perpendicular_velocity, mass, magneticfieldstrength)
    magneticmoment(vel, mass, electricfield, magneticfield)
    # interpolation object for direct of the Lorentz equation.
    magneticmoment(position, velocity, mass, emfields)
"""
function magneticmoment(
    perpendicular_velocity::Real,
    mass::Real,
    magneticfieldstrength::Real
)
    0.5mass * perpendicular_velocity^2 / magneticfieldstrength
end
function magneticmoment(
    vel::AbstractVector,
    mass::Real,
    electricfield::AbstractVector,
    magneticfield::AbstractVector,
)
    vperp = norm(perpendicular_velocity(
        vel,
        magneticfield,
        electricfield
    ))
    B = norm(magneticfield)
    return magneticmoment(vperp, mass, B)
end
function magneticmoment(
    position::AbstractVector,
    velocity::AbstractVector,
    mass::Real,
    emfields
)
    electricfield, magneticfield = emfields(
        position[1], position[2], position[3]
    )
    return magneticmoment(velocity, mass, electricfield, magneticfield)
end

"""
Return the characteristic length of a field.
# Methods
    characteristicfieldlength(fieldstrength, fieldstrengthgradient)
    characteristicfieldlength(position, emfields)
"""
function characteristicfieldlength(
    fieldstrength::Real,
    fieldstrengthgradient::AbstractVector
)
    return fieldstrength / norm(fieldstrengthgradient)
end
function characteristicfieldlength(
    position::AbstractVector,
    emfields,
)
    fieldstrength = norm(emfields(position[1], position[2], position[3])[2])
    grad = ForwardDiff.gradient(position) do x
        norm(emfields(x[1], x[2], x[3])[2])
    end
    return characteristicfieldlength(fieldstrength, grad)
end


"""
Calculates the ratio between the Larmor radius and the characteristic field
length of the magnetic field.

# Methods
    scalesratio(B, ∇B, vperp, q, m)
    # Using interpolation object
    scalesratio(R, mass, charge, magneticmoment, emfields)
    scalesratio(pos, vel, mass, charge, emfields)

# Arguments
- `R`, the guiding centre location.
- `pos`, the 3 component particle position.
- `vel`, the 3 component particle velocity.
- `∇B`, the magnitude of the gradient of the magnetic field
- `B`, the magnetic field strength.
- `vperp`, the particle velocity component perpendicular to the magnetic field.
- `q`, the particle charge.
- `m`, the particle mass.
- `emfields`, a functor that returns a tuple of the electric and magnetic field
  vectors at a given position, in the form
"""
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
    R::AbstractVector,
    mass::Real,
    charge::Real,
    magneticmoment::Real,
    emfields,
)
    B = norm(emfields(R[1], R[2], R[3])[2])
    vperp = perpendicular_velocity(magneticmoment, mass, B)
    ∇B = ForwardDiff.gradient(R) do x
        norm(emfields(x[1], x[2], x[3])[2])
    end
    scalesratio(B, ∇B, vperp, charge, mass)
end
function scalesratio(
    pos::AbstractVector,
    vel::AbstractVector,
    mass::Real,
    charge::Real,
    emfields
)
    electricfield, magneticfield = emfields(pos[1], pos[2], pos[3])
    vperp = norm(perpendicular_velocity(
        vel,
        magneticfield,
        electricfield
    ))
    B = norm(magneticfield)
    ∇B = ForwardDiff.gradient(pos) do x
        norm(emfields(x[1], x[2], x[3])[2])
    end
    return scalesratio(B, ∇B, vperp, charge, mass)
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
    exbdrift(pos, emfields)
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
function exbdrift(
    pos::AbstractVector,
    emfields,
)
    electricfield, magneticfield = emfields(pos[1], pos[2], pos[3])
    return exbdrift(magneticfield, electricfield)
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
        x::Real,
        y::Real,
        z::Real,
        time::Real,
        electromagneticfield::Any
    )
Calculate gradients of the electromagnetic field given by the function/functor
`electromagneticfield`. The gradient that are calculated are
- ∇b, ∇ExB, ∇B, ∂b, ∂ExB, ∂B, where '∇' is the 3D spatial gradient and '∂' is
the time derivative.

The magnetic and electric field vectors are created by the call
```
E, B = electromagneticfield(x, y, z, time)
```

Return the gradients as `StaticArrays` if `electromagneticfield` returns
`StaticArrays`.
"""
function fieldgradients(
    x::Real,
    y::Real,
    z::Real,
    time::Real,
    electromagneticfield::Any
)
    # Interpolate the electromagnetic field to the guiding centre position.
    E_vec, B_vec = electromagneticfield(x, y, z, time)

    spacetime = SVector{4}(x, y, z, time)

    # Calculate the field gradients.
    J = ForwardDiff.jacobian(spacetime) do x
        E_vec_fd, B_vec_fd = electromagneticfield(x[1], x[2], x[3], x[4])
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

    # The time derivatives
    ∂b = SVector{3}(J[22], J[23], J[24])
    ∂ExB = SVector{3}(J[25], J[26], J[27])
    ∂B = J[28]

    return ∇b, ∇ExB, ∇B, B_vec, E_vec, ∂b, ∂ExB, ∂B
end

"""
    drifts(R, vparal, q, m, μ, electromagneticfield)
Calculates and returns the perpendicular drifts in the guiding centre
approximation
"""
function drifts(
    x::Real,
    y::Real,
    z::Real,
    t::Real,
    vparal::Real,
    q::Real,
    m::Real,
    μ::Real,
    electromagneticfield::Any,
)
    ∇b, ∇ExB, ∇B, B_vec, E_vec, ∂b, ∂ExB, _ = fieldgradients(
        x, y, z, t, electromagneticfield
    )
    B = norm(B_vec)
    B_inv = 1 / B
    b = B_vec * B_inv
    q_inv = 1 / q
    ExBdrift = exbdrift(B_vec, E_vec)

    vparal_vec = vparal * b
    total_velocity = vparal_vec + ExBdrift
    dbdt = ∇b * total_velocity + ∂b
    dExBdt = ∇ExB * total_velocity + ∂ExB

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
electromagnetic fiel using the guiding centre approximation. The function
returns the total time derivative of the guiding centre position and parallel
velocity in the most computationally efficient way.
"""
function gca_drift_and_acceleration(
    ∇b::AbstractMatrix,
    ∇ExBdrift::AbstractMatrix,
    ∇B::AbstractVector,
    B_vec::AbstractVector,
    E_vec::AbstractVector,
    ∂b::AbstractVector,
    ∂ExB::AbstractVector,
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

    # Total time derivatives. Neglects other drifts than ExB
    vparal_vec = vparal * b
    total_velocity = vparal_vec + ExBdrift
    dbdt = ∇b * total_velocity + ∂b
    dExBdt = ∇ExBdrift * total_velocity + ∂ExB
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
    get_guidingcentre(
        pos          ::AbstractVector,
        vel          ::AbstractVector,
        magneticfield::AbstractVector,
        electricfield::AbstractVector,
        charge       ::Real,
        mass         ::Real
        )
Calculates and returns the guiding-centre position, parallel velocity, and
 magnetic moment of a charged particle in a electromagnetic field.
"""
function get_guidingcentre(
    pos::AbstractVector,
    vel::AbstractVector,
    magneticfield::AbstractVector,
    electricfield::AbstractVector,
    charge::Real,
    mass::Real,
)
    B = norm(magneticfield) # Magnetic field strength
    B_inv = 1 / B # Inverse of the magnetic field strength
    b_vec = magneticfield * B_inv # Magnetic field direction (unit vector)
    vel_in_E_frame = vel - exbdrift(b_vec, electricfield, B_inv)
    # Calculate the guiding centre posistion
    R = pos + mass / (charge * B) * (vel_in_E_frame × b_vec)
    # Calculate the velocity parallell to the magnetic field -- vparal
    # and the magnetic moment -- mu
    vparal = vel ⋅ b_vec
    vperp = perpendicular_velocity(vel_in_E_frame, b_vec, vparal)
    mu = magneticmoment(norm(vperp), mass, B)
    # Out-of-plece return
    return R, vparal, mu
end
function get_guidingcentre(
    pos::AbstractVector,
    vel::AbstractVector,
    electromagneticfield,
    args...
)
    electricfield, magneticfield = electromagneticfield(pos[1], pos[2], pos[3])
    return get_guidingcentre(pos, vel, magneticfield, electricfield, args...)
end

"""
    get_guidingcentre!(u::AbstractVector, args...)
In-place version of `get_guidingcentre`.
"""
function get_guidingcentre!(
    u::AbstractVector,
    args...
)
    R, vparal, mu = get_guidingcentre(args...)
    u[1:3] .= R
    u[4] = vparal
    return mu
end

"""
    get_fullorbit(
        magneticfield::AbstractVector,
        electricfield::AbstractVector,
        R            ::AbstractVector,
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
    magneticfield::AbstractVector,
    electricfield::AbstractVector,
    R::AbstractVector, # Guiding centre position
    vparal::Real, # Velocity parallel to the magnetic field
    μ::Real, # Magnetic moment of particle
    charge::Real,
    mass::Real,
    phaseangle::Real, # Arbitrary phase angle of gyration.
)
    B = norm(magneticfield) # Magnetic field strength
    b = magneticfield / B   # Magnetic field direction (unit vector)
    v_exb = exbdrift(magneticfield, electricfield) # get E cross B drift,
    # magnetic field strength and magnetic field direction
    vperp = perpendicular_velocity(μ, mass, B)
    r_L = larmorradius(mass, vperp, charge, B)
    e₁ = b × R
    e₁ = e₁ / norm(e₁) # Normalize the vector
    e₂ = b × e₁
    e₂ = e₂ / norm(e₂) # Normalize the vector
    sθ, cθ = sincos(phaseangle)
    position = R + r_L * (sθ * e₁ + cθ * e₂)
    velocity = sign(charge) * vperp * (cθ * e₁ - sθ * e₂) + vparal * b + v_exb
    return [position; velocity]
end
function get_fullorbit(
    electromagneticfield,
    R::AbstractVector, # Guiding centre position
    args...
)
    electricfield, magneticfield = electromagneticfield(R[1], R[2], R[3])
    return get_fullorbit(magneticfield, electricfield, R, args...)
end

"""
    get_fullorbit!(u::AbstractVector, args...)
In-place version of `get_fullorbit`.
"""
function get_fullorbit!(
    u::AbstractVector,
    args...
)
    u[:] .= get_fullorbit(args...)
    return nothing
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

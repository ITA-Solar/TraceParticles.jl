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
    R::Vector{<:Real},
    μ::Real,
    q::Real,
    m::Real,
    itpvec::Vector{<:AbstractInterpolation}
)
    B = norm([itp(R...) for itp in itpvec])
    vperp = perpendicular_velocity(μ, m, B)
    return larmorradius(m, vperp, q, B)
end
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
"""
function perpendicular_velocity(magnetic_moment, mass, magneticfieldstrength::Real)
    sqrt(2magnetic_moment * magneticfieldstrength / mass)
end
function perpendicular_velocity(
    vel::Vector{<:Real},
    magneticfield::Vector{<:Real},
    electricfield::Vector{<:Real},
)
    B = norm(magneticfield) # Magnetic field strength
    b_vec = magneticfield / B   # Magnetic field direction (unit vector)
    ExBdrift = exbdrift(magneticfield, electricfield) # get E cross B drift,
    vel_in_E_frame = vel - ExBdrift
    vparal = vel ⋅ b_vec
    # Calculate mangetic moment -- mu
    vperp = vel_in_E_frame - vparal * b_vec
    return vperp
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
    fieldstrengthgradient::Vector{<:Real}
)
    return fieldstrength / norm(fieldstrengthgradient)
end
function characteristicfieldlength(
    position::Vector{<:Real},
    itpvec::Vector{<:AbstractInterpolation}
)
    fieldstrength = norm([itpvec[i](position...) for i in 1:3])
    grad = ∇(position, itpvec)
    return characteristicfieldlength(fieldstrength, grad)
end


"""
Calculates the ratio between the Larmor radius and the characteristic field
length of the magnetic field.

# Methods
    scalesratio(R, itpvec, ∇B, params)
    scalesratio(pos, vel, itpvec, ∇B, params)
    scalesratio(statevector, time, params, switch)

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
    R::Vector{<:Real},
    itpvec::Vector{Any},
    ∇B::Vector{<:Real},
    params::Any
)
    μ = params.magneticmoment
    q = params.charge
    m = params.mass
    B = norm([itp(R...) for itp in itpvec[1:3]])
    L_B = characteristicfieldlength(B, ∇B)
    vperp = perpendicular_velocity(μ, m, B)
    r_L = larmorradius(m, vperp, q, B)
    return r_L / L_B
end
function scalesratio(
    pos::Vector{<:Real},
    vel::Vector{<:Real},
    itpvec::Vector{<:AbstractInterpolation},
    ∇B::Vector{<:Real},
    params::Any
)
    q = params.charge
    m = params.mass
    Bvec = [itp(pos...) for itp in itpvec[1:3]]
    Evec = [itp(pos...) for itp in itpvec[4:6]]
    B = norm(Bvec)
    L_B = characteristicfieldlength(B, ∇B)
    vperp = norm(perpendicular_velocity(vel, Bvec, Evec))
    r_L = larmorradius(m, vperp, q, B)
    return r_L / L_B
end
function scalesratio(statevector, time, params, switch)
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
function scalesratio_2Dxz(u, time, params)
    R = [u[1], u[3]]
    itpvec = params.fields[1:3]
    ∇B = ∇(R, itpvec)
    scalesratio(R, itpvec, ∇B, params)
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
        exbdrift         ::Vector{<:Real},
        ∇Bdrift          ::Vector{<:Real},
        curvaturedrift   ::Vector{<:Real},
        polarisationdrift::Vector{<:Real},
        mass             ::Real,
        )
Return the kinetic energy of charged particle given the particle mass,
its drifts, and `parallel_velocity` and perpendicular (`vperp`)
velocity of its *guiding centre*.
"""
function kineticenergy(
    parallel_velocity::Real,
    vperp::Real,
    exbdrift::Vector{<:Real},
    ∇Bdrift::Vector{<:Real},
    curvaturedrift::Vector{<:Real},
    polarisationdrift::Vector{<:Real},
    mass::Real,
)
    vdrift = norm(exbdrift + ∇Bdrift + curvaturedrift + polarisationdrift)
    return 0.5 * mass * (parallel_velocity^2 + vperp^2 + vdrift^2)
end


"""
    kineticenergy(
         parallel_velocity::Real,
         vperp            ::Real,
         exbdrift         ::Vector{<:Real},
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
    exbdrift::Vector{<:Real},
    mass::Real,
)
    velsquared = parallel_velocity^2 + vperp^2 + norm(exbdrift)^2
    return 0.5 * mass * velsquared
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


"""
    kineticenergy(
        state::Vector{<:Real}, # Guiding centre state
        q::Real, # Particle charge
        m::Real, # Particle mass
        μ::Real, # Particle magnetic moment
        itpvec::Vector{<:AbstractInterpolation}, 
        dimensionality="3D",
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
    itpvec::Vector{<:AbstractInterpolation}, # electromagnetic field interpolators
    time::Real
    ;
    dimensionality="2Dxz"
)
    kineticenergy(GCAState(
        state,
        q,
        m,
        μ,
        itpvec;
        time=time,
        dimensionality=dimensionality,
    ))
end


"""
Calculate the E cross B drift in a given electromagnetic field. 

# Methods
    exbdrift(
        magneticfield::Vector{<:Real}, 
        electricfield::Vector{<:Real}
    )
    exbdrift(
        magneticfield_direction::Vector{<:Real},
        electricfield::Vector{<:Real},
        inverted_magneticfieldstrength::Real,
    )
"""
function exbdrift(
    b̂::Vector{<:Real},
    E_vec::Vector{<:Real},
    B_inv::Real,
)
    return (E_vec × b̂) * B_inv
end
function exbdrift(magneticfield::Vector{<:Real}, electricfield::Vector{<:Real})
    B = norm(magneticfield) # Magnetic field strength
    B_inv = 1 / B
    b̂ = magneticfield * B_inv     # Magnetic field direction (unit vector)
    return exbdrift(b̂, electricfield, B_inv)
end



"""
    gradbdrift(
        b̂    ::Vector{<:Real},
        ∇B   ::Vector{<:Real},
        B_inv::Real,
        μ    ::Real,
        q_inv::Real,
        )
Calculate the ∇B-drift given the inverse of the magnetic field strength `B_inv`,
the magnetic field direction `b̂` (a unit vector), the gradient of the magnetic
field strength `∇B`, inverse of particle charge `q_inv` and magnetic moment `μ`.
"""
function gradbdrift(
    b̂::Vector{<:Real}, # The direction of the magnetic fielda
    ∇B::Vector{<:Real}, # The gradient of the magnetic field strength
    μ::Real, # the magnetic moment of the particle
    B_inv::Real, # The inverse of the magnetic field strength
    q_inv::Real  # The inverse of the charge of the particle
)
    return q_inv * B_inv * μ * (b̂ × ∇B)
end


"""
    curvaturedrift(
        b̂     ::Vector{<:Real},  # The direction of the magnetic fielda
        db̂dt  ::Matrix{<:Real},# Material derivative of the magn. field direction
        vparal::Real, # Particle velocity component parallel to the magnetic field
        B_inv ::Real, # The inverse of the magnetic field strength
        q_inv ::Real, # The inverse of the charge of the particle
        mass  ::Real  # The particle mass
        )
Calculate the curvature drift of a charged particle in an electromagnetic field.
"""
function curvaturedrift(
    b̂::Vector{<:Real},  # The direction of the magnetic fielda
    db̂dt::Vector{<:Real},# Material derivative of the magn. field direction
    vparal::Real, # Particle velocity component parallel to the magnetic field
    B_inv::Real, # The inverse of the magnetic field strength
    q_inv::Real, # The inverse of the charge of the particle
    mass::Real  # The particle mass
)
    return q_inv * B_inv * mass * b̂ × (vparal * db̂dt)
end


"""
    polarisationdrift(
        b̂     ::Vector{<:Real},  # The direction of the magnetic fielda
        dExBdt::Matrix{<:Real},# Material derivative of the E cross B drift
        B_inv ::Real, # The inverse of the magnetic field strength
        q_inv ::Real, # The inverse of the charge of the particle
        mass  ::Real  # The particle mass
        )
Calculate the polarisation drift of a charged particle in an electromagnetic
field.
"""
function polarisationdrift(
    b̂::Vector{<:Real},  # The direction of the magnetic fielda
    dExBdt::Vector{<:Real},# Material derivative of the E cross B drift
    B_inv::Real, # The inverse of the magnetic field strength
    q_inv::Real, # The inverse of the charge of the particle
    mass::Real  # The particle mass
)
    return q_inv * B_inv * mass * b̂ × dExBdt
end


"""
    magneticmirror_acceleration(
        b̂ ::Vector{<:Real},
        ∇B::Vector{<:Real},
        μ ::Real,
        m ::Real,
        )
Calculate the magnetic mirror force acting on a guiding centre particle in a
magnetic field, given by the magnetic field direction `b̂`, the gradient of the
magnetic field strength `∇B`, the magnetic moment `μ` and the particle mass `m`.
"""
function magneticmirror_acceleration(
    b̂::Vector{<:Real},  # The direction of the magnetic fielda
    ∇B::Vector{<:Real},  # The gradient of the magnetic field strength
    μ::Real, # The magnetic moment of the particle
    m::Real, # The particle mass
)
    return -μ * b̂ ⋅ ∇B / m
end


"""
    parallel_acceleration(
        b̂    ::Vector{<:Real},
        E_vec::Vector{<:Real},
        q    ::Real,
        m    ::Real,
        )
Calculate the acceleration of a charged particle due to electric fields
`E_vec` parallal to the magnetic field direction `b̂`.
"""
function parallel_acceleration(
    b̂::Vector{<:Real},  # The direction of the magnetic fielda
    E_vec::Vector{<:Real},# The electric field
    q::Real, # The charge of the particle
    m::Real, # The particle mass
)
    return q * (E_vec ⋅ b̂) / m
end


"""
    fermi_acceleration(
        ExBdrift::Vector{<:Real},
        ∇Bdrift::Vector{<:Real},
        dbdt   ::Vector{<:Real},
        )
    fermi_acceleration(
        ExBdrift::Vector{<:Real},
        dbdt   ::Vector{<:Real},
        )
Calculate the acceleration of a charged particle due to an E cross B drift
along the material derivative of the magnetic field direction. This
acceleration term is associated with a curvature drift along the electric
field. It is also termed 1st order Fermi acceleration by Birn et al. 2017.
"""
function fermi_acceleration(
    ExBdrift::Vector{<:Real},
    ∇Bdrift::Vector{<:Real},
    dbdt::Vector{<:Real},
)
    return (ExBdrift + ∇Bdrift) ⋅ dbdt
end
function fermi_acceleration(
    ExBdrift::Vector{<:Real},
    dbdt::Vector{<:Real},
)
    return ExBdrift ⋅ dbdt
end


"""
    fieldgradients(
        R::Vector{<:Real},
        itpvec
    )
Calculate gradients of the electromagnetic field given by the function `itpvec`.
The gradient that are calculated are
- ∇b, ∇ExB and ∇B

The magnetic and electric field vectors are created by the call
```
B_vec, E_vec = emfieldatpos(R, itpvec)
```
"""
function fieldgradients(
    R::Vector{<:Real},
    itpvec
)
    # Interpolate the electromagnetic field to the guiding centre position.
    B_vec, E_vec = emfieldatpos(R, itpvec)

    # Calculate the field gradients.
    jacobian_matrix = ForwardDiff.jacobian(R) do x
        if typeof(itpvec) <: Vector
            vec = [itp(x...) for itp in itpvec]
        else
            vec = itpvec(x...)
        end
        B_vec_fd = vec[1:3]
        E_vec_fd = vec[4:6]
        B_fd = norm(B_vec_fd)
        b_fd = B_vec_fd / B_fd
        return [b_fd; E_vec_fd × b_fd / B_fd; B_fd]
    end
    ∇b = jacobian_matrix[1:3, :]
    ∇ExB = jacobian_matrix[4:6, :]
    ∇B = jacobian_matrix[7, :]

    return ∇b, ∇ExB, ∇B, B_vec, E_vec
end


"""
    drifts(R, vparal, q, m, μ, itpvec)
Calculates and returns the perpendicular drifts in the guiding centre
approximation
"""
function drifts(
    R::Vector{<:Real},
    vparal::Real,
    q::Real,
    m::Real,
    μ::Real,
    itpvec
)
    ∇b, ∇ExB, ∇B, B_vec, E_vec = fieldgradients(R, itpvec)

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
        ∇b         ::Matrix{<:Real},
        ∇ExBdrift  ::Matrix{<:Real},
        ∇B         ::Vector{<:Real},
        B_vec      ::Vector{<:Real},
        E_vec      ::Vector{<:Real},
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
    ∇b::Matrix{<:Real},
    ∇ExBdrift::Matrix{<:Real},
    ∇B::Vector{<:Real},
    B_vec::Vector{<:Real},
    E_vec::Vector{<:Real},
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
    dRperpdt = ExBdrift + ∇Bdrift + q_inv * B_inv * m * b × (vparal * dbdt + dExBdt)

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
function get_guidingcentre(
    pos::Vector{<:Real},
    vel::Vector{<:Real},
    magneticfield::Vector{<:Real},
    electricfield::Vector{<:Real},
    charge::Real,
    mass::Real
)

    B = norm(magneticfield) # Magnetic field strength
    B_inv = 1 / B
    b_vec = magneticfield / B_inv     # Magnetic field direction (unit vector)
    ExBdrift = exbdrift(magneticfield, electricfield)
    vel_in_E_frame = vel - ExBdrift
    # Calculate the guiding centre posistion 
    R = pos + mass / (charge * B) * (vel_in_E_frame × b_vec)
    # Calculate the velocity parallell to the magnetic field -- vparal
    vparal = vel ⋅ b_vec
    # Calculate mangetic moment -- mu
    vperp = vel_in_E_frame - vparal * b_vec
    mu = magneticmoment(norm(vperp), mass, B)
    return R, vparal, mu
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
    magneticfield::Vector{<:Real},
    electricfield::Vector{<:Real},
    charge::Real,
    mass::Real
)
    B = norm(magneticfield) # Magnetic field strength
    b_vec = magneticfield / B     # Magnetic field direction (unit vector)
    ExBdrift = exbdrift(magneticfield, electricfield)
    vel_in_E_frame = vel - ExBdrift
    # Calculate the guiding centre posistion
    R = pos + mass / (charge * B) * (vel_in_E_frame × b_vec)
    # Calculate the velocity parallell to the magnetic field -- vparal
    vparal = vel ⋅ b_vec
    # Calculate mangetic moment -- mu
    vperp = vel_in_E_frame - vparal * b_vec
    # In place return
    mu = magneticmoment(norm(vperp), mass, B)
    u[1:3] = R
    u[4] = vparal
    return mu
end


"""
    get_guidingcentre_2Dxzitp(
        pos          ::Vector{<:Real},
        vel          ::Vector{<:Real},
        emfields_itpvec::Vector{<:AbstractInterpolation},
        charge       ::Real,
        mass         ::Real
        )
Same as `get_guidingcentre` but with interpolation vector for the
electromagnetic fields.
"""
function get_guidingcentre_2Dxzitp(
    pos::Vector{<:Real},
    vel::Vector{<:Real},
    emfields_itpvec::Vector{<:AbstractInterpolation},
    charge::Real,
    mass::Real
)
    posx, poz = pos[1], pos[3]
    magneticfield = [emfields_itpvec[i](posx, poz) for i = 1:3]
    electricfield = [emfields_itpvec[i](posx, poz) for i = 4:6]
    get_guidingcentre(pos, vel, magneticfield, electricfield, charge, mass)
end


function get_gca_velocities(
    vel::Vector{<:Real},
    magneticfield::Vector{<:Real},
    electricfield::Vector{<:Real},
    mass::Real
)
    B = norm(magneticfield) # Magnetic field strength
    b_vec = magneticfield / B     # Magnetic field direction (unit vector)
    ExBdrift = exbdrift(magneticfield, electricfield)
    vel_in_E_frame = vel - ExBdrift
    # Calculate the velocity parallell to the magnetic field -- vparal
    vparal = vel ⋅ b_vec
    # Calculate mangetic moment -- mu
    vperp = vel_in_E_frame - vparal * b_vec
    mu = magneticmoment(norm(vperp), mass, B)
    return vparal, mu, vel_in_E_frame
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
    itpvec
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
    velocity = vparal*b + v_exb + sign(charge)*vperp*(cθ*e₂ + sθ*e₁)
    # In-place return
    u[:] = [position; velocity]
end

"""
    cosineof_pitchangle(
        magneticfield    ::Vector{<:Real},
        parallel_velocity::Real,
        mass             ::Real,
        magneticmoment   ::Real,
        )
Calculates the cosine of the pitch angle of a charge particle in a magnetic
field, given the particle's parallel velocity, mass and magnetic moment.
"""
function cosineof_pitchangle(
    magneticfield::Vector{<:Real},
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
    return 1 / √(1 - speed^2 * tp.cSqrdInv)
end

"""
    kineticspeed(kineticenergy::Real, mass::Real)
Calculate the speed of a particle from its `kineticenergy` and `mass`.
"""
function kineticspeed(kineticenergy::Real, mass::Real)
    return √(2kineticenergy / mass)
end

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
function larmorradius(mass, vperp, charge, magneticfieldstrength)
    return abs(mass*vperp/(charge*magneticfieldstrength))
end


"""
    gyrofrequency(mass, charge, magneticfieldstrength)
Return the gyrofrequency of a charged particle.
"""
function gyrofrequency(mass, charge, magneticfieldstrength)
    return abs(charge*magneticfieldstrength/mass)
end


"""
    characteristicfieldlength(fieldstrength, fieldstrengthgradient)
Return the characteristic length of a field.
"""
function characteristicfieldlength(fieldstrength, fieldstrengthgradient)
    return fieldstrength / norm(fieldstrengthgradient)
end


"""
    kineticenergy(velocity, mass)
Return the kinetic energy of a particle with `velocity` and `mass`.
"""
function kineticenergy(velocity, mass)
    0.5*mass*velocity^2
end


"""
    kineticenergy(velocity::Vector, mass)
Return the kinetic energy of a particle given the velocity *vector* and mass.
"""
function kineticenergy(velocity::Vector, mass)
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
velocity of its *guiding centre*. Other drifts are neglected
"""
function kineticenergy(
    parallel_velocity::Real,
    vperp            ::Real,
    exbdrift         ::Vector{<:Real},
    ∇Bdrift          ::Vector{<:Real},
    curvaturedrift   ::Vector{<:Real},
    polarisationdrift::Vector{<:Real},
    mass             ::Real,
    )
    vdrift = norm(exbdrift + ∇Bdrift + curvaturedrift + polarisationdrift)
    velsquared = parallel_velocity^2 + vperp^2 + norm(vdrift)^2
    return 0.5*mass*velsquared
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
    vperp            ::Real,
    exbdrift         ::Vector{<:Real},
    mass             ::Real,
    )
    velsquared = parallel_velocity^2 + vperp^2 + norm(exbdrift)^2
    return 0.5*mass*velsquared
end


"""
    kineticenergy(
        bfield           ::Vector{<:Real},
        efield           ::Vector{<:Real},
        parallel_velocity::Real,
        magneticmoment   ::Real,
        mass             ::Real,
        )
Return the kinetic energy of charged particle given the `parallel_velocity`
of its *guiding centre*, its `magnetic_moment`, its `mass`, and the external
electromagnetic field. Drifts other than the E cross B drift are neglected.
"""
function kineticenergy(
    bfield           ::Vector{<:Real},
    efield           ::Vector{<:Real},
    parallel_velocity::Real,
    magneticmoment   ::Real,
    mass             ::Real,
    )
    velsquared = parallel_velocity^2 + perpendicular_velocity(
        magneticmoment,
        mass,
        norm(bfield)
        )^2 + norm(exbdrift(bfield, efield)[1])^2
    return 0.5*mass*velsquared
end


"""
Kinetic energy of a charged particle based on its guiding centre, parallel
velocity, magnetic moment, and the electromagnetic field in which it is
embedded. Accounts for E cross B drift, magnetic gradient drift, curvature
drift, and polarisationdrift.
"""
function kineticenergy(
    R     ::Vector{<:Real},
    vparal::Real,
    q     ::Real, # Particle charge
    m     ::Real, # Particle mass
    μ     ::Real, # Particle magnetic moment
    itpvec::Vector{<:AbstractInterpolation}, # electromagnetic field interpolators
    )
    vperp, drifts = drifts(R, vparal, q, m, μ, itpvec)
    kineticenergy(vparal, vperp, drifts...)
end


"""
    perpendicular_velocity(magnetic_moment, mass, magneticfieldstrength)
Return the perpendicular velocity of a charged particle in a magnetic field.
"""
function perpendicular_velocity(magnetic_moment, mass, magneticfieldstrength)
    sqrt(2magnetic_moment*magneticfieldstrength/mass)
end


"""
    magneticmoment(perpendicular_velocity, mass, magneticfieldstrength)
Return the magnetic moment of a chargred particle in a magnetic field.
"""
function magneticmoment(perpendicular_velocity, mass, magneticfieldstrength)
    0.5mass*perpendicular_velocity^2/magneticfieldstrength
end


"""
exbdrift(
    b̂::Vector{<:Real},
    E_vec::Vector{<:Real},
    B_inv::Real,
    )
"""
function exbdrift(
    b̂::Vector{<:Real},
    E_vec::Vector{<:Real},
    B_inv::Real,
    )
    return (E_vec × b̂)*B_inv
end


"""
    exbdrift(magneticfield::Vector{<:Real}, electricfield::Vector{<:Real})
Calculate the E cross B drift in a given electromagnetic field. Also return
the magnetic field strength and mangnetic field direction.
"""
function exbdrift(magneticfield::Vector{<:Real}, electricfield::Vector{<:Real})
    B = norm(magneticfield) # Magnetic field strength
    B_inv = 1/B
    b̂ = magneticfield*B_inv     # Magnetic field direction (unit vector)
    return exbdrift(b̂, electricfield, B_inv), B, B_inv, b̂
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
    b̂    ::Vector{<:Real}, # The direction of the magnetic fielda
    ∇B   ::Vector{<:Real}, # The gradient of the magnetic field strength
    μ    ::Real, # the magnetic moment of the particle
    B_inv::Real, # The inverse of the magnetic field strength
    q_inv::Real  # The inverse of the charge of the particle
    )
    return q_inv*B_inv*μ*(b̂ × ∇B)
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
    b̂     ::Vector{<:Real},  # The direction of the magnetic fielda
    db̂dt  ::Vector{<:Real},# Material derivative of the magn. field direction
    vparal::Real, # Particle velocity component parallel to the magnetic field
    B_inv ::Real, # The inverse of the magnetic field strength
    q_inv ::Real, # The inverse of the charge of the particle
    mass  ::Real  # The particle mass
    )
    return q_inv*B_inv*mass*b̂ × (vparal*db̂dt)
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
    b̂     ::Vector{<:Real},  # The direction of the magnetic fielda
    dExBdt::Vector{<:Real},# Material derivative of the E cross B drift
    B_inv ::Real, # The inverse of the magnetic field strength
    q_inv ::Real, # The inverse of the charge of the particle
    mass  ::Real  # The particle mass
    )
    return q_inv*B_inv*mass*b̂ × dExBdt
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
    b̂ ::Vector{<:Real},  # The direction of the magnetic fielda
    ∇B::Vector{<:Real},  # The gradient of the magnetic field strength
    μ ::Real, # The magnetic moment of the particle
    m ::Real, # The particle mass
    )
    return -μ*b̂⋅∇B/m
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
    b̂    ::Vector{<:Real},  # The direction of the magnetic fielda
    E_vec::Vector{<:Real},# The electric field
    q    ::Real, # The charge of the particle
    m    ::Real, # The particle mass
    )
    return q*(E_vec⋅b̂)/m
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
    pos          ::Vector{<:Real},
    vel          ::Vector{<:Real},
    magneticfield::Vector{<:Real},
    electricfield::Vector{<:Real},
    charge       ::Real,
    mass         ::Real
    )
    ExBdrift, B, _, b_vec = exbdrift(magneticfield, electricfield)
    vel_in_E_frame = vel - ExBdrift
    # Calculate the guiding centre posistion 
    R = pos - mass/(charge*B) * (vel_in_E_frame × b_vec)
    # Calculate the velocity parallell to the magnetic field -- vparal
    vparal = vel ⋅ b_vec
    # Calculate mangetic moment -- mu
    vperp = vel_in_E_frame - vparal*b_vec
    mu = magneticmoment(norm(vperp), mass, B)
    return R, vparal, mu
end


"""
    get_guidingcentre_2Dxz(
        pos          ::Vector{<:Real},
        vel          ::Vector{<:Real},
        emfields_itpvec::Vector{<:AbstractInterpolation},
        charge       ::Real,
        mass         ::Real
        )
"""
function get_guidingcentre_2Dxz(
    pos          ::Vector{<:Real},
    vel          ::Vector{<:Real},
    emfields_itpvec::Vector{<:AbstractInterpolation},
    charge       ::Real,
    mass         ::Real
    )
    posx, poz = pos[1], pos[3]
    magneticfield = [emfields_itpvec[i](posx, poz) for i = 1:3]
    electricfield = [emfields_itpvec[i](posx, poz) for i = 4:6]
    ExBdrift, B, b_vec = exbdrift(magneticfield, electricfield)
    vel_in_E_frame = vel - ExBdrift
    # Calculate the guiding centre posistion
    R = pos - mass/(charge*B) * (vel_in_E_frame × b_vec)
    # Calculate the velocity parallell to the magnetic field -- vparal
    vparal = vel ⋅ b_vec
    # Calculate mangetic moment -- mu
    vperp = vel_in_E_frame - vparal*b_vec
    mu = magneticmoment(norm(vperp), mass, B)
    return R, vparal, mu
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
    R            ::Vector{<:Real}, # Guiding centre position
    vparal       ::Real,           # Velocity parallel to the magnetic field
    μ            ::Real,           # Magnetic moment of particle
    mass         ::Real,
    charge       ::Real,
    phaseangle   ::Real,           # Arbitrary phase angle of gyration.
    )
    v_exb, B, b = exbdrift(magneticfield, electricfield) # get E cross B drift,
    # magnetic field strength and magnetic field direction
    vperp = perpendicular_velocity(μ, mass, B)
    r_L = larmorradius(mass, vperp, charge, B)
    e₁ = (R × b)/norm(R)
    e₂ = e₁ × b
    cθ = cos(phaseangle)
    sθ = sin(phaseangle)
    position = R + r_L*(cθ*e₁ + sθ*e₂)
    velocity = vparal*b + v_exb + vperp*(cθ*e₂ - sθ*e₁)
    return [position; velocity]
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
    magneticfield    ::Vector{<:Real},
    parallel_velocity::Real,
    mass             ::Real,
    magneticmoment   ::Real,
    )
    B = norm(magneticfield)
    return sqrt(1/(2B*magneticmoment/(mass*parallel_velocity^2) + 1))
end


"""
    maxwellboltzmanndistr(
    v          ::AbstractArray{<:Real},
    temperature::Real,
    mass       ::Real
    )
Return the probability of the velocity `v` of a Maxwell-Boltzmann 
distribution with given particle `temperature` and `mass`.
"""
function maxwellboltzmanndistr(
    v          ::AbstractArray{<:Real},
    temperature::Real,
    mass       ::Real
    )
    σ = √(k_B*temperature/mass) # Standard deviation of velocity
    return @.  (1/(2π))^(3/2) *σ^(-3) * exp(-0.5(v/σ)^2) * 4π*v^2
end


"""
    magneticmirrorfield(
        x::Real,
        y::Real,
        z::Real,
        B0::Real,
        L ::Real
        )
Return the magnetic field at (`x`, `y`, `z`) in a magnetic mirror/bottle with 
length `L` and strength `B0`.

B⃗(x, y, z) = -xzB₀/L² x̂  - yzB₀/L² ŷ  + (B₀ + zB₀/L²)ẑ
"""
function magneticmirrorfield(
    x::Real,
    y::Real,
    z::Real,
    B0::Real,
    L ::Real,
    )
    a = B0*z/L^2
    return [-x*a, -y*a, B0 + z*a]
end # mirroringfield


"""
    magneticdipolefield(
        x::Real,
        y::Real,
        z::Real,
        M ::Real
        )
Return the magnetic field at (`x`, `y`, `z`) in a magnetic dipole field with 
mangetic moment `M`.

B⃗(x, y, z) = 3αzx x̂  + 3αzy ŷ  + α(2z² - x² - y²)ẑ, 
where 
α = M/(x² + y² + z²)^(2.5).
"""
function magneticdipolefield(
    x::Real,
    y::Real,
    z::Real,
    M ::Real,
    )
    a = M/(x^2 + y^2 + z^2)^(5/2)
    return [3a*z*x, 3a*z*y, a*(2z^2 - x^2 - y^2)]
end

"""
    magneticfield_2Dplasmoid(
        x::Real,
        y::Real,
        z::Real,
        μ_x::Real,
        μ_y::Real,
        σ_x::Real,
        σ_y::Real,
        normfactor::Real,
        B0::Vector{<:Real},
        )
Return the magnetic field at (`x`, `y`, `z`) in a plasmoid configuration with
background field `B0`. The field is analytical with a vector potential

A⃗(x,y,z) = normfactor * N_x(μ_x, σ_x)*N_y(μ_y, σ_y) * ẑ - B0_x * zŷ + B0_y * zx̂

where N_x is normal distribution of x around mean `μ_x` with standard deviation
`σ_x`.
"""
function magneticfield_2Dplasmoid(
    x::Real,
    y::Real,
    z::Real,
    μ_x::Real,
    μ_y::Real,
    σ_x::Real,
    σ_y::Real,
    normfactor::Real,
    B0::Vector{<:Real},
    )
    A_z = normfactor/(σ_x*σ_y*2π) * exp(-0.5*(((x - μ_x)/σ_x)^2 + ((y - μ_y)/σ_y)^2))
    return [-A_z*(y - μ_y)/σ_y^2 + B0[1], A_z*(x - μ_x)/σ_x^2 + B0[2], B0[3]]
end


"""
    magneticfieldstrength_2Dplasmoid(args...)
Return the magnetic field strength of the 2D plasmoid configuration defined in
@magneticfield_2Dplasmoid.
"""
function magneticfieldstrength_2Dplasmoid(args...)
    B = magneticfield_2Dplasmoid(args...)
    return norm(B)
end

"""
    magneticfieldstrengthgradient_2Dplasmoid(
        x::Real,
        y::Real,
        z::Real,
        μ_x::Real,
        μ_y::Real,
        σ_x::Real,
        σ_y::Real,
        normfactor::Real,
        B0::Vector{<:Real},
        )
Return the gradient magnetic field strength of the 2D plasmoid configuration
defined in @magneticfield_2Dplasmoid.
"""
function magneticfieldstrengthgradient_2Dplasmoid(
    x::Real,
    y::Real,
    z::Real,
    μ_x::Real,
    μ_y::Real,
    σ_x::Real,
    σ_y::Real,
    normfactor::Real,
    B0::Vector{<:Real},
    )
    frac_x = (x - μ_x)/σ_x
    frac_y = (y - μ_y)/σ_y
    A_z = normfactor/(2π*σ_x*σ_y) * exp(-0.5*(frac_x^2 + frac_y^2))
    frac_3 = 2A_z*frac_x*frac_y/(σ_x*σ_y)
    α = -A_z * frac_y/σ_y + B0[1]
    β = A_z * frac_x/σ_x + B0[2]
    γ = sqrt(α^2 + β^2 + B0[3]^2)
    dβ2dx = 2β*A_z/σ_x^2 * (1 - frac_x^2)
    dβ2dy = -β*frac_3
    dα2dx = α*frac_3
    dα2dy = -2α*A_z/σ_y^2 * (1 - frac_y^2)

    ∇B_x = dα2dx + dβ2dx
    ∇B_y = dα2dy + dβ2dy
    #∇B_x = (x - μ_x)/σ_x^2 * (1/(γ*σ_x^2) - γ)
    #∇B_y = (y - μ_y)/σ_y^2 * (1/(γ*σ_y^2) - γ)
    return 1/(2γ)*[∇B_x, ∇B_y, 0]
end


"""
    fadeevEquilibrium(
        (x0, y0, z0)::Tuple{Real, Real, Real},
        (xf, yf, zf)::Tuple{Real, Real, Real},
        (nx, ny, nz)::Tuple{Integer, Integer, Integer},
        λ           ::Real,
        ϵ           ::Real,
        B0          ::Real
        )
"""
function fadeevEquilibrium(
    (x0, y0, z0)::Tuple{Real, Real, Real},
    (xf, yf, zf)::Tuple{Real, Real, Real},
    (nx, ny, nz)::Tuple{Integer, Integer, Integer},
    λ           ::Real,
    ϵ           ::Real,
    B0          ::Real
    )
    # Create spatial axes and find the grid sizes
    xx, yy, zz, dx, dy, dz = createaxes((x0, y0, z0),
                                        (xf, yf, zf),
                                        (nx, ny, nz)
                                        )
    # Initialise the vector field
    ndims = 3
    A = zeros(ndims, nx, ny, nz)
    for i = 1:nx
        for j = 1:ny
            A[3,i,j,:] .= B0 * λ * log2( ϵ * cos(xx[i]/λ) + cosh(yy[j]/λ) )
        end
    end
    return (xx, yy, zz), (dx, dy, dz), A
    
end # function FadeevEquilibrium


"""
    gcagradients_from_emfield(
        bField ::Array{T, 4} where {T<:Real},
        eField ::Array{T, 4} where {T<:Real},
        xCoords::Vector{T} where {T<:Real},
        yCoords::Vector{T} where {T<:Real},
        zCoords::Vector{T} where {T<:Real},
        scheme ::Function
        )
Compute the gradients used in the guiding centre approximation from an
electromagnetic field. The gradients are:
- The gradient of the magnetic field strength
- the gradient of the magnetic field direction
- the gradient of the E cross B drift.
"""
function gcagradients_from_emfield(
    bField ::Array{T, 4} where {T<:Real},
    eField ::Array{T, 4} where {T<:Real},
    xCoords::Vector{T} where {T<:Real},
    yCoords::Vector{T} where {T<:Real},
    zCoords::Vector{T} where {T<:Real},
    scheme ::Function
    )
    wfp = typeof(bField[1])
    nx, ny, nz, ncomp = size(bField)
    println("---> Creating mesh of B")
    BB = norm4(bField, axis=4)
    b̂ = zeros(wfp, nx, ny, nz, ncomp)
    ExBdrift = zeros(wfp, nx, ny, nz, ncomp)
    # I guess this could be done much more efficiently
    println("---> Creating mesh of b, and ExB")
    for i = 1:nx
        for j= 1:ny
            for k = 1:nz
                B⃗ = bField[i,j,k,:]
                E⃗ = eField[i,j,k,:]
                B = BB[i,j,k]
                b̂[i,j,k,:]  .= B⃗ ./ B
                ExBdrift[i,j,k,:] .= (E⃗ × B⃗) ./ B^2
            end
        end
    end
    println("---> Differentiating B")
    ∇B = ∇(BB, xCoords, yCoords, zCoords, scheme)
    println("---> Differentiating b")
    ∇b̂ = ∇(b̂,  xCoords, yCoords, zCoords, scheme)
    println("---> Differentiating ExB")
    ∇ExBdrift= ∇(ExBdrift, xCoords, yCoords, zCoords, scheme)
    return ∇B, ∇b̂, ∇ExBdrift
end

"""
    lorentzfactor(speed::Real)
The relativistic Lorentz factor in SI units.
"""
function lorentzfactor(speed::Real)
    return 1/√(1-speed^2*tp.cSqrdInv)
end

"""
    kineticspeed(kineticenergy::Real, mass::Real)
Calculate the speed of a particle from its `kineticenergy` and `mass`.
"""
function kineticspeed(kineticenergy::Real, mass::Real)
    return √(2kineticenergy/mass)
end

#-------------------------------------------------------------------------------
# Created 02.12.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 equations_of_motion.jl
#
#-------------------------------------------------------------------------------
# Contains various equations of motion.
#-------------------------------------------------------------------------------

"""
    fieldlinetracing_forward(u, p, _)
The equation of motion for tracing a field line forward. The statevector `u` 
contains the position in 3D, and `p` contains an interpolation functor giving 
the field at the position.
using Base: Forward
"""
function fieldlinetracing_forward(u, p, _)
    field_itp = p
    x, y, z = u
    field = field_itp(x, y, z)
    fieldstrength = norm(field)
    fielddirection = field ./ fieldstrength
    dsdt = fielddirection
    return dsdt
end


"""
    fieldlinetracing_backward(u, p, _)
The equation of motion for tracing a field line backward. The statevector `u` 
contains the position in 3D, and `p` contains an interpolation functor giving 
the field at the position.
"""
function fieldlinetracing_backward(u, p, _)
    field_itp = p
    x, y, z = u
    field = field_itp(x, y, z)
    fieldstrength = norm(field)
    fielddirection = field ./ fieldstrength
    dsdt = -fielddirection
    return dsdt
end


"""
    emfieldatpos(
        pos::Vector{<:Real},
        fields_itp,
    )::Tuple(Vector{<:Real, <:Real}
Interpolate the electromagetic field to `pos`. Returns the magnetic and
electric field vectors. Methods for both a vector of interpolationobjects,
and just a single interpolationobject.
"""
function emfieldatpos(
    pos::Vector{<:Real},
    fields_itp::Vector{<:AbstractInterpolation}
)
    fields = [itp(pos...) for itp in fields_itp]
    B = fields[1:3]
    E = fields[4:6]
    return B, E
end
function emfieldatpos(
    pos::Vector{<:Real},
    fields_itp::AbstractInterpolation
)
    fields = fields_itp(pos...)
    B = fields[1:3]
    E = fields[4:6]
    return B, E
end
"""
    lorentzforce!(du, u, p, _)
Equations of motion described by the Lorentz Force in 3 dimensions. The 
parameters `p` contains the charge and mass of the particle, and an 
interpolation functor giving the magnetic and electric field by passing the 
position coordinates.
"""
function lorentzforce!(du, u, p, _)
    q, m = p.charge, p.mass
    x = u[1:3] # The position vector
    v = u[4:6] # The velocity vector
    B, E = emfieldatpos(x, p.fields)
    dvdt = q / m * (E + v × B)
    dxdt = v
    du[:] = [dxdt; dvdt]
end
function lorentzforce(u, p, _)
    q, m, field_itp = p
    x = u[1:3] # The position vector
    v = u[4:6] # The velocity vector
    fields = field_itp(x...)
    B = fields[1:3]
    E = fields[4:6]
    dvdt = q / m * (E + v × B)
    dxdt = v
    return [dxdt; dvdt]
end


"""
    guidingcentreapproximation!(du, u, p, _)
The equations of motion of a charged particle in collisionless plasma under 
the guiding centre approximation (GCA). Compared to the full motion (described 
by the Lorentz force), the statevector `u` is reduced from 6 DoF to 4, namely 
the guiding centre position and the guiding centre velocity parallel to the 
magnetic field, `[Rx, Ry, Rz, vparal]`, respectively.

The parameters `p` are the particle charge, mass and magnetic moment (which is 
assumed to be constant in the GCA), along with the interpolation functor 
from Interpolations.jl giving the magnetic and electric field at an arbitrary 
position. The functor should give a 6-component vector, where the first 3 
components represents the magnetic field and the last 3 the electric field.

The function uses ForwardDiff and Interpolations to calculate gradients.
"""
function guidingcentreapproximation!(du, u, p, _)
    R = u[1:3]
    vparal = u[4] # Particle velocity parallel to the magnetic field
    # Extract parameters
    q, m, μ, itpvec = p.charge, p.mass, p.magneticmoment, p.fields
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

    # We index du with 1:4 to allow aribitrary size of u.
    # (This is necessary when doing hybrid swtich method).
    du[1:4] = gca_drift_and_acceleration(
        ∇b, ∇ExB, ∇B, B_vec, E_vec, vparal, q, m, μ
    )
end


"""
    hybridgcafo!(du, u, p, _)
A hybrid EoM. Runs either the guiding centre approximation or full orbit
integration depending on the parameter `switch`.
"""
function hybridgcafo!(du, u, p, _)
    if p.switch == 1
        guidingcentreapproximation!(du, u, p, nothing)
    elseif p.switch == 2
        lorentzforce!(du, u, p, nothing)
    end
end


"""
    gca_2Dxz!(du, u, p, _)
The equations of motion of a charged particle in collisionless plasma under 
the guiding centre approximation (GCA) in 2.5D xz-plane. Compared to the full 
motion (described by the Lorentz force), the statevector `u` is reduced from 6 
DoF to 4, namely the guiding centre position and the guiding centre velocity 
parallel to the magnetic field, `[Rx, Ry, Rz, vparal]`, respectively. In 
practice, `Ry` is necessary.

The parameters `p` are the particle charge, mass and magnetic moment (which is 
assumed to be constant in the GCA), along with the interpolation functor 
from Interpolations.jl giving the magnetic and electric field at an arbitrary 
position. The functor should give a 6-component vector, where the first 3 
components represents the magnetic field and the last 3 the electric field.

The function uses ForwardDiff and Interpolations to calculate gradients.
"""

function gca_2Dxz!(du, u, p, _)
    # NOTE: Comments are copied from 3D implementation. Running times are
    # not correct in this case since this is 2D version.
    #
    Rx, Rz = u[1], u[3]
    vparal = u[4] # Particle velocity parallel to the magnetic field
    # Extract parameters
    q, m, μ, itpvec = p.charge, p.mass, p.magneticmoment, p.fields
    # Use the gyrocentre position interpolate the vectors
    # scalars from the interpolation objects.
    B_vec = [itpvec[i](Rx, Rz) for i in 1:3]
    E_vec = [itpvec[i](Rx, Rz) for i in 4:6]

    # Calculate the field gradients.
    jacobian_matrix = ForwardDiff.jacobian([Rx, Rz]) do x
        vec = [itp(x...) for itp in itpvec]
        B_vec_fd = vec[1:3]
        E_vec_fd = vec[4:6]
        B_fd = norm(B_vec_fd)
        b_fd = B_vec_fd / B_fd
        return [b_fd; E_vec_fd × b_fd / B_fd; B_fd]
    end
    jacobian_matrix = [
        jacobian_matrix[:, 1];;
        zeros(typeof(Rx), 7);;
        jacobian_matrix[:, 2]
    ]
    ∇b = jacobian_matrix[1:3, :]
    ∇ExB = jacobian_matrix[4:6, :]
    ∇B = jacobian_matrix[7, :]

    du[:] = gca_drift_and_acceleration(
        ∇b, ∇ExB, ∇B, B_vec, E_vec, vparal, q, m, μ
    )
end


"""
    gca_highmemory(u, p, _)
Same as `guidingcentreapproximation!`, but does not calculate gradients
on the fly. Uses interpolation and gradients stored in memory instead.
"""
function gca_highmemory!(du, u, p, _)
    R = u[1:3]    # Position of the gyrocentre
    vparal = u[4] # Particle velocity parallel to the magnetic field
    beta = u[5]   # Cosine of pitch angle

    # Extract parameters:
    #
    # Particle charge and particle mass and the Coulomb logarithm
    q, m, mu = p[1:3]
    # interpolation object for:
    # the magnetic vector field, the electric vector field, the magnetic
    # gradient vector field, the gradient of the magnetic direction
    # (a 3rd order tensor field) and the gradient of the ExB-drift (a 3rd
    # order tensor field) 
    field_itp = p[4]
    fields = field_itp(R...)

    # Use the gyrocentre position interpolate the vectors, matrices and
    # scalars from the interpolation objects.
    B_vec = fields[1:3]
    E_vec = fields[4:6]
    gradB_vec = fields[7:9]
    gradb = fields[10:18]
    gradExB = fields[19:27]

    # Calculate some quantities
    B = norm(B_vec)   # The magnetic field strength
    b_vec = B_vec / B   # An unit vector pointing in the direction of the
    #  magnetic field
    Eparal = E_vec ⋅ b_vec # Electric field component parallel to the
    # magnetic field

    # Calculate drifts
    ExBdrift = (E_vec × b_vec) / B
    ∇Bdrift = mu / (q * B) * (b_vec × gradB_vec)

    # Material derivatives of the magnetic field strength, the magnetic field
    # direction, and the ExB-drift. Assumes ∂/∂t = 0
    # See Ripperda et al. (2018) and notes.
    dbdt = vparal * (gradb * b_vec) + gradb * ExBdrift
    dExBdt = vparal * (gradExB * b_vec) + gradExB * ExBdrift

    # Compute the perpendicular velcoity
    dRperpdt = ExBdrift + ∇Bdrift + m * b_vec / (q * B) × (vparal * dbdt + dExBdt)

    # Compute the acceleration along the magnetic field
    dvparaldt = (q * Eparal - mu * b_vec ⋅ gradB_vec) / m
    # With correction proposed by Birn et al., 2004:
    #dvparaldt = (q*Eparal - μ*b̂⋅∇B)/m + ExBdrift⋅db̂dt + ∇Bdrift⋅db̂dt
    #dRperpdt = b̂/B × (-c*E⃗ + μ*c/q * ∇B) #old

    # Compute the total velocity of the guiding centre
    dRdt = vparal * b_vec + dRperpdt

    # Update the statevector
    du[:] = [dRdt; dvparaldt]

end


"""
    gca_highmemory_2Dxz!(du, u, p, _)
Same as @gca_2Dxz!, but does not calculate gradients on the fly. Uses
interpolation and gradients stored in memory instead.
"""
function gca_highmemory_2Dxz!(du, u, p, _)
    Rx, _, Rz = u[1:3]    # Position of the gyrocentre
    vparal = u[4] # Particle velocity parallel to the magnetic field

    # Extract parameters:
    #
    # Particle charge and particle mass and the Coulomb logarithm
    q = p.charge
    m = p.mass
    mu = p.magneticmoment
    itpvec = p.fields
    # interpolation object for:
    # the magnetic vector field, the electric vector field, the magnetic
    # gradient vector field, the gradient of the magnetic direction
    # (a 3rd order tensor field) and the gradient of the ExB-drift (a 3rd
    # order tensor field)
    B_vec = [itpvec[i](Rx, Rz) for i in 1:3]
    E_vec = [itpvec[i](Rx, Rz) for i in 4:6]
    gradB = [itpvec[i](Rx, Rz) for i in 7:9]
    gradb = reshape([itpvec[i](Rx, Rz) for i in 10:18], 3, 3)
    gradExB = reshape([itpvec[i](Rx, Rz) for i in 19:27], 3, 3)

    # Calculate some quantities
    B = norm(B_vec)   # The magnetic field strength
    b = B_vec / B   # An unit vector pointing in the direction of the
    #  magnetic field
    Eparal = E_vec ⋅ b  # Electric field component parallel to the
    # magnetic field

    # Calculate drifts
    ExBdrift = (E_vec × b) / B
    ∇Bdrift = mu / (q * B) * (b × gradB)

    # Material derivatives of the magnetic field strength, the magnetic field
    # direction, and the ExB-drift. Assumes ∂/∂t = 0
    # See Ripperda et al. (2018) and notes.
    #dBdt = vparal * b_vec ⋅ gradB_vec + ExBdrift ⋅ gradB_vec
    dbdt = vparal * (gradb * b) + gradb * ExBdrift
    dExBdt = vparal * (gradExB * b) + gradExB * ExBdrift

    # Compute the perpendicular velcoity
    dRperpdt = ExBdrift + ∇Bdrift + m * b / (q * B) × (vparal * dbdt + dExBdt)

    # Compute the acceleration
    #dvparaldt = (q*Eparal - μ*b⋅∇B)/m + ExBdrift⋅dbdt
    # With correction proposed by Birn et al., 2004:
    dvparaldt = (q * Eparal - μ * b ⋅ ∇B) / m + (ExBdrift + ∇Bdrift) ⋅ dbdt
    # Old version:
    #dvparaldt = (q*Eparal - μ*b⋅∇B)/m # along the magnetic field lines

    # Compute the total velocity of the guiding centre
    dRdt = vparal * b + dRperpdt

    # Update the statevector
    return du[:] = [dRdt; dvparaldt]

end

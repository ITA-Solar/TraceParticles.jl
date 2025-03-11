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
    fields_itp::Vector{<:Any}
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
    # Calculate the field gradients.
    ∇b, ∇ExB, ∇B, B_vec, E_vec = fieldgradients(R, itpvec)

    # We index `du` with 1:4 to allow aribitrary size of u.
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
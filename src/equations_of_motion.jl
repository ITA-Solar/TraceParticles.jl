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
function fieldlinetracing_forward!(du, u, p, _)
    time = p.time
    x, y, z = u[1], u[2], u[3]
    field = p.field(x, y, z, time)
    fieldstrength = norm(field)
    fielddirection = field ./ fieldstrength
    dsdt = fielddirection
    du[1:3] .= dsdt
end
function fieldlinetracing_backward!(du, u, p, _)
    time = p.time
    x, y, z = u[1], u[2], u[3]
    field = p.field(x, y, z, time)
    fieldstrength = norm(field)
    fielddirection = field ./ fieldstrength
    dsdt = -fielddirection
    du[1:3] .= dsdt
end

"""
    lorentzforce!(du, u, p, t)
Equations of motion described by the Lorentz Force in 3 dimensions. The 
parameters `p` contains the charge and mass of the particle, and an 
interpolation functor giving the magnetic and electric field by passing the 
position coordinates.
"""
function lorentzforce!(du, u, p, t)
    q, m = p.charge, p.mass
    v = SVector(u[4], u[5], u[6]) # The velocity vector
    E, B = p.electromagneticfield(u[1], u[2], u[3], t)
    dvdt = q / m * (E + v × B)
    dxdt = v
    du[1:6] .= [dxdt; dvdt]
end


"""
    guidingcentreapproximation!(du, u, p, t)
The equations of motion of a charged particle in collisionless plasma under 
the guiding centre approximation (GCA). Compared to the full motion (described 
by the Lorentz force), the statevector `u` is reduced from 6 DoF to 4, namely 
the guiding centre position and the guiding centre velocity parallel to the 
magnetic field, `[Rx, Ry, Rz, vparal]`, respectively.

The parameters `p` are the particle charge, mass and magnetic moment (which is
assumed to be constant in the GCA), along with an interpolation functor/functor
taking the arguments `(x, y, z, t)` and giving the magnetic and electric field
at an arbitrary position. The function/functor should return two
`AbstractVector`s, representing the 3-component electric and magnetic field,
respectievely. If the function/functor returns `StaticArrays`, `fieldgradients`
will also return `StaticArrays`.

The function uses ForwardDiff and Interpolations to calculate gradients.
"""
function guidingcentreapproximation!(du, u, p, t)
    x, y, z = u[1], u[2], u[3] # The guiding centre position vector
    vparal = u[4] # Particle velocity parallel to the magnetic field
    q, m, μ = p.charge, p.mass, p.magneticmoment
    # Calculate the field gradients.
    ∇b, ∇ExB, ∇B, B_vec, E_vec, ∂b, ∂ExB, _ = fieldgradients(
        x, y, z, t, p.electromagneticfield
    )
    # We index `du` with 1:4 to allow arbitrary size of u.
    # (This is necessary when doing hybrid swtich method).
    du[1:4] .= gca_drift_and_acceleration(
        ∇b, ∇ExB, ∇B, B_vec, E_vec, ∂b, ∂ExB, vparal, q, m, μ
    )
    return nothing
end


"""
    hybridgcafo!(du, u, p, t)
A hybrid EoM. Runs either the guiding centre approximation or full orbit
integration depending on the parameter `switch`.
"""
function hybridgcafo!(du, u, p, t)
    if p.eomid == 1
        guidingcentreapproximation!(du, u, p, t)
        du[5:6] .= 0.0
    elseif p.eomid == 2
        lorentzforce!(du, u, p, t)
    end
end

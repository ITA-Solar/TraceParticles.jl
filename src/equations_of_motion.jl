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
    field = field_itp(x,y,z)
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
    field = field_itp(x,y,z)
    fieldstrength = norm(field)
    fielddirection = field ./ fieldstrength
    dsdt = -fielddirection
    return dsdt
end
    

"""
    lorentzforce!(du, u, p, _)
Equations of motion described by the Lorentz Force in 3 dimensions. The 
parameters `p` contains the charge and mass of the particle, and an 
interpolation functor giving the magnetic and electric field by passing the 
position coordinates.
"""
function lorentzforce!(du, u, p, _)
    q, m, field_itp = p
    x = u[1:3] # The position vector
    v = u[4:6] # The velocity vector
    fields = field_itp(x...)
    B = fields[1:3]
    E = fields[4:6]
    dvdt = q/m * (E + v × B) 
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
    dvdt = q/m * (E + v × B) 
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
    q, m, μ, fields_itp = p
    q_inv = 1/q # Inverse of q - to replace division with multiplication
    # Use the gyrocentre position interpolate the vectors
    # scalars from the interpolation objects.
    fields = fields_itp(R...)
    B_vec = fields[1:3]
    E_vec = fields[4:6]
    B = norm(B_vec)   # The magnetic field strength
    B_inv = 1/B       # Inverse of B - to replace divition with multiplication
    b = B_vec*B_inv   # An unit vector pointing in the direction of the
                      #  magnetic field
    ExBdrift = (E_vec × b)/B # The E cross B-drift

    # Calculate the gradient of the magnetic field strength
    ∇B = ForwardDiff.gradient(R) do x
        norm(fields_itp(x...)[1:3])
    end
    #
    # Calculate the gradient of the magnetic field direction
    #____
    # The following two lines is according to a @Benchmark test with 201^3 
    # grid points 1.17 faster than the ForwardDiff.jacobian method (commented 
    # below). The results of both methods where exactly the same.  
    #
    # Edit: B_itp replaced by field_itp, which means that calculating the 
    # gradient also gives the jacobian matric of the electric field at the same
    # time. Will save time in total but this step might take som extra time.
    #
    #∇b = ForwardDiff.jacobian(R) do x
    #    B_vec = B_itp(x...)
    #    return B_vec/norm(B_vec)
    #end
    #____
    jacobian_matrix = stack(Interpolations.gradient(fields_itp, R...))
    ∇B_vec = jacobian_matrix[1:3,:]
    ∇b = (∇B_vec - b * ∇B')*B_inv
    #
    # Calculate the Jacobian matrix of the ExB-drift
    #____ 
    # The following five lines is according to a @Benchmark test with 201^3
    # grid points 1.53 faster than the ForwardDiff.jacobian method (commented 
    #  below). The results of both methods where exactly the same.  
    #ForwardDiff.jacobian(R) do x
    #    B_vec = B_vec(x...)
    #    E_vec = E_vec(x...)
    #    return (E_vec × B_vec)/(norm(B_vec)^2)
    #end
    #____
    ∇E_vec = jacobian_matrix[4:6,:]
    skewE = skewsymmetric_matrix(E_vec)
    skewb = skewsymmetric_matrix(b)
    ∇ExB = (-skewb*∇E_vec + skewE*∇b - ExBdrift * ∇B')*B_inv

    # Electric field component parallel to the magnetic field
    Eparal = E_vec⋅b
    # Calculate ∇B-drift
    ∇Bdrift = q_inv*B_inv*μ*(b × ∇B)
    # Total time derivatives. Assumes ∂/∂t = 0,
    dbdt = vparal * (∇b * b) + ∇b*ExBdrift
    dExBdt = vparal * (∇ExB * b) + ∇ExB*ExBdrift
    
    # Compute the perpendicular velcoity
    dRperpdt = ExBdrift + ∇Bdrift + q_inv*B_inv*m*b × (vparal*dbdt + dExBdt)
    #dRperpdt = b̂/B × (-E⃗ + μ/q*∇B + m/q*(vparal*db̂dt + dExBdt))  # old

    # Compute the acceleration 
    dvparaldt = (q*Eparal - μ*b⋅∇B)/m # along the magnetic field lines
    # With correction proposed by Birn et al., 2004:
    #dvparaldt = (q*Eparal - μ*b̂⋅∇B)/m + ExBdrift⋅db̂dt + ∇Bdrift⋅db̂dt
    #dRperpdt = b̂/B × (-c*E⃗ + μ*c/q * ∇B) #old

    # Compute the velocity
    dRdt = vparal*b + dRperpdt
    
    # How to store the perpendicular velocity? Would need to know b̂ at
    #   each R calculate vperp at each R. Could be interesting to store this
    #   as an auxiliary variable somehow, to se how the drift evolves.
    du[:] = [dRdt; dvparaldt]
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
    q, m, μ, itpvec = p
    q_inv = 1/q # Inverse of q - to replace division with multiplication
    # Use the gyrocentre position interpolate the vectors
    # scalars from the interpolation objects.
    B_vec = [itpvec[i](Rx, Rz) for i in 1:3]
    E_vec = [itpvec[i](Rx, Rz) for i in 4:6]
    B = norm(B_vec)   # The magnetic field strength
    B_inv = 1/B       # Inverse of B - to replace divition with multiplication
    b = B_vec*B_inv   # An unit vector pointing in the direction of the
                      #  magnetic field
    ExBdrift = (E_vec × b)/B # The E cross B-drift

    # Calculate the field gradients.
    jacobian_matrix = ForwardDiff.jacobian([Rx, Rz]) do x
        vec = [itp(x...) for itp in itpvec]
        B_vec_fd = vec[1:3]
        E_vec_fd = vec[4:6]
        B_fd = norm(B_vec_fd)
        b̂_fd = B_vec_fd/B_fd
        return [b̂_fd; E_vec_fd × b̂_fd/B_fd; B_fd]
    end
    jacobian_matrix = [
        jacobian_matrix[:,1];;
        zeros(typeof(Rx), 7);; 
        jacobian_matrix[:,2]
        ]
    ∇b = jacobian_matrix[1:3,:]
    ∇ExB = jacobian_matrix[4:6,:]
    ∇B = jacobian_matrix[7,:]

    # Electric field component parallel to the magnetic field
    Eparal = E_vec⋅b
    # Calculate drifts
    ∇Bdrift = q_inv*B_inv*μ*(b × ∇B)
    # Total time derivatives. Assumes ∂/∂t = 0,
    dbdt = vparal * (∇b * b) + ∇b*ExBdrift
    dExBdt = vparal * (∇ExB * b) + ∇ExB*ExBdrift
    
    # Compute the perpendicular velcoity
    dRperpdt = ExBdrift + ∇Bdrift + q_inv*B_inv*m*b × (vparal*dbdt + dExBdt)
    #dRperpdt = b̂/B × (-E⃗ + μ/q*∇B + m/q*(vparal*db̂dt + dExBdt))  # old

    # Compute the acceleration 
    #dvparaldt = (q*Eparal - μ*b⋅∇B)/m + ExBdrift⋅dbdt
    # With correction proposed by Birn et al., 2004:
    dvparaldt = (q*Eparal - μ*b⋅∇B)/m + (ExBdrift + ∇Bdrift) ⋅ dbdt
    # Old version:
    #dvparaldt = (q*Eparal - μ*b⋅∇B)/m # along the magnetic field lines

    # Compute the velocity
    dRdt = vparal*b + dRperpdt
    
    # How to store the perpendicular velocity? Would need to know b̂ at
    #   each R calculate vperp at each R. Could be interesting to store this
    #   as an auxiliary variable somehow, to se how the drift evolves.
    du[:] = [dRdt; dvparaldt]
end


"""
    gca_highmemory(u, p, _)
Same as @gca, but does not calculate gradients on the fly.
Uses interpolation and gradients stored in memory instead.
"""
function gca_highmemory(u, p, _)
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
    b_vec = B_vec/B   # An unit vector pointing in the direction of the
                      #  magnetic field
    Eparal = E_vec⋅b_vec # Electric field component parallel to the
                         # magnetic field

    # Calculate drifts
    ExBdrift = (E_vec × b_vec)/B
    ∇Bdrift = mu/(q*B)*(b_vec × gradB_vec)

    # Material derivatives of the magnetic field strength, the magnetic field
    # direction, and the ExB-drift. Assumes ∂/∂t = 0
    # See Ripperda et al. (2018) and notes.
    dBdt = vparal * b_vec ⋅ gradB_vec + ExBdrift ⋅ gradB_vec
    dbdt = vparal * (gradb * b_vec) + gradb*ExBdrift
    dExBdt = vparal * (gradExB * b_vec) + gradExB*ExBdrift

    # Compute the perpendicular velcoity
    dRperpdt = ExBdrift + ∇Bdrift + m*b_vec/(q*B) × (vparal*dbdt + dExBdt)

    # Compute the acceleration along the magnetic field
    dvparaldt = (q*Eparal - mu*b_vec⋅gradB_vec)/m
    # With correction proposed by Birn et al., 2004:
    #dvparaldt = (q*Eparal - μ*b̂⋅∇B)/m + ExBdrift⋅db̂dt + ∇Bdrift⋅db̂dt
    #dRperpdt = b̂/B × (-c*E⃗ + μ*c/q * ∇B) #old

    # Compute the total velocity of the guiding centre
    dRdt = vparal*b_vec + dRperpdt

    # Compute the pitch angle rate of change
    dbetadt = (1/vparal*dvparaldt - 1/2B*dBdt)*beta*(1-beta^2)

    # Update the statevector
    return [dRdt; dvparaldt; dbetadt] 

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
    b = B_vec/B   # An unit vector pointing in the direction of the
                      #  magnetic field
    Eparal = E_vec⋅b  # Electric field component parallel to the
                         # magnetic field

    # Calculate drifts
    ExBdrift = (E_vec × b)/B
    ∇Bdrift = mu/(q*B)*(b × gradB)

    # Material derivatives of the magnetic field strength, the magnetic field
    # direction, and the ExB-drift. Assumes ∂/∂t = 0
    # See Ripperda et al. (2018) and notes.
    #dBdt = vparal * b_vec ⋅ gradB_vec + ExBdrift ⋅ gradB_vec
    dbdt = vparal * (gradb * b) + gradb*ExBdrift
    dExBdt = vparal * (gradExB * b) + gradExB*ExBdrift

    # Compute the perpendicular velcoity
    dRperpdt = ExBdrift + ∇Bdrift + m*b/(q*B) × (vparal*dbdt + dExBdt)

    # Compute the acceleration
    #dvparaldt = (q*Eparal - μ*b⋅∇B)/m + ExBdrift⋅dbdt
    # With correction proposed by Birn et al., 2004:
    dvparaldt = (q*Eparal - μ*b⋅∇B)/m + (ExBdrift + ∇Bdrift) ⋅ dbdt
    # Old version:
    #dvparaldt = (q*Eparal - μ*b⋅∇B)/m # along the magnetic field lines

    # Compute the total velocity of the guiding centre
    dRdt = vparal*b + dRperpdt

    # Update the statevector
    return du[:] = [dRdt; dvparaldt]

end


"""
The following EoMs where constructed for use in an introductory course in 
numerical solutions to stochastic differential equations.
"""

function constantchange(du, _, p, _)
    du .= p
end
function constantchange(_, p, _)
    return p
end

function proportionalchange(du, u, p, _)
    du .= p*u
end
function proportionalchange(u, p, _)
    return p*u
end

function duffingvanderpoloscillator_drift(du, u, p, _)
    X, V = u
    dXdt = V
    dVdt = (X*(p - X^2) - V)
    du .= [dXdt; dVdt]
end
function duffingvanderpoloscillator_drift(u, p, _)
    X, V = u
    dXdt = V
    dVdt = (X*(p - X^2) - V)
    return [dXdt; dVdt]
end

function duffingvanderpoloscillator_diffusion(du, u, p, _)
    X, _ = u
    dXdW = 0.0
    dVdW = p*X
    du .= [dXdW; dVdW]
end
function duffingvanderpoloscillator_diffusion(u, p, _)
    X, _ = u
    dXdW = 0.0
    dVdW = p*X
    return [dXdW; dVdW]
end

function ornsteinuhlenbeck_drift(du, u, p, _)
    _, V = u
    dXdt = V
    dVdt = -V/p
    du .= [dXdt, dVdt]
end
function ornsteinuhlenbeck_drift(u, p, _)
    _, V = u
    dXdt = V
    dVdt = -V/p
    return [dXdt, dVdt]
end

function ornsteinuhlenbeck_diffusion(du, _, p, _)
    dXdW = 0.0
    dVdW = p
    du .= [dXdW, dVdW]
end
function ornsteinuhlenbeck_diffusion(_, p, _)
    dXdW = 0.0
    dVdW = p
    return [dXdW, dVdW]
end


function GCAPitchAngleFriction_lowmemory_2Dxz(u, p, _)
    # NOTE: Comments are copied from 3D implementation. Running times are
    # not correct in this case since this is 2D version.
    #
    Rx, Rz = u[1], u[3]
    vparal = u[4] # Particle velocity parallel to the magnetic field
    beta = u[5]   # Cosine of pitch angle
    # Extract parameters
    q, m, kf, fields_itp = p
    q_inv = 1/q # Inverse of q - to replace division with multiplication
    # Use the gyrocentre position interpolate the vectors
    # scalars from the interpolation objects.
    fields = fields_itp(Rx, Rz)
    B_vec = fields[1:3]
    E_vec = fields[4:6]
    n = fields[7] # Numver density
    T = fields[8] # gas temperature
    B = norm(B_vec)   # The magnetic field strength
    B_inv = 1/B       # Inverse of B - to replace divition with multiplication
    b = B_vec*B_inv   # An unit vector pointing in the direction of the
                      #  magnetic field
    ExBdrift = (E_vec × b)/B # The E cross B-drift
    μ = m*vparal^2/2B*(1/beta^2 - 1) # The magnetic moment

    # Calculate the gradient of the magnetic field strength
    ∇B = ForwardDiff.gradient([Rx, Rz]) do x
        norm(fields_itp(x...)[1:3])
    end
    ∇B = [∇B[1], 0f0, ∇B[2]]
    #
    # Calculate the gradient of the magnetic field direction
    #____
    # The following two lines is according to a @Benchmark test with 201^3 
    # grid points 1.17 faster than the ForwardDiff.jacobian method (commented 
    # below). The results of both methods where exactly the same.  
    #
    # Edit: B_itp replaced by field_itp, which means that calculating the 
    # gradient also gives the jacobian matric of the electric field at the same
    # time. Will save time in total but this step might take som extra time.
    #
    #∇b = ForwardDiff.jacobian(R) do x
    #    B_vec = B_itp(x...)
    #    return B_vec/norm(B_vec)
    #end
    #____
    jacobian_matrix = stack(Interpolations.gradient(fields_itp, Rx, Rz))
    # Add zeros-column representing derivatives along the y-axis
    jacobian_matrix = [
        jacobian_matrix[1:6,1];;
        zeros(typeof(Rx), 6);; 
        jacobian_matrix[1:6,2]
        ]
    ∇B_vec = jacobian_matrix[1:3,:]
    ∇b = (∇B_vec - b * ∇B')*B_inv
    #
    # Calculate the Jacobian matrix of the ExB-drift
    #____ 
    # The following five lines is according to a @Benchmark test with 201^3
    # grid points 1.53 faster than the ForwardDiff.jacobian method (commented 
    #  below). The results of both methods where exactly the same.  
    #ForwardDiff.jacobian(R) do x
    #    B_vec = B_vec(x...)
    #    E_vec = E_vec(x...)
    #    return (E_vec × B_vec)/(norm(B_vec)^2)
    #end
    #____
    ∇E_vec = jacobian_matrix[4:6,:]
    skewE = skewsymmetric_matrix(E_vec)
    skewb = skewsymmetric_matrix(b)
    ∇ExB = (-skewb*∇E_vec + skewE*∇b - ExBdrift * ∇B')*B_inv

    # Electric field component parallel to the magnetic field
    Eparal = E_vec⋅b
    # Calculate drifts
    ∇Bdrift = q_inv*B_inv*μ*(b × ∇B)
    # Total time derivatives. Assumes ∂/∂t = 0,
    dBdt = vparal * b ⋅ ∇B + ExBdrift ⋅ ∇B
    dbdt = vparal * (∇b * b) + ∇b*ExBdrift
    dExBdt = vparal * (∇ExB * b) + ∇ExB*ExBdrift
    
    # Compute the perpendicular velcoity
    dRperpdt = ExBdrift + ∇Bdrift + q_inv*B_inv*m*b × (vparal*dbdt + dExBdt)
    #dRperpdt = b̂/B × (-E⃗ + μ/q*∇B + m/q*(vparal*db̂dt + dExBdt))  # old

    # Compute the acceleration 
    dvparaldt = (q*Eparal - μ*b⋅∇B)/m # along the magnetic field lines
    # With correction proposed by Birn et al., 2004:
    #dvparaldt = (q*Eparal - μ*b̂⋅∇B)/m + ExBdrift⋅db̂dt + ∇Bdrift⋅db̂dt
    #dRperpdt = b̂/B × (-c*E⃗ + μ*c/q * ∇B) #old

    # Compute the velocity
    dRdt = vparal*b + dRperpdt

    # Compute the pitch angle rate of change
    dbetadt = (1/vparal*dvparaldt - 1/2B*dBdt)*beta*(1-beta^2)

    # Compute the electron collision frequency
    coulomb_logarithm = 0.0
    eta = kf*2.91 * n[1] * coulomb_logarithm * T[1]^(-1.5)
    beta_friction = -eta*beta


    # Update the statevector
    return [dRdt; dvparaldt; dbetadt + beta_friction] 
end


function GCAPitchAngleDiffusion_lowmemory_2Dxz(u, params, _)
    Rx, Rz = u[1], u[3]    # Position of the gyrocentre
    beta = u[5]   # Cosine of pitch angle

    # Extract parameters
    # coulomb logarithm
    kd = params[3]
    # interpolation-object to get number density and gas temperature
    fields_itp = params[4]
    # Interpolate density to guiding centre location
    fields = fields_itp(Rx, Rz) 
    n = fields[7]
    T = fields[8]

    # Compute the electron collision frequency
    coulomb_logarithm = 20
    eta = kd * 2.91 * n[1] * coulomb_logarithm * T[1]^(-1.5)

    # Compute the diffusion coefficients
    dRdW = [0,0,0]
    dvparaldW = 0
    dbetadW = √(eta*(1-beta^2))

    return [dRdW; dvparaldW; dbetadW]
end

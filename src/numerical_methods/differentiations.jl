#-------------------------------------------------------------------------------
# Created 19.01.24. Code from back in  02.12.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 differentiations.jl
#
#-------------------------------------------------------------------------------
# Contains differentiation schemes
#-------------------------------------------------------------------------------


"""
    dervateupwind(
        field::Array{T, 3} where {T<:Real},
        dx   ::Vector{T} where {T<:Real},
        axis ::Tuple{Integer, Integer, Integer}
    )
Differentiates a 3D `field` with respect to a specified `axis` using the upwind
scheme. The grid size may be variable, hence given as the vector `dx`. End point
of result will be ill-calculated and the derivative will be defined at half grid
point higher than the input field.
"""
function derivateupwind(
    field::Array{T, 3} where {T<:Real},
    xx   ::Vector{T} where {T<:Real},
    yy   ::Vector{T} where {T<:Real},
    zz   ::Vector{T} where {T<:Real},
    ;
    wfp::DataType=typeof(field[1])
    )
    ni, nj, nk = size(field)

    if length(xx) > 1
        df = circshift(field, (-1,0,0)) - field
        # Pad the Δx array with a copy of the first value at the end. I.e. assume
        # periodic boundary conditions
        dx = circshift(xx, -1) - xx
        ddx = df ./ dx
    else
        ddx = zeros(wfp, ni,nj,nk)
    end

    if length(yy) > 1
        df = circshift(field, (0,-1,0)) - field
        dy = circshift(yy, -1) - yy
        ddy = df ./ dy'
    else
        ddy = zeros(wfp, ni, nj, nk)
    end

    if length(zz) > 1
        df = circshift(field, (0,0,-1)) - field
        dz = circshift(zz, -1) - zz
        # I don't know of a fast method for the 3rd dimension
        ddz = zeros(wfp, ni, nj, nk)
        for k = 1:nk
            ddz[:,:,k] = df[:,:,k] / dz[k]
        end
    end
    
    return ddx, ddy, ddz
end #function derivateUpwind


"""
    derivatecentral(field, dx)
First and last grid point are ill calculated.
"""
function derivatecentral(field::Vector{T} where {T<:Real},
                         dx
                         )
    return (circshift(field, -1) - circshift(field, 1))/2dx
end # function derivateCentral
#|
function derivatecentral(field      ::Array{T, 3} where {T<:Real},
                         gridSpacing::Real,
                         axis       ::Tuple{Int64, Int64, Int64}
                         )
    ax1 = -1 .* axis
    return (circshift(field, ax1) - circshift(field, axis))/2gridSpacing
end # function derivateCentral


"""
    derivate4thorder(field, dx)
First and last two grid point are ill calculated.
"""
function derivate4thorder(field      ::Array{T, 3} where {T<:Real},
                          gridSpacing::Real,
                          axis       ::Tuple{Int64, Int64, Int64}
                          )
    ax1 = -2 .* axis
    ax2 = -1 .* axis
    ax3 =  2 .* axis
    return (-circshift(field, ax1) + 8.0 .* circshift(field, ax2) - 
        8.0 .* circshift(field, axis) + circshift(field, ax3))/12.0gridSpacing
end # function derivate4thOrder


"""
    ∇(field, dx, dy, dz, scheme)
The gradient operator. Calucaletes the gradient of a 3- or 2-dimensional scalar
field and the Jacobian of a 3D vector field. Requires grid spacing on all three 
axis and the numerical scheme as arguments. The scheme is given as a function 
type, e.g. Schemes.derivateCentral.
"""
function ∇(field::Array{T, 4} where {T<:Real},
           dx   ::Real,
           dy   ::Real,
           dz   ::Real,
           scheme::Function
           ;
           wfp::DataType=typeof(field[1])
           )
    fx = field[1,:,:,:]
    fy = field[2,:,:,:]
    fz = field[3,:,:,:]
    #
    ∂fx∂x = scheme(fx, dx, (1,0,0))
    ∂fx∂y = scheme(fx, dy, (0,1,0))
    ∂fx∂z = scheme(fx, dz, (0,0,1))
    #
    ∂fy∂x = scheme(fy, dx, (1,0,0))
    ∂fy∂y = scheme(fy, dy, (0,1,0))
    ∂fy∂z = scheme(fy, dz, (0,0,1))
    #
    ∂fz∂x = scheme(fz, dx, (1,0,0))
    ∂fz∂y = scheme(fz, dy, (0,1,0))
    ∂fz∂z = scheme(fz, dz, (0,0,1))
    #
    _, nx, ny, nz = size(field)
    jacobian = zeros(wfp, 3, 3, nx, ny, nz)
    #
    jacobian[1, 1, :,:,:] = ∂fx∂x
    jacobian[1, 2, :,:,:] = ∂fx∂y
    jacobian[1, 3, :,:,:] = ∂fx∂z
    #
    jacobian[2, 1, :,:,:] = ∂fy∂x
    jacobian[2, 2, :,:,:] = ∂fy∂y
    jacobian[2, 3, :,:,:] = ∂fy∂z
    #
    jacobian[3, 1, :,:,:] = ∂fz∂x
    jacobian[3, 2, :,:,:] = ∂fz∂y
    jacobian[3, 3, :,:,:] = ∂fz∂z
    #
    return jacobian
end # function ∇ 
#|
function ∇(field::Array{T, 3} where {T<:Real},
           dx   ::Real,
           dy   ::Real,
           dz   ::Real,
           scheme::Function
           ;
           wfp::DataType=typeof(field[1])
           )
    ∂f∂x = scheme(field, dx, (1,0,0))
    ∂f∂y = scheme(field, dy, (0,1,0))
    ∂f∂z = scheme(field, dz, (0,0,1))
    nx, ny, nz = size(field)
    gradient = zeros(wfp, 3, nx, ny, nz)
    gradient[1, :, :, :] = ∂f∂x
    gradient[2, :, :, :] = ∂f∂y
    gradient[3, :, :, :] = ∂f∂z
    return gradient
end # function ∇ 
#|
function ∇(field::Array{T, 2} where {T<:Real},
           dx   ::Real,
           dy   ::Real,
           scheme::Function
           ;
           wfp::DataType=typeof(field[1])
           )
    dfdx = scheme(field, dx, (1,0,0))
    dfdy = scheme(field, dy, (0,1,0))
    nx, ny = size(field)
    gradient = zeros(wfp, 3, nx, ny)
    gradient[1, :, :] = ∂f∂x
    gradient[2, :, :] = ∂f∂y
    return gradient
end # function ∇ 
    
"""
Newer versions of the gradient, allowing for non-uniform structured grid.
"""
function ∇(field::Array{T, 4} where {T<:Real},
           xx   ::Vector{T} where {T<:Real},
           yy   ::Vector{T} where {T<:Real},
           zz   ::Vector{T} where {T<:Real},
           scheme::Function
           ;
           wfp::DataType=typeof(field[1])
           )
    #
    fx = field[:,:,:,1]
    fy = field[:,:,:,2]
    fz = field[:,:,:,3]
    #
    ∂fx∂x, ∂fx∂y, ∂fx∂z = scheme(fx, xx, yy, zz)
    ∂fy∂x, ∂fy∂y, ∂fy∂z = scheme(fy, xx, yy, zz)
    ∂fz∂x, ∂fz∂y, ∂fz∂z = scheme(fz, xx, yy, zz)
    #
    return [∂fx∂x, ∂fy∂x, ∂fz∂x, ∂fx∂y, ∂fy∂y, ∂fz∂y, ∂fx∂z, ∂fy∂z, ∂fz∂z]
end # function ∇
#|
function ∇(field::Array{T, 3} where {T<:Real},
           xx   ::Vector{T} where {T<:Real},
           yy   ::Vector{T} where {T<:Real},
           zz   ::Vector{T} where {T<:Real},
           scheme::Function
           ;
           wfp::DataType=typeof(field[1])
           )
        
    ∂f∂x, ∂f∂y, ∂f∂z = scheme(field, xx, yy, zz)
    return [∂f∂x, ∂f∂y, ∂f∂z]
end # function ∇


"""
    curlfield, gridsizes, derivscheme)
The curl operator. Calculates the curl of a 3-dimensional vector
field. Requires uniform grid spacing on all three axis. 
arguments. The scheme is given as a function type, e.g. Schemes.derivateCentral.
"""
function curl(field      ::Array{T, 4} where {T<:Real},
              gridsizes  ::Tuple{Real, Real, Real},
              derivscheme::Function
              ;
              wfp::DataType=typeof(field[1])
              )
    dx, dy, dz = gridsizes
    fx = field[1,:,:,:]
    fy = field[2,:,:,:]
    fz = field[3,:,:,:]
    derivx = (1,0,0)
    derivy = (0,1,0)
    derivz = (0,0,1)
    ∂fx∂y = derivscheme(fx, dy, derivy)
    ∂fx∂z = derivscheme(fx, dz, derivz)
    ∂fy∂x = derivscheme(fy, dx, derivx)
    ∂fy∂z = derivscheme(fy, dz, derivz)
    ∂fz∂x = derivscheme(fz, dx, derivx)
    ∂fz∂y = derivscheme(fz, dy, derivy)
    _, nx, ny, nz = size(field)
    result = zeros(wfp, 3, nx, ny, nz)
    result[1, :, :, :] = ∂fz∂y .- ∂fy∂z
    result[2, :, :, :] = ∂fx∂z .- ∂fz∂x
    result[3, :, :, :] = ∂fy∂x .- ∂fx∂y
    return result
end # functin curl

"""
    curl(field, dx, dy, dz, derivscheme)
The curl operator. Calculates the curl of a 3-dimensional vector
field. Requires grid spacing on all three axis (may be non-uniform)
The scheme is given as a function type, e.g. Schemes.derivateCentral.
"""
function curl(
    field      ::Array{T, 4} where {T<:Real},
    dx         ::Vector{T} where {T<:Real},
    dy         ::Vector{T} where {T<:Real},
    dz         ::Vector{T} where {T<:Real},
    derivscheme::Function
    ;
    wfp::DataType=typeof(field[1])
    )
    fx = field[1,:,:,:]
    fy = field[2,:,:,:]
    fz = field[3,:,:,:]
    derivx = (1,0,0)
    derivy = (0,1,0)
    derivz = (0,0,1)
    ∂fx∂y = derivscheme(fx, dy, derivy)
    ∂fx∂z = derivscheme(fx, dz, derivz)
    ∂fy∂x = derivscheme(fy, dx, derivx)
    ∂fy∂z = derivscheme(fy, dz, derivz)
    ∂fz∂x = derivscheme(fz, dx, derivx)
    ∂fz∂y = derivscheme(fz, dy, derivy)
    _, nx, ny, nz = size(field)
    result = zeros(wfp, 3, nx, ny, nz)
    result[1, :, :, :] = ∂fz∂y .- ∂fy∂z
    result[2, :, :, :] = ∂fx∂z .- ∂fz∂x
    result[3, :, :, :] = ∂fy∂x .- ∂fx∂y
    return result
end # function curl
#|
function curl(
    fx         ::Array{T, 3} where {T<:Real},
    fy         ::Array{T, 3} where {T<:Real},
    fz         ::Array{T, 3} where {T<:Real},
    xx         ::Vector{T} where {T<:Real},
    yy         ::Vector{T} where {T<:Real},
    zz         ::Vector{T} where {T<:Real},
    derivscheme::Function
    ;
    wfp::DataType=typeof(field[1])
    )
    ∂fx∂x, ∂fx∂y, ∂fx∂z = derivscheme(fx, xx, yy, zz)
    ∂fy∂x, ∂fy∂y, ∂fy∂z = derivscheme(fy, xx, yy, zz)
    ∂fz∂x, ∂fz∂y, ∂fz∂z = derivscheme(fz, xx, yy, zz)
    nx, ny, nz = size(fx)
    result = zeros(wfp, 3, nx, ny, nz)
    result[1, :, :, :] = ∂fz∂y .- ∂fy∂z
    result[2, :, :, :] = ∂fx∂z .- ∂fz∂x
    result[3, :, :, :] = ∂fy∂x .- ∂fx∂y
    return result
end # function curl


"""
    LinearAlgebra.cross(f, g)
Method for computing the cross product between two 3D vector-fields.
"""
function LinearAlgebra.cross(
    f::Array{T, 4} where {T<:Real},
    g::Array{T, 4} where {T<:Real}
    ;
    wfp::DataType=typeof(field[1])
    )
    _, ni, nj, nk = size(f)
    crossproduct = zeros(wfp, 3, ni, nj, nk)
    for i = 1:ni
        for j = 1:nj
            for k = 1:nk
                crossproduct[:,i,j,k] = f[:,i,j,k] × g[:,i,j,k]
            end
        end
    end
    return crossproduct
end # function cross

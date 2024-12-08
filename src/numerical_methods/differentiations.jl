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
    upwind1storder(
        field::Array{T, 3} where {T<:Real},
        dx   ::Vector{T} where {T<:Real},
        axis ::Tuple{Integer, Integer, Integer}
    )
Differentiates a 3D `field` with respect to a specified `axis` using the upwind
scheme. The grid size may be variable, hence given as the vector `dx`. End point
of result will be ill-calculated and the derivative will be defined at half grid
point higher than the input field.
"""
function upwind1storder(
    field::Array{T,3} where {T<:Real},
    xx::Vector{T} where {T<:Real},
    yy::Vector{T} where {T<:Real},
    zz::Vector{T} where {T<:Real},
    ;
    wfp::DataType=typeof(field[1])
)
    ni, nj, nk = size(field)

    if length(xx) > 1
        df = circshift(field, (-1, 0, 0)) - field
        # Pad the Δx array with a copy of the first value at the end. I.e. assume
        # periodic boundary conditions
        dx = circshift(xx, -1) - xx
        ddx = df ./ dx
    else
        ddx = zeros(wfp, ni, nj, nk)
    end

    if length(yy) > 1
        df = circshift(field, (0, -1, 0)) - field
        dy = circshift(yy, -1) - yy
        ddy = df ./ dy'
    else
        ddy = zeros(wfp, ni, nj, nk)
    end

    if length(zz) > 1
        df = circshift(field, (0, 0, -1)) - field
        dz = circshift(zz, -1) - zz
        # I don't know of a fast method for the 3rd dimension
        ddz = zeros(wfp, ni, nj, nk)
        for k = 1:nk
            ddz[:, :, k] = df[:, :, k] / dz[k]
        end
    end

    return ddx, ddy, ddz
end #function derivateUpwind


function central4thorder(
    field::Array{T,3} where {T<:Real},
    xx::Vector{T} where {T<:Real},
    yy::Vector{T} where {T<:Real},
    zz::Vector{T} where {T<:Real},
    ;
    wfp::DataType=typeof(field[1])
)
    ni, nj, nk = size(field)
    a = 0.08333333333333333f0
    b = 0.6666666666666666f0
    #
    if length(xx) > 1
        dx = diff(xx)[1]
        dfdx = (-a * circshift(field, (-2, 0, 0)) + b * circshift(field, (-1, 0, 0)) -
                b * circshift(field, (1, 0, 0)) + a * circshift(field, (2, 0, 0))) / dx
    else
        dfdx = zeros(wfp, ni, nj, nk)
    end
    if length(yy) > 1
        dy = diff(yy)[1]
        dfdy = (-a * circshift(field, (0, -2, 0)) + b * circshift(field, (0, -1, 0)) -
                b * circshift(field, (0, 1, 0)) + a * circshift(field, (0, 2, 0))) / 12.0dy
    else
        dfdy = zeros(wfp, ni, nj, nk)
    end
    if length(zz) > 1
        dz = diff(zz)[1]
        dfdz = (-a * circshift(field, (0, 0, -2)) + b * circshift(field, (0, 0, -1)) -
                b * circshift(field, (0, 0, 1)) + a * circshift(field, (0, 0, 2))) / dz
    else
        dfdz = zeros(wfp, ni, nj, nk)
    end
    return dfdx, dfdy, dfdz
end

"""
    bifrostscheme(
        field::Array{T, 3} where {T<:Real},
        xx   ::Vector{T} where {T<:Real},
        yy   ::Vector{T} where {T<:Real},
        zz   ::Vector{T} where {T<:Real},
        ;
        wfp::DataType=typeof(field[1])
    )
Differentiates a 3D `field` with respect to a specified `axis` using a 6th
order scheme, the same as the one used in the Bifrost code.

The grid size may be variable, hence given as the vector `dx`. End point
of result will be ill-calculated and the derivative will be defined at half grid
point higher than the input field.
"""
function bifrostscheme(
    field::Array{T,3} where {T<:Real},
    xx::Vector{T} where {T<:Real},
    yy::Vector{T} where {T<:Real},
    zz::Vector{T} where {T<:Real},
    ;
    wfp::DataType=typeof(field[1])
)
    ni, nj, nk = size(field)
    # 6th order upwind finite difference scheme coefficients:
    c = 3 / 640
    b = (-1 - 120c) / 24
    a = 1 - 3b - 5c
    # Differentiate along x-axis:
    if length(xx) > 1
        df1 = a * (-field + circshift(field, (-1, 0, 0)))
        df2 = b * (-circshift(field, (1, 0, 0)) + circshift(field, (-2, 0, 0)))
        df3 = c * (-circshift(field, (2, 0, 0)) + circshift(field, (-3, 0, 0)))
        # Pad the Δx array with a copy of the first value at the end. I.e. assume
        # periodic boundary conditions
        dx = circshift(xx, -1) - xx
        dfdx = (df1 + df2 + df3) ./ dx
    else
        dfdx = zeros(wfp, ni, nj, nk)
    end
    # Differentiate along y-axis:
    if length(yy) > 1
        df1 = a * (-field + circshift(field, (0, -1, 0)))
        df2 = b * (-circshift(field, (0, 1, 0)) + circshift(field, (0, -2, 0)))
        df3 = c * (-circshift(field, (0, 2, 0)) + circshift(field, (0, -3, 0)))
        # Pad the Δx array with a copy of the first value at the end. I.e. assume
        # periodic boundary conditions
        dy = circshift(yy, -1) - yy
        dfdy = (df1 + df2 + df3) ./ dy'
    else
        dfdy = zeros(wfp, ni, nj, nk)
    end
    # Differentiate along z-axis:
    if length(zz) > 1
        df1 = a * (-field + circshift(field, (0, 0, -1)))
        df2 = b * (-circshift(field, (0, 0, 1)) + circshift(field, (0, 0, -2)))
        df3 = c * (-circshift(field, (0, 0, 2)) + circshift(field, (0, 0, -3)))
        # I don't know of a fast method for the 3rd dimension
        dz = circshift(zz, -1) - zz
        dfdz = Array{wfp,3}(undef, ni, nj, nk)
        for k = 1:nk
            dfdz[:, :, k] = (df1[:, :, k] + df2[:, :, k] + df3[:, :, k]) / dz[k]
        end
    else
        dfdz = zeros(wfp, ni, nj, nk)
    end

    return dfdx, dfdy, dfdz
end #function bifrostscheme


"""
    derivatecentral(field, dx)
First and last grid point are ill calculated.
"""
function derivatecentral(field::Vector{T} where {T<:Real},
    dx
)
    return (circshift(field, -1) - circshift(field, 1)) / 2dx
end # function derivateCentral
#|
function derivatecentral(field::Array{T,3} where {T<:Real},
    gridSpacing::Real,
    axis::Tuple{Int64,Int64,Int64}
)
    ax1 = -1 .* axis
    return (circshift(field, ax1) - circshift(field, axis)) / 2gridSpacing
end # function derivateCentral


"""
    derivate4thorder(field, dx)
First and last two grid point are ill calculated.
"""
function derivate4thorder(field::Array{T,3} where {T<:Real},
    gridSpacing::Real,
    axis::Tuple{Int64,Int64,Int64}
)
    ax1 = -2 .* axis
    ax2 = -1 .* axis
    ax3 = 2 .* axis
    return (-circshift(field, ax1) + 8.0 .* circshift(field, ax2) -
            8.0 .* circshift(field, axis) + circshift(field, ax3)) / 12.0gridSpacing
end # function derivate4thOrder


"""
    ∇(field, dx, dy, dz, scheme)
The gradient operator. Calucaletes the gradient of a 3- or 2-dimensional scalar
field and the Jacobian of a 3D vector field. Requires grid spacing on all three 
axis and the numerical scheme as arguments. The scheme is given as a function 
type, e.g. Schemes.derivateCentral.
"""
function ∇(field::Array{T,4} where {T<:Real},
    dx::Real,
    dy::Real,
    dz::Real,
    scheme::Function
    ;
    wfp::DataType=typeof(field[1])
)
    fx = field[1, :, :, :]
    fy = field[2, :, :, :]
    fz = field[3, :, :, :]
    #
    ∂fx∂x = scheme(fx, dx, (1, 0, 0))
    ∂fx∂y = scheme(fx, dy, (0, 1, 0))
    ∂fx∂z = scheme(fx, dz, (0, 0, 1))
    #
    ∂fy∂x = scheme(fy, dx, (1, 0, 0))
    ∂fy∂y = scheme(fy, dy, (0, 1, 0))
    ∂fy∂z = scheme(fy, dz, (0, 0, 1))
    #
    ∂fz∂x = scheme(fz, dx, (1, 0, 0))
    ∂fz∂y = scheme(fz, dy, (0, 1, 0))
    ∂fz∂z = scheme(fz, dz, (0, 0, 1))
    #
    _, nx, ny, nz = size(field)
    jacobian = zeros(wfp, 3, 3, nx, ny, nz)
    #
    jacobian[1, 1, :, :, :] = ∂fx∂x
    jacobian[1, 2, :, :, :] = ∂fx∂y
    jacobian[1, 3, :, :, :] = ∂fx∂z
    #
    jacobian[2, 1, :, :, :] = ∂fy∂x
    jacobian[2, 2, :, :, :] = ∂fy∂y
    jacobian[2, 3, :, :, :] = ∂fy∂z
    #
    jacobian[3, 1, :, :, :] = ∂fz∂x
    jacobian[3, 2, :, :, :] = ∂fz∂y
    jacobian[3, 3, :, :, :] = ∂fz∂z
    #
    return jacobian
end # function ∇ 
#|
function ∇(field::Array{T,3} where {T<:Real},
    dx::Real,
    dy::Real,
    dz::Real,
    scheme::Function
    ;
    wfp::DataType=typeof(field[1])
)
    ∂f∂x = scheme(field, dx, (1, 0, 0))
    ∂f∂y = scheme(field, dy, (0, 1, 0))
    ∂f∂z = scheme(field, dz, (0, 0, 1))
    nx, ny, nz = size(field)
    gradient = zeros(wfp, 3, nx, ny, nz)
    gradient[1, :, :, :] = ∂f∂x
    gradient[2, :, :, :] = ∂f∂y
    gradient[3, :, :, :] = ∂f∂z
    return gradient
end # function ∇ 
#|
function ∇(field::Array{T,2} where {T<:Real},
    dx::Real,
    dy::Real,
    scheme::Function
    ;
    wfp::DataType=typeof(field[1])
)
    dfdx = scheme(field, dx, (1, 0, 0))
    dfdy = scheme(field, dy, (0, 1, 0))
    nx, ny = size(field)
    gradient = zeros(wfp, 3, nx, ny)
    gradient[1, :, :] = ∂f∂x
    gradient[2, :, :] = ∂f∂y
    return gradient
end # function ∇ 

"""
NEW VERSIONS of the gradient and curl, allowing for non-uniform structured grid.
"""
function ∇(field::Array{T,4} where {T<:Real},
    xx::Vector{T} where {T<:Real},
    yy::Vector{T} where {T<:Real},
    zz::Vector{T} where {T<:Real},
    scheme::Function
    ;
    wfp::DataType=typeof(field[1])
)
    #
    fx = field[:, :, :, 1]
    fy = field[:, :, :, 2]
    fz = field[:, :, :, 3]
    #
    ∂fx∂x, ∂fx∂y, ∂fx∂z = scheme(fx, xx, yy, zz)
    ∂fy∂x, ∂fy∂y, ∂fy∂z = scheme(fy, xx, yy, zz)
    ∂fz∂x, ∂fz∂y, ∂fz∂z = scheme(fz, xx, yy, zz)
    #
    return [∂fx∂x, ∂fy∂x, ∂fz∂x, ∂fx∂y, ∂fy∂y, ∂fz∂y, ∂fx∂z, ∂fy∂z, ∂fz∂z]
end # function ∇
#|
function ∇(field::Array{T,3} where {T<:Real},
    xx::Vector{T} where {T<:Real},
    yy::Vector{T} where {T<:Real},
    zz::Vector{T} where {T<:Real},
    scheme::Function
    ;
    wfp::DataType=typeof(field[1])
)

    ∂f∂x, ∂f∂y, ∂f∂z = scheme(field, xx, yy, zz)
    return [∂f∂x, ∂f∂y, ∂f∂z]
end # function ∇


"""
    ∇(
        coords::Vector{<:Real},
        itpvec::Vector{AbstractInterpolation},
        )
Compute the local gradient of the strength of a vector field given by an vector
of interpolation objects.
"""
function ∇(
    coords::Vector{<:Real},
    itpvec::Vector{<:AbstractInterpolation},
)
    res = ForwardDiff.gradient(coords) do x
        norm([itp(x...) for itp in itpvec])
    end
end
function ∇(
    coords::Vector{<:Real},
    itp::AbstractInterpolation,
)
    res = ForwardDiff.gradient(coords) do x
        norm(itp(x...)[1:3])
    end
end



"""
    curl(field, gridsizes, derivscheme)
The curl operator. Calculates the curl of a 3-dimensional vector
field. Requires uniform grid spacing on all three axis. 
arguments. The scheme is given as a function type, e.g. Schemes.derivateCentral.
"""
function curl(field::Array{T,4} where {T<:Real},
    gridsizes::Tuple{Real,Real,Real},
    derivscheme::Function
    ;
    wfp::DataType=typeof(field[1])
)
    dx, dy, dz = gridsizes
    fx = field[1, :, :, :]
    fy = field[2, :, :, :]
    fz = field[3, :, :, :]
    derivx = (1, 0, 0)
    derivy = (0, 1, 0)
    derivz = (0, 0, 1)
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
    field::Array{T,4} where {T<:Real},
    xx::Vector{T} where {T<:Real},
    yy::Vector{T} where {T<:Real},
    zz::Vector{T} where {T<:Real},
    derivscheme::Function
    ;
    wfp::DataType=typeof(field[1])
)
    fx = field[1, :, :, :]
    fy = field[2, :, :, :]
    fz = field[3, :, :, :]
    return curl(fx, fy, fz, xx, yy, zz, derivscheme; wfp=wfp)
end
#|
function curl(
    fx::Array{T,3} where {T<:Real},
    fy::Array{T,3} where {T<:Real},
    fz::Array{T,3} where {T<:Real},
    xx::Vector{T} where {T<:Real},
    yy::Vector{T} where {T<:Real},
    zz::Vector{T} where {T<:Real},
    derivscheme::Function
    ;
    wfp::DataType=eltype(fx)
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
    f::Array{T,4} where {T<:Real},
    g::Array{T,4} where {T<:Real}
    ;
    wfp::DataType=typeof(field[1])
)
    _, ni, nj, nk = size(f)
    crossproduct = zeros(wfp, 3, ni, nj, nk)
    for i = 1:ni
        for j = 1:nj
            for k = 1:nk
                crossproduct[:, i, j, k] = f[:, i, j, k] × g[:, i, j, k]
            end
        end
    end
    return crossproduct
end # function cross


"""
    ∇(
        coords::Vector{<:Real},
        itpvec::Vector{AbstractInterpolation},
        )
Compute the local gradient of the strength of a vector field given by an vector
of interpolation objects.
"""
function ∇(
    coords::Vector{<:Real},
    itpvec::Vector{<:AbstractInterpolation},
    )
    res = ForwardDiff.gradient(coords) do x
        norm([itp(x...) for itp in itpvec])
    end
end


"""
    jacobian(
        coords::Vector{<:Real},
        itpvec::Vector{AbstractInterpolation},
        )
Compute the local Jacobian of a vector field given by an vector
of interpolation objects.
"""
function jacobian(
    coords::Vector{<:Real},
    itpvec::Vector{<:AbstractInterpolation},
)
    res = ForwardDiff.jacobian(coords) do x
        [itp(x...) for itp in itpvec]
    end
end


function local_grid_fdm(
    fdm::FiniteDifferences.FiniteDifferenceMethod,
    A::Vector{<:Real},
    idx::Integer,
    dx::Real,
)
    diff = 0.0
    for (i, coef) in zip(fdm.grid, fdm.coefs)
        #diff += coef*A[idx + i]
        diff += coef * (circshift(A, -i)[idx])
    end
    return diff / dx
end


function local_grid_fdm(
    fdm::FiniteDifferences.FiniteDifferenceMethod,
    A::Array{<:Real,3},
    x::Real,
    y::Real,
    z::Real,
    xaxis::Vector{<:Real},
    yaxis::Vector{<:Real},
    zaxis::Vector{<:Real},
)
    wfp = Float64
    i = locate_cell(xaxis, x)
    j = locate_cell(yaxis, y)
    k = locate_cell(zaxis, z)
    dx = diff(xaxis)[1]
    dy = diff(yaxis)[1]
    dz = diff(zaxis)[1]
    #
    gradx = Array{wfp,3}(undef, 2, 2, 2)
    gradx[1] = local_grid_fdm(fdm, A[:, j, k], i, dx)
    gradx[2] = local_grid_fdm(fdm, A[:, j, k], i + 1, dx)
    gradx[3] = local_grid_fdm(fdm, A[:, j+1, k], i, dx)
    gradx[4] = local_grid_fdm(fdm, A[:, j+1, k], i + 1, dx)
    gradx[5] = local_grid_fdm(fdm, A[:, j, k+1], i, dx)
    gradx[6] = local_grid_fdm(fdm, A[:, j, k+1], i + 1, dx)
    gradx[7] = local_grid_fdm(fdm, A[:, j+1, k+1], i, dx)
    gradx[8] = local_grid_fdm(fdm, A[:, j+1, k+1], i + 1, dx)
    #
    grady = Array{wfp,3}(undef, 2, 2, 2)
    grady[1] = local_grid_fdm(fdm, A[i, :, k], j, dy)
    grady[2] = local_grid_fdm(fdm, A[i+1, :, k], j, dy)
    grady[3] = local_grid_fdm(fdm, A[i, :, k], j + 1, dy)
    grady[4] = local_grid_fdm(fdm, A[i+1, :, k], j + 1, dy)
    grady[5] = local_grid_fdm(fdm, A[i, :, k+1], j, dy)
    grady[6] = local_grid_fdm(fdm, A[i+1, :, k+1], j, dy)
    grady[7] = local_grid_fdm(fdm, A[i, :, k+1], j + 1, dy)
    grady[8] = local_grid_fdm(fdm, A[i+1, :, k+1], j + 1, dy)
    #
    gradz = Array{wfp,3}(undef, 2, 2, 2)
    gradz[1] = local_grid_fdm(fdm, A[i, j, :], k, dz)
    gradz[2] = local_grid_fdm(fdm, A[i+1, j, :], k, dz)
    gradz[3] = local_grid_fdm(fdm, A[i, j+1, :], k, dz)
    gradz[4] = local_grid_fdm(fdm, A[i+1, j+1, :], k, dz)
    gradz[5] = local_grid_fdm(fdm, A[i, j, :], k + 1, dz)
    gradz[6] = local_grid_fdm(fdm, A[i+1, j, :], k + 1, dz)
    gradz[7] = local_grid_fdm(fdm, A[i, j+1, :], k + 1, dz)
    gradz[8] = local_grid_fdm(fdm, A[i+1, j+1, :], k + 1, dz)
    #
    xx = xaxis[i]:dx:xaxis[i+1]
    yy = yaxis[j]:dy:yaxis[j+1]
    zz = zaxis[k]:dz:zaxis[k+1]
    axes = (xx, yy, zz)
    itpx = linear_interpolation(axes, gradx)
    itpy = linear_interpolation(axes, grady)
    itpz = linear_interpolation(axes, gradz)
    return [itpx(x, y, z), itpy(x, y, z), itpz(x, y, z)]
end

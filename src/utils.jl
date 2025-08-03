#-------------------------------------------------------------------------------
# Created 17.01.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                Utilities.jl
#
#-------------------------------------------------------------------------------
# Module containing the utility functions.
#-------------------------------------------------------------------------------

"""
    ElectromagneticFieldInterpolator{T<:AbstractVector}(interpolator::T)
A convenience struct for a vector of interpolators representing the components
of the electromagnetic field. Following the convention of storing the
vector of components that are stored as a JLD2 object [bx, by, bz, ex, ey, ez].
Calling an instance of this struct with spatial coordinates returns the electric
and magnetic fields at that point as two `SVector{3}` vectors.
"""
struct ElectromagneticFieldInterpolator{T<:AbstractVector}
    interpolator::T
end
function (self::ElectromagneticFieldInterpolator{T})(x, y, z, t) where {T}
    ex = self.interpolator[4](x, y, z, t)
    ey = self.interpolator[5](x, y, z, t)
    ez = self.interpolator[6](x, y, z, t)
    bx = self.interpolator[1](x, y, z, t)
    by = self.interpolator[2](x, y, z, t)
    bz = self.interpolator[3](x, y, z, t)
    return SVector{3}(ex, ey, ez), SVector{3}(bx, by, bz)
end
function (self::ElectromagneticFieldInterpolator{T})(x, y, z) where {T}
    ex = self.interpolator[4](x, y, z)
    ey = self.interpolator[5](x, y, z)
    ez = self.interpolator[6](x, y, z)
    bx = self.interpolator[1](x, y, z)
    by = self.interpolator[2](x, y, z)
    bz = self.interpolator[3](x, y, z)
    return SVector{3}(ex, ey, ez), SVector{3}(bx, by, bz)
end
function (self::ElectromagneticFieldInterpolator{T})(x, y) where {T}
    ex = self.interpolator[4](x, y)
    ey = self.interpolator[5](x, y)
    ez = self.interpolator[6](x, y)
    bx = self.interpolator[1](x, y)
    by = self.interpolator[2](x, y)
    bz = self.interpolator[3](x, y)
    return SVector{3}(ex, ey, ez), SVector{3}(bx, by, bz)
end

"""
    StaticInterpolation
An interpolator that calls `(x, y, z)` if called with `(x, y, z, t)`.
## Fields
- itp::AbstractInterpolation
"""
struct StaticInterpolation{T}
    itp::T
end
function (self::StaticInterpolation)(x::Real, y::Real, z::Real, t::Real)
    self.itp(x, y, z)
end
function (self::StaticInterpolation)(x::Real, y::Real, z::Real)
    self.itp(x, y, z)
end
function Base.ndims(itp::StaticInterpolation)
    return ndims(itp.itp)
end
function Base.size(itp::StaticInterpolation)
    return size(itp.itp)
end
"""
    XZInterpolation
An interpolator that calls `(x, z, t)` if called with `(x, y, z, t)`.
## Fields
- itp::AbstractInterpolation
"""
struct XZInterpolation{T}
    itp::T
end
function (self::XZInterpolation)(x::Real, y::Real, z::Real, t::Real)
    self.itp(x, z, t)
end
function (self::XZInterpolation)(x::Real, z::Real, t::Real)
    self.itp(x, z, t)
end
function Base.ndims(itp::XZInterpolation)
    return ndims(itp.itp)
end
function Base.size(itp::XZInterpolation)
    return size(itp.itp)
end

"""
    StaticXZInterpolation
An interpolator that calls `(x, z)` if called with `(x, y, z, t)`.
## Fields
- itp::AbstractInterpolation
"""
struct StaticXZInterpolation{T}
    itp::T
end
function (self::StaticXZInterpolation)(x::Real, y::Real, z::Real, t::Real)
    self.itp(x, z)
end
function (self::StaticXZInterpolation)(x::Real, y::Real, z::Real)
    self.itp(x, z)
end
function (self::StaticXZInterpolation)(x::Real, z::Real)
    self.itp(x, z)
end
function Base.ndims(itp::StaticXZInterpolation)
    return ndims(itp.itp)
end
function Base.size(itp::StaticXZInterpolation)
    return size(itp.itp)
end

"""
    Base.dropdims(arr::AbstractArray)
Extend `dropdims` to find and drop all axes with one point when `dims` keyword 
is excluded.
"""
function Base.dropdims(arr::AbstractArray)
    meshsize = size(arr)
    idxs = findall(x -> x == 1, meshsize)
    for dim in idxs
        arr = dropdims(arr, dims=dim)
    end
    return arr
end


"""
    Base.dropdims(axes::Tuple{Vararg{Vector})
Extend `dropdims` to find and drop all axes with one point when `dims` keyword 
is excluded.
"""
function Base.dropdims(axes::Tuple{Vararg{Vector}})
    mask = findall(x -> length(x) != 1, axes)
    return axes[mask]
end


#____/\_____/\_________________________________________________________________

"""
Convenience structure for storing the maximum value of some quantity during
a `Differentialequations` simulation.
"""
mutable struct MaxValue{T1,T2,T3}
    max::T1
    max_u::T2
    max_t::T3
end


"""
Convenience structure for storing the minimum value of some quantity during
a `Differentialequations` simulation.
"""
mutable struct MinValue{T1,T2,T3}
    min::T1
    min_u::T2
    min_t::T3
end

function findslice(axis, lim)
    if lim[1] > axis[end] || lim[2] < axis[1]
        error("Limit did not match any points in the mesh. Wrong units?")
    end
    i = findlast(axis -> axis < lim[1], axis)
    j = findfirst(axis -> axis > lim[2], axis)
    if isnothing(j)
        j = length(axis)
    end
    if isnothing(i)
        i = 1
    end
    return range(i, j)
end
"""
    findslice(x, y, z, xlim, ylim, zlim)
    findslice(x, z, xlim, zlim)
Find the ranges which slices the `BifrostExperiment` according to the given
limits.
"""
function findslice(
    x::AbstractVector,
    y::AbstractVector,
    z::AbstractVector,
    xlim::Any,
    ylim::Any,
    zlim::Any
    ;
    verbose=true
)
    if isnothing(xlim)
        slicex = 1:length(x)
    else
        slicex = findslice(x, xlim)
    end
    if isnothing(zlim)
        slicez = 1:length(z)
    else
        slicez = findslice(z, zlim)
    end
    if isnothing(ylim)
        slicey = 1:length(y)
    else
        slicey = findslice(y, ylim)
    end
    if verbose
        @info "Range of x-axis: $(slicex[1]:slicex[end])"
        @info "Range of y-axis: $(slicey[1]:slicey[end])"
        @info "Range of z-axis: $(slicez[1]:slicez[end])"
    end
    return slicex, slicey, slicez
end
function findslice(
    x::AbstractVector,
    z::AbstractVector,
    xlim::Any,
    zlim::Any
    ;
    verbose=true
)
    if isnothing(xlim)
        slicex = 1:length(x)
    else
        slicex = findslice(x, xlim)
    end
    if isnothing(zlim)
        slicez = 1:length(z)
    else
        slicez = findslice(z, zlim)
    end
    if verbose
        @info "Range of x-axis: $(slicex[1]:slicex[end])"
        @info "Range of z-axis: $(slicez[1]:slicez[end])"
    end
    return slicex, slicez
end

#----------------------------------------#
# Vector potential generation            #
# Mesh generation from analytical fields #
#-------------------------------------------------------------------------------

"""
    create3Daxes(
        (x0, y0, z0)::Tuple{Real,Real,Real},
        (xf, yf, zf)::Tuple{Real,Real,Real},
        (nx, ny, nz)::Tuple{Integer,Integer,Integer}
    )
Create axes from axis limits and number of grid points.
"""
function create3Daxes(
    (x0, y0, z0)::Tuple{Real,Real,Real},
    (xf, yf, zf)::Tuple{Real,Real,Real},
    (nx, ny, nz)::Tuple{Integer,Integer,Integer}
)
    xx = range(x0, xf, length=nx)
    dx = xx[2] - xx[1]
    yy = range(y0, yf, length=ny)
    dy = yy[2] - yy[1]
    # Account for single point in the z-axis
    if nz == 1
        zz = [z0]
        dz = 0.0
    else
        zz = range(z0, zf, length=nz)
        dz = zz[2] - zz[1]
    end
    return xx, yy, zz, dx, dy, dz
end # function createaxes


"""
    discretise(
        func::Function,
        xx  ::Vector{T} where {T<:Real},
        yy  ::Vector{T} where {T<:Real},
        zz  ::Vector{T} where {T<:Real},
        args...
        )
Discretise a `func`tion onto a mesh with grid points given by
`xx`, `yy` and `zz`.
"""
function discretise(
    func::Function,
    xx::AbstractVector,
    yy::AbstractVector,
    zz::AbstractVector,
    args...
)
    f111 = func(xx[1], yy[1], zz[1])
    field = Array{Vector{eltype(f111)},3}(
        undef, length(xx), length(yy), length(zz)
    )
    for i in eachindex(xx)
        for j in eachindex(yy)
            for k in eachindex(zz)
                field[i, j, k] = func(xx[i], yy[j], zz[k])
            end
        end
    end
    return field
end # function discretise

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
of the electromagnetic field. Assumes the vector of components stored as
[bx, by, bz, ex, ey, ez].
Calling an instance of this struct with spatial coordinates returns the electric
and magnetic fields at that point as two `SVector{3}` vectors.
"""
struct ElectromagneticFieldInterpolator{T}
    interpolator::T
end
function (self::ElectromagneticFieldInterpolator{T})(
    x, y, z, t
) where {T<:AbstractVector}
    ex = self.interpolator[4](x, y, z, t)
    ey = self.interpolator[5](x, y, z, t)
    ez = self.interpolator[6](x, y, z, t)
    bx = self.interpolator[1](x, y, z, t)
    by = self.interpolator[2](x, y, z, t)
    bz = self.interpolator[3](x, y, z, t)
    return SVector{3}(ex, ey, ez), SVector{3}(bx, by, bz)
end
function (self::ElectromagneticFieldInterpolator{T})(
    x, y, z
) where {T<:AbstractVector}
    ex = self.interpolator[4](x, y, z)
    ey = self.interpolator[5](x, y, z)
    ez = self.interpolator[6](x, y, z)
    bx = self.interpolator[1](x, y, z)
    by = self.interpolator[2](x, y, z)
    bz = self.interpolator[3](x, y, z)
    return SVector{3}(ex, ey, ez), SVector{3}(bx, by, bz)
end
function (self::ElectromagneticFieldInterpolator{T})(
    x, y
) where {T<:AbstractVector}
    ex = self.interpolator[4](x, y)
    ey = self.interpolator[5](x, y)
    ez = self.interpolator[6](x, y)
    bx = self.interpolator[1](x, y)
    by = self.interpolator[2](x, y)
    bz = self.interpolator[3](x, y)
    return SVector{3}(ex, ey, ez), SVector{3}(bx, by, bz)
end
function (self::ElectromagneticFieldInterpolator{T})(
    x, y, z, t
) where {T}
    bx, by, bz, ex, ey, ez = self.interpolator(x, y, z, t)
    return SVector(ex, ey, ez), SVector(bx, by, bz)
end
function (self::ElectromagneticFieldInterpolator{T})(
    x, y, z
) where {T}
    bx, by, bz, ex, ey, ez = self.interpolator(x, y, z)
    return SVector(ex, ey, ez), SVector(bx, by, bz)
end
function (self::ElectromagneticFieldInterpolator{T})(
    x, y
) where {T}
    bx, by, bz, ex, ey, ez = self.interpolator(x, y)
    return SVector(ex, ey, ez), SVector(bx, by, bz)
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
function Base.maximum(itp::StaticInterpolation)
    return maximum(Interpolations.coefficients(itp.itp.itp))
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
function Base.maximum(itp::XZInterpolation)
    return maximum(Interpolations.coefficients(itp.itp.itp))
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
function Base.maximum(itp::StaticXZInterpolation)
    return maximum(Interpolations.coefficients(itp.itp.itp))
end

#____/\_____/\_________________________________________________________________
#
# Various utility functions
#
"""
    squeeze(arr::AbstractArray)
Find and drop all axes with one point.
"""
function squeeze(arr::AbstractArray)
    singular_axes = findall(x -> x == 1, size(arr))
    while !(isempty(singular_axes))
        arr = dropdims(arr, dims=first(singular_axes))
        singular_axes = findall(x -> x == 1, size(arr))
    end
    return arr
end


"""
    squeeze(axes::Tuple{Vararg{Vector})
Find and drop all axes with one point.
"""
function squeeze(axes::Tuple{Vararg{Vector}})
    mask = findall(x -> length(x) != 1, axes)
    return axes[mask]
end

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


function logtimed(time_taken, msg)
    @info """$(msg):
    Elapsed time: $(time_taken.time) seconds
    Bytes allocated: $(@sprintf "%.1E" time_taken.bytes)"""
end

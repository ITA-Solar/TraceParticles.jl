#-------------------------------------------------------------------------------
# Created 17.01.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                utils.jl
#
#-------------------------------------------------------------------------------
# Module containing the utility functions.
#-------------------------------------------------------------------------------

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

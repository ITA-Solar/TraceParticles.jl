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
    itpvec::Vector
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

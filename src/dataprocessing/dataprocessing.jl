#-------------------------------------------------------------------------------
# Created 29.01.24
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 dataprocessing.jl
#
#-------------------------------------------------------------------------------
# Contains functions for processing results.
#-------------------------------------------------------------------------------
"""
    getinitialstate(
        u::Vector{<:ODESolution},
        )
Return the initial state (or initial condition) of a vector of `ODEsolution`'s.
"""
function getinitialstate(
    u::Vector{<:ODESolution},
    )
    npart = length(u)[1]
    ndof = length(u[1].u[1])
    RealT = typeof(u[1].u[1][1])
    u0 = Array{RealT}(undef, ndof, npart)
    for i = 1:npart
        u0[:,i] .= u[i].u[1]
    end
    return u0
end


"""
    getinitialstate(
        u::Array{<:Real, 3},:
        )
Return the initial state (or initial condition) of a solution stored as a
3D-array.
"""
function getinitialstate(
    u::Array{<:Real, 3},
    )
    return u[1,:,:]
end


"""
    getfinalstate(
        u::Vector{<:ODESolution},
        )
Return the final state of a vector of `ODEsolution`'s.
"""
function getfinalstate(
    u::Vector{<:ODESolution},
    )
    npart = length(u)[1]
    ndof = length(u[1].u[1])
    RealT = typeof(u[1].u[1][1])
    uf = Array{RealT}(undef, ndof, npart)
    for i = 1:npart
        uf[:,i] .= u[i].u[end]
    end
    return uf
end


"""
    getfinalstate(
        u::Array{<:Real, 3},:
        )
Return the final state of a solution stored as a 3D-array.
"""
function getfinalstate(
    u::Array{<:Real, 3},
    )
    return u[end,:,:]
end


"""
    getfinaltime(
        u::Vector{<:ODESolution},
        )
Return the final times of a vector of `ODEsolution`'s.
"""
function getfinaltime(
    u::Vector{<:ODESolution},
    )
    return [u.t[end] for u in u.u]
end


function theoretical_energygain_2Dxz(
    solution::ODESolution
    )
    error("Incomplete function. Need a way to calculate the electric potential.")
    nsteps = length(solution.t)
    charge = solution.prob.p[1]
    x0 = solution.u[1][1]
    z0 = solution.u[1][3]
    E_0 = norm([solution.prob.p[4][i](x0, z0) for i = 4:6])
    E_gain = Vector{eltype(E_0)}(undef, nsteps - 1)
    for i = 1:nsteps-1
        x = solution.u[i][1]
        z = solution.u[i][3]
        E_gain = norm([solution.prob.p[4][i](x, z) for i = 4:6]) - E_0
    end
    return charge*E_gain
end

function avgnumsteps(
    solution::Vector{Tuple{
        Any, Any, Any, Vector, Vector, Vector, Vector, Any, Any
        }}
    )
    sum = 0
    for particle in solution
        sum += particle[8]
    end
    return sum/length(solution)
end


function findinittemp(
    solution::Vector{Tuple{
        Any, Any, Any, Vector, Vector, Vector, Vector, Any, Any
        }},
    tempitp::AbstractInterpolation
    )
    nofparticles = length(solution)
    T = Vector{typeof(solution[1][1][1])}(undef, nofparticles)
    for (particle, i) in zip(solution, 1:nofparticles)
        T[i] = tempitp(particle[1][1], particle[1][3])
    end
    return T
end

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
Convenience struct for storing user data in the parameters of a
`DifferentialEquations` problem.

For now contains only a return message `retmsg` to be edited by user-defined
callbacks that terminate the simulation.
"""
mutable struct UserData
    retmsg::String
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
    return range(i,j)
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


"""
    normal3Donlyz((x0, y0, z0), 
                  (xf, yf, zf), 
                  (nx, ny, nz), 
                  (μx, μy), 
                  (σx, σy),
                  amplitude
                  )
According to given spatial domain, resolution, expectationvalue, standard
deviation and amplitude, creates a vector field who's z-component is normally
distributed in x and y according to the formula fz(i, j) = amplitude *
fx(i)fy(j), where fx and fy are normal distributions in x, and y with
expectation value and std equal to μx, μy, σx, σy, respetively. 
"""
function normal3Donlyz(
    (x0, y0, z0)::Tuple{Real,Real,Real},
    (xf, yf, zf)::Tuple{Real,Real,Real},
    (nx, ny, nz)::Tuple{Integer,Integer,Integer},
    (μx, μy)::Tuple{Real,Real},
    (σx, σy)::Tuple{Real,Real},
    amplitude::Real=1.0
)
    # Create spatial axes and find the grid sizes
    xx, yy, zz, dx, dy, dz = create3Daxes(
        (x0, y0, z0),
        (xf, yf, zf),
        (nx, ny, nz)
    )
    # Initialise the vector field
    ndims = 3
    A = zeros(ndims, nx, ny, nz)
    # Evaluate the z-component of the vecor field to be normally distributed in
    # the x and y dimensions.
    gaussx = normaldistr(xx, μx, σx)
    gaussy = normaldistr(yy, μy, σy)
    for i = 1:nx
        for j = 1:ny
            A[3, i, j, :] .= amplitude * gaussx[i] * gaussy[j]
        end
    end
    return (xx, yy, zz), (dx, dy, dz), A
end # function normal3donlyz


#-------------------------#
# Particle initialisation #
#-------------------------------------------------------------------------------
"""
    initparticlesuniform(numparticles, pos0, posf, vel0, velf, seed)
Distribute positions and velocities uniformly according to limits.
"""
function initparticlesuniform(
    numparticles::Integer,
    pos0::Vector{T} where {T<:Real},
    posf::Vector{T} where {T<:Real},
    vel0::Vector{T} where {T<:Real},
    velf::Vector{T} where {T<:Real},
    seed::Integer=0  # random-seed
)
    numdims = 3
    spatialextent = posf .- pos0
    velocityrange = velf .- vel0
    r = rand(MersenneTwister(seed),
        typeof(pos0[1]), (Int64(2numdims), Int64(numparticles)))
    # 2 for both position and
    # velocity 
    positions = pos0 .+ (spatialextent .* r[1:numdims, :])
    velocities = vel0 .+ (velocityrange .* r[4:2numdims, :])
    return positions, velocities
end # function initparticlesuniform


"""
    initparticlesmaxwellian(numparticles, pos0, posf, vel0, velf, seed)
Distribute positions uniformly and velocities Maxwellian according to limits.
"""
function initparticlesmaxwellian(
    numparticles::Integer,
    pos0::Vector{T} where {T<:Real},
    posf::Vector{T} where {T<:Real},
    temperature::Real, # temperature of the Maxwellian distribution
    mass::Real, # mass of particles
    seed::Integer=0  # random-seed
    ;
    wfp::DataType=typeof(pos0[1])
)
    numdims = 3
    #
    # Velocities
    σ = √(k_B * temperature / mass) # Standard deviation of velocity
    # components 
    μ = 0.0 # Expectation-value of velocity distributions. 
    # Generate velocities from a normal distribution
    velocities = μ .+ σ .* randn(MersenneTwister(seed),
        wfp, (numdims, numparticles))
    #
    # Posistions: Generate from a uniform distribution
    spatialextent = posf .- pos0
    positions = pos0 .+
                (spatialextent .* rand(MersenneTwister(seed),
        wfp, (numdims, numparticles)))
    #
    return positions, velocities
end # function initparticlesmaxwellian
# Created 22.01.24
# Author: e.s.oyre@astro.uio.no
#
#           statistics.jl
#-------------------------------------------------------------------------------
# Contains statistical distributions, random number generators and samplers.
#-------------------------------------------------------------------------------


#---------------#
# Distributions #
#---------------#---------------------------------------------------------------
"""
    normaldistr(
        x::Vector{T} where {T<:Real},
        μ::Real,
        σ::Real
        )

Function returning the values of `x` on a 1D normalised normal distribution with
expectation value `μ` and standard deviation `σ`.

See also [`uniformdistr`](@ref).
"""
function normaldistr(
    x::AbstractVector,
    μ::Real,
    σ::Real
)
    return @. 1 / (σ * √(2π)) * exp(-0.5((x - μ) / σ)^2)
end # normaldistr
function normaldistr(
    x::Real,
    μ::Real,
    σ::Real
)
    return 1 / (σ * √(2π)) * exp(-0.5((x - μ) / σ)^2)
end

#----------------#
# Random numbers #
#----------------#--------------------------------------------------------------
"""
    Base.randn(μ, σ, dims)
    
Method which return random variables from a normal distribution with 
expectation value `μ` and standard deviation `σ`.
"""
function Base.randn(
    rng::AbstractRNG,
    μ::AbstractFloat,
    σ::AbstractFloat,
)
    return μ + σ * randn(rng)
end
function Base.randn(
    rng::AbstractRNG,
    μ::AbstractFloat,
    σ::AbstractFloat,
    dims...
)
    return μ .+ σ .* randn(rng, typeof(μ), dims)
end
function randn!(
    r::AbstractVector,
    rng::AbstractRNG,
    μ::AbstractFloat,
    σ::AbstractFloat,
)
    for i in eachindex(r)
        r[i] = randn(rng, μ, σ)
    end
end
"""
    Base.rand(rng, a::AbstractFloat, b::AbstractFloat)
Method which return random variables from a uniform distribution on the interval
[a, b] with DataType `precision` and random seed generator `rng`.
"""
function Base.rand(
    rng::AbstractRNG,
    a::AbstractFloat,
    b::AbstractFloat,
)
    return a + rand(rng, typeof(a)) * (b - a)
end
function Base.rand(
    rng::AbstractRNG,
    a::AbstractFloat,
    b::AbstractFloat,
    dims...
)
    return a + rand(rng, typeof(a), dims...) * (b - a)
end
function Base.rand(
    rng::AbstractRNG,
    domain::Vector{<:Tuple{Real,Real}}
)
    numaxes = length(domain)
    r = zeros(numaxes)
    for i = 1:numaxes
        a = domain[i][1]
        b = domain[i][2]
        r[i] = rand(rng, a, b)
    end
    return r
end # function rand
function rand!(
    r::AbstractVector,
    rng::AbstractRNG,
    domain::Vector{<:Tuple{Real,Real}}
)
    numaxes = length(domain)
    for i in eachindex(1:numaxes)
        a = domain[i][1]
        b = domain[i][2]
        r[i] = rand(rng, a, b)
    end
    return nothing
end
function rand!(
    r::AbstractVector,
    rng::AbstractRNG,
    bounds::NTuple{N,Tuple{Real,Real}} where {N}
)
    numaxes = length(bounds)
    for i in eachindex(1:numaxes)
        a = bounds[i][1]
        b = bounds[i][2]
        r[i] = rand(rng, a, b)
    end
    return nothing
end
function rand!(
    r::AbstractVector,
    rng::AbstractRNG,
    a::AbstractFloat,
    b::AbstractFloat,
    n::Integer
)
    for i in 1:n
        r[i] = rand(rng, a, b)
    end
    return nothing
end

#----------#
# Sampling #
#----------#--------------------------------------------------------------------

"""
    maxwellianvelocitysample(
        rng::AbstractRNG,
        temperature_itp::Any,
        mass::Real,
        args...
        ;
        precision::DataType=Float64,
    )
Draw a velocity component sample from a Maxwellian distribution with temperature
given by a function call with aribtrary arguments, and by a mass.
"""
function maxwellianvelocitysample(
    rng::AbstractRNG,
    T::Real,
    mass::Real,
)
    μ = 0.0
    σ = sqrt(TraceParticles.k_B * T / mass) # Standard deviation of the Maxwell
    # distribution at this temperature
    return randn(rng, μ, σ)
end


"""
    rejectionsample(
        target  ::Any,
        maxvalue::Real,
        domain  ::Vector{<:Tuple{Real,Real}},
        rng     ::AbstractRNG
        )
Sample a point from the `target`-distribution using rejection sampling. The
`target` function must be callable with dimension (the number of arguments)
equal to the length of `domain`. The `maxvalue` determines the upper bound
of the uniform "proposal" distribution. Make sure sure that this value
is equal to or larger than the maximum possible value of the target
distribution. If not, the point will not be a true sample of the
target distribution.

The function returns the sampled point and the number of rejections.

The rejection sample algorithm draws independent and identically distributed
random variables, at the cost of a potentially large number of rejections.
The number of rejectiosn icreases with dimension.
"""
function rejectionsample(
    target::Any,
    maxvalue::Real,
    domain::Vector{<:Tuple{Real,Real}},
    rng::AbstractRNG
)
    numrejections = 0
    pos = rand(rng, domain)
    yguess = rand(rng, 0.0, maxvalue)
    ytarget = target(pos...)[1] # Target should return a single float, but in
    while yguess > ytarget
        numrejections += 1
        pos = rand(rng, domain)
        yguess = rand(rng, 0.0, maxvalue)
        ytarget = target(pos...)[1]
    end
    return pos, numrejections
end
function rejectionsample!(
    pos::AbstractVector,
    target::Any,
    maxvalue::Real,
    domain::Vector{<:Tuple{Real,Real}},
    rng::AbstractRNG
)
    numrejections = 0
    rand!(pos, rng, domain)
    yguess = rand(rng, 0.0, maxvalue)
    ytarget = target(pos...)[1] # Target should return a single float, but in
    while yguess > ytarget
        numrejections += 1
        rand!(pos, rng, domain)
        yguess = rand(rng, 0.0, maxvalue)
        ytarget = target(pos...)[1]
    end
    return pos, numrejections
end
function rejectionsample(
    rng::AbstractRNG,
    target::Any,
    maxvalue::Real,
    x0::Real,
    xf::Real,
    y0::Real,
    yf::Real,
)
    numrejections = 0
    x = rand(rng, x0, xf)
    y = rand(rng, y0, yf)
    u = rand(rng, 0.0, maxvalue)
    t = target(x, y)
    while u > t
        numrejections += 1
        x = rand(rng, x0, xf)
        y = rand(rng, y0, yf)
        u = rand(rng, 0.0, maxvalue)
        t = target(x, y)
    end
    return x, y, numrejections
end
function rejectionsample(
    rng::AbstractRNG,
    target::Any,
    maxvalue::Real,
    x0::Real,
    xf::Real,
    y0::Real,
    yf::Real,
    z0::Real,
    zf::Real,
)
    numrejections = 0
    x = rand(rng, x0, xf)
    y = rand(rng, y0, yf)
    z = rand(rng, z0, zf)
    u = rand(rng, 0.0, maxvalue)
    t = target(x, y, z)
    while u > t
        numrejections += 1
        x = rand(rng, x0, xf)
        y = rand(rng, y0, yf)
        z = rand(rng, z0, zf)
        u = rand(rng, 0.0, maxvalue)
        t = target(x, y, z)
    end
    return x, y, z, numrejections
end

function rejectionsample(
    rng::AbstractRNG,
    target::Any,
    maxvalue::Real,
    x0::Real,
    xf::Real,
    y0::Real,
    yf::Real,
    z0::Real,
    zf::Real,
    t0::Real,
    tf::Real,
)
    numrejections = 0
    x = rand(rng, x0, xf)
    y = rand(rng, y0, yf)
    z = rand(rng, z0, zf)
    t = rand(rng, t0, tf)
    u = rand(rng, 0.0, maxvalue)
    v = target(x, y, z, t)
    while u > v
        numrejections += 1
        x = rand(rng, x0, xf)
        y = rand(rng, y0, yf)
        z = rand(rng, z0, zf)
        t = rand(rng, t0, tf)
        u = rand(rng, 0.0, maxvalue)
        v = target(x, y, z, t)
    end
    return x, y, z, t, numrejections
end

"""
    rejectionsample(
        target    ::Function,
        maxvalue  ::Real,
        domain    ::Vector{<:Tuple{Real,Real}},
        rng       ::AbstractRNG,
        numsamples::Integer,
        )
Draws `numsamples` points from the `target`-distribution using the rejection
algorithm. See [`rejectionsample`](@ref) for more information.
"""
function rejectionsample(
    target,
    maxvalue::Real,
    domain::Vector{<:Tuple{Real,Real}},
    rng::AbstractRNG,
    numsamples::Integer,
)
    numdims = size(domain)[1]
    positions = zeros((numdims, numsamples))
    numaccepted = 0
    numrejected = 0
    while numaccepted < numsamples
        pos, rejections = rejectionsample(target, maxvalue, domain, rng)
        numaccepted += 1
        numrejected += rejections
        positions[:, numaccepted] .= pos
    end
    acceptancerate = numaccepted / (numaccepted + numrejected)
    if numdims == 1
        return positions[1, :], acceptancerate
    else
        return positions, acceptancerate
    end
end # function rejection sampling


"""
    binindex(values::AbstractVector, edges::AbstractVector)
Bins `values` into the bins defined by `edges`. Returns the indices of the bins.
"""
function binindex(values::AbstractVector, edges::AbstractVector)
    indices = missings(Int, length(values))
    nbins = length(edges) - 1
    Threads.@threads for i in eachindex(values)
        if values[i] == edges[end]
            indices[i] = nbins
            continue
        end
        for j in 1:nbins
            if edges[j] <= values[i] < edges[j+1]
                indices[i] = j
                break
            end
        end
    end
    return indices
end


"""
    binmap(
        data,
        ;
        weights=ones(length(data)),
        mapfunc::Function = (x,w) -> sum(w),
        xvalues=nothing,
        xedges=nothing,
        nbins= 100,
        logx=false,
        )
Maps a binning of the `data` to the function `mapfunc`. The default mapping
is the sum of the weights of the data in each bin. The defualt weighing is 1.

If no edges or data-corresponding x-values are given, we bin the data in a
linear, uniform range between its minimum and maximum value.

Returns the midpoints of the bins and the mapped values of the bins.
"""
function binmap(
    data,
    ;
    weights=ones(length(data)),
    mapfunc::Function=(x, w) -> sum(w),
    xvalues=nothing,
    xedges=nothing,
    nbins=100,
    logx=false,
)
    if isnothing(xedges) && isnothing(xvalues)
        # If no edges or data-corresponding x-values are given, we bin the
        # data in a linear, uniform range between its minimum and maximum
        # value.
        if logx
            log10data = log10.(data)
            minx = minimum(log10data)
            maxx = maximum(log10data)
            xedges = 10 .^ LinRange(minx, maxx, nbins + 1)
        else
            minx = minimum(data)
            maxx = maximum(data)
            xedges = LinRange(minx, maxx, nbins + 1)
        end
    elseif isnothing(xedges)
        # If data-corresponding x-values are given, we use these to bin the
        # data.
        if logx
            log10data = log10.(xvalues)
            minx = minimum(log10data)
            maxx = maximum(log10data)
            xedges = 10 .^ LinRange(minx, maxx, nbins + 1)
        else
            minx = minimum(xvalues)
            maxx = maximum(xvalues)
            xedges = LinRange(minx, maxx, nbins + 1)
        end
    else
        # If edges are given, we use these to bin the data.
        nbins = length(xedges) - 1
    end
    # Bin the data
    binindex_x = binindex(data, xedges)

    # Here comes maping of each bin.
    # First,
    # construct a DataFrame with the data and its corresponding binindex
    # and weight as columns.
    df = DataFrame(
        x=data,
        binindex_x=binindex_x,
        weight=weights
    )
    # Second,
    # group the data by binindex
    gdf = groupby(df, :binindex_x)
    # Third,
    # Apply the mapfunc to each exisiting group (some bins might be empty and
    # hence not present in the groupby object).
    binvalues = missings(eltype(data), nbins)
    k = keys(gdf)
    Threads.@threads for i in 1:length(gdf)
        if !ismissing(k[i].binindex_x)
            groupdf = gdf[i]
            binvalues[k[i].binindex_x] = mapfunc(
                groupdf.x, groupdf.weight
            )
        end
    end

    emptybins = sum(ismissing.(binvalues))
    if emptybins != 0
        @warn @sprintf("%.2f %% of the bins are empty.", emptybins / (nbins) * 100)
    end

    return xedges, binvalues
end


"""
    binmap(
        xvalues,
        yvalues,
        ;
        data=ones(length(xvalues)),
        mapfunc::Function = x -> length(x),
        nbinsx=100,
        nbinsy=100,
        logx=false,
        logy=false,
        )
Maps a 2D binning of the `data` to the function `mapfunc`. Does not support
weighing of the `data`.

Returns the edges of the bins and their mapped values.
"""
function binmap(
    xvalues,
    yvalues,
    ;
    data=ones(length(xvalues)),
    mapfunc::Function=(d, w) -> sum(w),
    weights=ones(length(xvalues)),
    nbinsx=100,
    nbinsy=100,
    logx=false,
    logy=false,
)
    if logx
        log10data = log10.(xvalues)
        minx = minimum(log10data)
        maxx = maximum(log10data)
        xedges = 10 .^ LinRange(minx, maxx, nbinsx + 1)
    else
        minx = minimum(xvalues)
        maxx = maximum(xvalues)
        xedges = LinRange(minx, maxx, nbinsx + 1)
    end
    if logy
        log10data = log10.(yvalues)
        miny = minimum(log10data)
        maxy = maximum(log10data)
        yedges = 10 .^ LinRange(miny, maxy, nbinsy + 1)
    else
        miny = minimum(yvalues)
        maxy = maximum(yvalues)
        yedges = LinRange(miny, maxy, nbinsy + 1)
    end

    binindex_x = binindex(xvalues, xedges)
    binindex_y = binindex(yvalues, yedges)
    df = DataFrame(
        x=xvalues,
        y=yvalues,
        d=data,
        binindex_x=binindex_x,
        binindex_y=binindex_y,
        weight=weights,
    )
    gdf = groupby(df, [:binindex_x, :binindex_y])
    #giddf = select(gdf, groupindices => :gid)

    binvalues = missings(eltype(data), nbinsx, nbinsy)
    k = keys(gdf)
    Threads.@threads for i in eachindex(k)
        if !ismissing(k[i].binindex_x) && !ismissing(k[i].binindex_y)
            groupdf = gdf[i]
            binvalues[k[i].binindex_x, k[i].binindex_y] = mapfunc(
                groupdf.d,
                groupdf.weight,
            )
        end
    end
    emptybins = sum(ismissing.(binvalues))
    if emptybins != 0
        @warn @sprintf("%.2f %% of the bins are empty.", emptybins / (nbinsx * nbinsy) * 100)
    end
    return binvalues, (xedges, yedges)
end


"""
    binmap(
        xvalues,
        yvalues,
        zvalues,
        ;
        data=ones(length(xvalues)),
        mapfunc::Function = x -> length(x),
        nbins=100
        logx=false,
        logy=false,
        )
Maps a 2D binning of the `data` to the function `mapfunc`. Does not support
weighing of the `data`.

Returns the edges of the bins and their mapped values.
"""
function binmap(
    xvalues,
    yvalues,
    zvalues
    ;
    data=ones(length(xvalues)),
    mapfunc::Function=(d, w) -> sum(w),
    weights=ones(length(xvalues)),
    nbinsx=100,
    nbinsy=100,
    nbinsz=100,
    logx=false,
    logy=false,
    logz=false,
)
    if logx
        log10data = log10.(xvalues)
        minx = minimum(log10data)
        maxx = maximum(log10data)
        xedges = 10 .^ LinRange(minx, maxx, nbinsx + 1)
    else
        minx = minimum(xvalues)
        maxx = maximum(xvalues)
        xedges = LinRange(minx, maxx, nbinsx + 1)
    end
    if logy
        log10data = log10.(yvalues)
        miny = minimum(log10data)
        maxy = maximum(log10data)
        yedges = 10 .^ LinRange(miny, maxy, nbinsy + 1)
    else
        miny = minimum(yvalues)
        maxy = maximum(yvalues)
        yedges = LinRange(miny, maxy, nbinsy + 1)
    end
    if logz
        log10data = log10.(zvalues)
        minz = minimum(log10data)
        maxz = maximum(log10data)
        zedges = 10 .^ LinRange(minz, maxz, nbinsz + 1)
    else
        minz = minimum(zvalues)
        maxz = maximum(zvalues)
        zedges = LinRange(minz, maxz, nbinsz + 1)
    end

    binindex_x = binindex(xvalues, xedges)
    binindex_y = binindex(yvalues, yedges)
    binindex_z = binindex(zvalues, zedges)
    df = DataFrame(
        x=xvalues,
        y=yvalues,
        z=zvalues,
        d=data,
        binindex_x=binindex_x,
        binindex_y=binindex_y,
        binindex_z=binindex_z,
        weight=weights,
    )
    gdf = groupby(df, [:binindex_x, :binindex_y, :binindex_z])
    #giddf = select(gdf, groupindices => :gid)

    binvalues = missings(eltype(eltype(data)), nbinsx, nbinsy, nbinsz)
    k = keys(gdf)
    Threads.@threads for i in eachindex(k)
        if !ismissing(k[i].binindex_x) &&
           !ismissing(k[i].binindex_y) &&
           !ismissing(k[i].binindex_z)
            groupdf = gdf[i]
            binvalues[k[i].binindex_x, k[i].binindex_y, k[i].binindex_z] = mapfunc(
                groupdf.d,
                groupdf.weight,
            )
        end
    end
    emptybins = sum(ismissing.(binvalues))
    if emptybins != 0
        @info @sprintf("%.2f %% of the bins are empty.", emptybins / (nbinsx * nbinsy * nbinsz) * 100)
    end
    return binvalues, (xedges, yedges, zedges)
end

"""
    normfactor_2Duniformmesh(var, axes; xlim=nothing, zlim=nothing)
Calculate the normalisation factor of the 2D, gridded scalar-field `var` in the
region defined by `xlim` and `zlim`. If `xlim` and `zlim` are not given, the
normalisation factor is calculated over the entire grid.
"""
function normfactor_2Duniformmesh(
    var::Matrix{<:Real},
    axes::Tuple{AbstractVector,AbstractVector};
    xlim=nothing,
    ylim=nothing,
)
    if !(isnothing(xlim) && isnothing(ylim))
        slicex, slicey = findslice(
            axes[1],
            axes[2],
            xlim,
            ylim,
        )
        var = var[slicex, slicey]
        axes = (axes[1][slicex], axes[2][slicey])
    end
    dxvec = [diff(ax)[1] for ax in axes]
    dv = prod(dxvec)
    return sum(var) * dv
end

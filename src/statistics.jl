# Created 22.01.24
# Author: e.s.oyre@astro.uio.no
#
#           statistics.jl
#-------------------------------------------------------------------------------
# Contains random number generators, samplers, and functions for binning
# particle data.
#-------------------------------------------------------------------------------

#----------------#
# Random numbers #
#----------------#--------------------------------------------------------------
"""
    Base.randn(μ, precision, σ, dims)
Method which return random variables from a normal distribution with 
expectation value `μ` and standard deviation `σ`.
"""
function Base.randn(
    rng::AbstractRNG,
    precision::DataType,
    μ::AbstractFloat,
    σ::AbstractFloat,
)
    return μ + σ * randn(rng, precision)
end

"""
    Base.rand(rng, precision, a::AbstractFloat, b::AbstractFloat)
Method which return random variables from a uniform distribution on the interval
[a, b] with DataType `precision` and random seed generator `rng`.
"""
function Base.rand(
    rng::AbstractRNG,
    precision::DataType,
    a::AbstractFloat,
    b::AbstractFloat,
)
    return a + rand(rng, precision) * (b - a)
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
    precision::DataType,
    T::Real,
    mass::Real,
)
    μ = precision(0)
    σ = sqrt(TraceParticles.k_B * T / mass) # Standard deviation of the Maxwell
    # distribution at this temperature
    return randn(rng, precision, μ, σ)
end


"""
    rejectionsample(
        rng::AbstractRNG,
        precision::DataType,
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
Sample a point from the 4D `target`-distribution using rejection sampling. The
`target` function must be callable with the signature (x, y, z, t). The
`maxvalue` determines the upper bound of the uniform "proposal" distribution.
Make sure sure that this value is equal to or larger than the maximum possible
value of the target distribution. If not, the point will not be a true sample
of the target distribution.

The function returns the sampled point and the number of rejections.

The rejection sample algorithm draws independent and identically distributed
random variables, at the cost of a potentially large number of rejections.
The number of rejections icreases with the number of dimensions.
"""
function rejectionsample(
    rng::AbstractRNG,
    precision::DataType,
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
    x = rand(rng, precision, x0, xf)
    y = rand(rng, precision, y0, yf)
    z = rand(rng, precision, z0, zf)
    t = rand(rng, precision, t0, tf)
    u = rand(rng, precision, 0.0, maxvalue)
    v = target(x, y, z, t)
    while u > v
        numrejections += 1
        x = rand(rng, precision, x0, xf)
        y = rand(rng, precision, y0, yf)
        z = rand(rng, precision, z0, zf)
        t = rand(rng, precision, t0, tf)
        u = rand(rng, precision, 0.0, maxvalue)
        v = target(x, y, z, t)
    end
    return x, y, z, t, numrejections
end

"""
    mhdsample(
        mass::Real,
        rng::AbstractRNG,
        proposal_distr,
        target_distr,
        tg_itp,
        xbounds::Tuple{Real,Real},
        ybounds::Tuple{Real,Real},
        zbounds::Tuple{Real,Real},
        temporalbounds::Tuple{Real,Real},
        max_value::Real,
    )
Draw `x`, `y`, `z` positions and time `t` from a `proposal_distr`ibution using
rejection sampling. Draws the 3D velocity vector from a Maxwellian distribution
with a temperature given by `tg_itp(x,y,z)`. Evaluates the statistical
importance weight of the particle using the `target_distr`ibution.
"""
function mhdsample(
    mass::Real,
    rng::AbstractRNG,
    proposal_distr,
    target_distr,
    tg_itp,
    xbounds::T,
    ybounds::T,
    zbounds::T,
    temporalbounds::T,
    max_value::Real;
    precision=typeof(mass),
) where {T<:Tuple{Real,Real}}
    # Extract spatial and temporal limits
    x0, xf = xbounds
    y0, yf = ybounds
    z0, zf = zbounds
    tmin, tmax = temporalbounds
    # Sample osition and time from the proposal distrubution using rejection
    # sampling.
    x, y, z, t, nrejections = rejectionsample(
        rng,
        precision,
        proposal_distr,
        max_value,
        x0, xf,
        y0, yf,
        z0, zf,
        tmin, tmax,
    )
    weight = precision(
        target_distr(x, y, z, t) / proposal_distr(x, y, z, t)
    )
    # Sample velocity components from a Maxwell-Boltzmann distribution.
    temperature = tg_itp(x, y, z, t)
    vx = maxwellianvelocitysample(rng, precision, temperature, mass)
    vy = maxwellianvelocitysample(rng, precision, temperature, mass)
    vz = maxwellianvelocitysample(rng, precision, temperature, mass)
    return x, y, z, vx, vy, vz, t, weight, nrejections
end

"""
    MHDSampler(
        target_distr
        proposal_distr
        tg_itp
        xbounds
        ybounds
        zbounds
        t0bounds
        max_value
        precision
    )
A convenience functor for doing repeated [`mhdsample`](@ref)s from the same MHD
environment.

# Methods
    (::MHDSampler)(mass, rng)
Returns `x`, `y`, `z`, `vx`, `vy`, `vz`, `t`, `weight` and `nrejections`.
"""
struct MHDSampler{T1,T2,T3,T4<:AbstractFloat,T5<:Real,T6<:DataType}
    target_distr::T1
    proposal_distr::T2
    tg_itp::T3
    xbounds::Tuple{T4,T4}
    ybounds::Tuple{T4,T4}
    zbounds::Tuple{T4,T4}
    t0bounds::Tuple{T4,T4}
    max_value::T5
    precision::T6
end
function (self::MHDSampler)(mass, rng)
    return mhdsample(
        mass,
        rng,
        self.proposal_distr,
        self.target_distr,
        self.tg_itp,
        self.xbounds,
        self.ybounds,
        self.zbounds,
        self.t0bounds,
        self.max_value;
        precision=self.precision
    )
end


#---------#
# Binning #
#---------#---------------------------------------------------------------------
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

    binvalues = missings(Any, nbinsx, nbinsy)
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
        xedges = range(minx, maxx, length=nbinsx + 1)
    end
    if logy
        log10data = log10.(yvalues)
        miny = minimum(log10data)
        maxy = maximum(log10data)
        yedges = 10 .^ LinRange(miny, maxy, nbinsy + 1)
    else
        miny = minimum(yvalues)
        maxy = maximum(yvalues)
        yedges = range(miny, maxy, length=nbinsy + 1)
    end
    if logz
        log10data = log10.(zvalues)
        minz = minimum(log10data)
        maxz = maximum(log10data)
        zedges = 10 .^ LinRange(minz, maxz, nbinsz + 1)
    else
        minz = minimum(zvalues)
        maxz = maximum(zvalues)
        zedges = range(minz, maxz, length=nbinsz + 1)
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

    binvalues = missings(Any, nbinsx, nbinsy, nbinsz)
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
    hist(
        data,
        weight;
        nbins=20,
        logx=false,
        maxval=nothing,
        minval=nothing,
        histmode=:pdf
    )
Calculate a histogram of `data` with `nbins` bins, weighted by `weight`,
and normalised according to `histmode`. If `logx` is true, the bins are
logarithmically spaced. If `maxval` or `minval` are given, these values
are used as the maximum and minimum values of the histogram bins.
"""
function hist(
    data,
    weight;
    nbins=20,
    logx=false,
    maxval=nothing,
    minval=nothing,
    histmode=:pdf
)
    loweredge = isnothing(minval) ? minimum(data) : minval
    upperedge = isnothing(maxval) ? maximum(data) : maxval
    if logx
        binedges = 10.0 .^ range(
            log10(loweredge),
            log10(upperedge),
            length=nbins + 1
        )
    else
        binedges = range(loweredge, upperedge, length=nbins + 1)
    end
    histcounts = StatsBase.fit(Histogram, data, weights(weight), binedges)
    histcounts = StatsBase.normalize(histcounts, mode=histmode)
    x = histcounts.edges[1]
    y = histcounts.weights
    return x, y
end

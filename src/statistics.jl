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

See also [`Utilities.uniformdistr`](@ref).
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

function bivariate_normaldistr(
    x::Real,
    y::Real,
    μx::Real,
    μy::Real,
    σx::Real,
    σy::Real,
    corr::Real # correlation
)
    normfactor = 1 / (2π * σx * σy * √(1 - corr^2))
    return normfactor * exp(
        -1 / (2 * (1 - corr^2)) * (
            (x - μx)^2 / σx^2 +
            (y - μy)^2 / σy^2 -
            2 * corr * (x - μx) * (y - μy) / (σx * σy)
        )
    )
end

"""
    uniformdistr(
        x::Vector{T} where {T<:Real},
        a::Real,
        b::Real
        )
Function returning the values of `x` on a 1D normalised unifrom distribution on
the interval [`a, `b`].

See also [`Utilities.normaldistr`](@ref).
"""
function uniformdistr(
    x::Array{T} where {T<:Real},
    a::Real,
    b::Real
)
    mask = a .<= x .<= b
    stepheight = 1.0 / (b - a)
    prob = zeros(typeof(x[1]), size(x))
    prob[mask] .= stepheight
    return prob
end # function uniformdistr


"""
    bivarate_uniformdistr(x,y,a,b,c,d)
Return the probability of (`x`,`y`) in the bivariate uniform distribution
f(x,y) = U(a,b) x U(c,d)
"""
function bivariate_uniformdistr(
    x::Real,
    y::Real
    ;
    a::Real=0.0,
    b::Real=1.0,
    c::Real=0.0,
    d::Real=1.0
)
    Δx = b - a
    Δy = d - c
    probability = abs(1 / (Δx * Δy))
    return a < x <= b && c < y <= d ? probability : 0.0
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
    μ::AbstractFloat,
    σ::AbstractFloat,
    dims...
)
    return μ .+ σ .* randn(dims)
end
function Base.randn(
    precision::DataType,
    μ::AbstractFloat,
    σ::AbstractFloat,
    dims...
)
    μ = convert.(precision, μ)
    σ = convert.(precision, σ)
    return μ .+ σ .* randn(precision, dims)
end
function Base.randn(
    rng::AbstractRNG,
    precision::DataType,
    μ::AbstractFloat,
    σ::AbstractFloat,
    dims...
)
    μ = convert.(precision, μ)
    σ = convert.(precision, σ)
    return μ .+ σ .* randn(rng, precision, dims)
end


"""
    Base.rand(rng, precision, a, b, dims)
Method which return random variables from a uniform distribution on the interval
[a, b] with DataType `precision` and random seed generator `rng`.
"""
function Base.rand(
    a::AbstractFloat,
    b::AbstractFloat,
    dims...
)
    return a .+ rand(dims...) .* (b - a)
end
function Base.rand(
    precision::DataType,
    a::AbstractFloat,
    b::AbstractFloat,
    dims...
)
    a = convert(precision, a)
    b = convert(precision, b)
    return a .+ rand(precision, dims...) .* (b - a)
end
function Base.rand(
    rng::AbstractRNG,
    a::AbstractFloat,
    b::AbstractFloat,
    dims...
)
    return a .+ rand(rng, dims...) .* (b - a)
end
function Base.rand(
    rng::AbstractRNG,
    precision::DataType,
    a::AbstractFloat,
    b::AbstractFloat,
    dims...
)
    a = convert(precision, a)
    b = convert(precision, b)
    return a .+ rand(rng, precision, dims...) .* (b - a)
end
function Base.rand(
    a::AbstractFloat,
    b::AbstractFloat,
    dims::Tuple{Vararg{Int}}
)
    return a .+ rand(dims...) .* (b - a)
end
function Base.rand(
    precision::DataType,
    a::AbstractFloat,
    b::AbstractFloat,
    dims::Tuple{Vararg{Int}}
)
    a = convert(precision, a)
    b = convert(precision, b)
    return a .+ rand(precision, dims...) .* (b - a)
end
function Base.rand(
    rng::AbstractRNG,
    a::AbstractFloat,
    b::AbstractFloat,
    dims::Tuple{Vararg{Int}}
)
    return a .+ rand(rng, dims...) .* (b - a)
end
function Base.rand(
    rng::AbstractRNG,
    precision::DataType,
    a::AbstractFloat,
    b::AbstractFloat,
    dims::Tuple{Vararg{Int}}
)
    a = convert(precision, a)
    b = convert(precision, b)
    return a .+ rand(rng, precision, dims...) .* (b - a)
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
        r[i] = a .+ rand(rng) .* (b - a)
    end
    return r
end # function rand


#----------#
# Sampling #
#----------#--------------------------------------------------------------------
function maxwellianvelocitysample(
    rng::AbstractRNG,
    temperature_itp::Any,
    mass::Real,
    args...
    ;
    precision::DataType=Float64,
)
    T = temperature_itp.(args...)
    μ = 0.0
    σ = sqrt.(tp.k_B * T / mass) # Standard deviation of the Maxwell
    # distribution at this temperature
    return randn.(rng, precision, μ, σ, 3)
end


"""
    importancesampling(
        target  ::Function, 
        proposal::Function, 
        randgen ::Function, 
        N       ::Integer,    
        )
Sample `N` points from the `proposal`-distribution and compute the importance
weights with respect to the `target`-distribution.
"""
function importancesampling(
    target::Function, # Target distribution
    proposal::Function, # Proposal distruv
    randgen::Function, # Random variable generator. Following proposal pdf.
    dims::Tuple{Vararg{Integer}}, # Number of samples
)
    samples = randgen(dims)
    weights = target(samples) ./ proposal(samples)
    return samples, weights
end # function importancesampling


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
        args...
        ;
        weights=ones(length(data)),
        mapfunc::Function = (x,w) -> sum(w),
        xvalues=nothing,
        xedges=nothing,
        nbins= 100,
        logx=false,
        logf=false,
        normalise=false,
        )
Maps a binning of the `data` to the function `mapfunc`. The default mapping
is the sum of the weights of the data in each bin. The defualt weighing is 1.

If no edges or data-corresponding x-values are given, we bin the data in a
linear, uniform range between its minimum and maximum value.

Returns the midpoints of the bins and the mapped values of the bins.
"""
function binmap(
    data,
    args...
    ;
    weights=ones(length(data)),
    mapfunc::Function=(x, w) -> sum(w),
    xvalues=nothing,
    xedges=nothing,
    nbins=100,
    logx=false,
    logf=false,
    normalise=false,
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
                groupdf.x, groupdf.weight, args...
            )
        end
    end

    emptybins = sum(ismissing.(binvalues))
    if emptybins != 0
        @warn @sprintf("%.2f %% of the bins are empty.", emptybins / (nbins) * 100)
    end

    # Some post-processing. Normalise the binvalues if requested.
    if normalise
        if logx
            error("Normalisation for non-uniform binning has a bug.")
        end
        mask = ismissing.(binvalues)
        bm = binvalues[.!mask]
        dx = diff(xedges)[.!mask]
        normfactor = sum(dx .* bm)
        binvalues = binvalues ./ normfactor
    end
    if logf
        binvalues = log10.(binvalues)
    end
    return xedges, binvalues
end


"""
    binmap(
        xvalues,
        yvalues,
        zvalues,
        mapfunc::Function = x -> length(x),
        args...
        ;
        nbinsx=100,
        nbinsy=100,
        logaxes=false,
        logx=false,
        logy=false,
        logz=false,
        normalise=false,
        )
Maps a 2D binning of the `zvalues` to the function `mapfunc`. Does not support
weighing of the `zvalues`.

Returns the 2D midpoints of the bins and their mapped values.
"""
function binmap(
    xvalues,
    yvalues,
    zvalues,
    ;
    mapfunc::Function=(z, w) -> sum(w),
    weights=ones(length(xvalues)),
    nbinsx=100,
    nbinsy=100,
    logx=false,
    logy=false,
    logz=false,
    normalise=false,
)
    #if logaxes
    #    xvalues = log10.(xvalues)
    #    yvalues = log10.(yvalues)
    #elseif logx
    #    xvalues = log10.(xvalues)
    #elseif logy
    #    yvalues = log10.(yvalues)
    #end
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

    #xedges = LinRange(minimum(xvalues), maximum(xvalues), nbinsx+1)
    #yedges = LinRange(minimum(yvalues), maximum(yvalues), nbinsy+1)
    binindex_x = binindex(xvalues, xedges)
    binindex_y = binindex(yvalues, yedges)
    df = DataFrame(
        x=xvalues,
        y=yvalues,
        z=zvalues,
        binindex_x=binindex_x,
        binindex_y=binindex_y,
        weight=weights,
    )
    gdf = groupby(df, [:binindex_x, :binindex_y])
    #giddf = select(gdf, groupindices => :gid)

    binvalues = missings(eltype(zvalues), nbinsx, nbinsy)
    k = keys(gdf)
    Threads.@threads for i in eachindex(k)
        if !ismissing(k[i].binindex_x) && !ismissing(k[i].binindex_y)
            groupdf = gdf[i]
            binvalues[k[i].binindex_x, k[i].binindex_y] = mapfunc(
                groupdf.z,
                groupdf.weight,
            )
        end
    end
    emptybins = sum(ismissing.(binvalues))
    if emptybins != 0
        @warn @sprintf("%.2f %% of the bins are empty.", emptybins / (nbinsx * nbinsy) * 100)
    end
    if normalise
        if logx || logy
            error("Normalisation for non-uniform binning is not implemented.")
        end
        dx = diff(xedges)[1]
        dy = diff(yedges)[1]
        dA = dx * dy
        binvalues = binvalues ./ (dA * sum(binvalues))
    end
    if logz
        binvalues = log10.(binvalues)
    end
    return midpoints(xedges), midpoints(yedges), binvalues
end


function binmap_witherror(
    data,
    nsplits,
    args...
    ;
    kwargs...
)

    xx, mean_bm = binmap(data, args...; kwargs...)
    bm_arr = Matrix{Float64}(undef, nsplits, length(mean_bm))
    ni = length(data)
    i = 1
    stride = floor(Int, ni / nsplits)
    for k in 1:nsplits
        if k != nsplits
            _, bm = binmap(data[i:i+stride-1], args...; kwargs...)
        else
            _, bm = binmap(data[i:end], args...; kwargs...)
        end
        i += stride
        bm_arr[k, :] .= bm
    end
    std_bm = mapslices(std, bm_arr; dims=1)
    error_bm = std_bm ./ sqrt(nsplits)
    return xx, mean_bm, error_bm
end


weighted_average(value, weight) = sum(value .* weight) / sum(weight)
weighted_count(value, weight) = sum(weight)

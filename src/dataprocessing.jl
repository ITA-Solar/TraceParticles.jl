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
    find_nonthermals(
        fname::String,
        nparticles::Int;
        relativegain::Bool=false,
        )
Find the `nparticles` most energetic particles from the test particle ensemble
stored as HDF5 files with the filename `fname`.
"""
function find_nonthermals(
    fname::String,
    nparticles::Int;
    relativegain::Bool=false,
)
    expname, _ = splitext(fname)
    x0, y0, z0, vparal0, mu0 = h5_getinitialstate(fname)
    e0, ef = h5_getenergies(expname)

    sortby = relativegain ? (ef .- e0) ./ e0 : ef
    idxs = sortperm(sortby)[end-nparticles+1:end]
    u0 = [[x0[i], y0[i], z0[i], vparal0[i]] for i in idxs]
    mu0 = mu0[idxs]

    return u0, mu0, idxs
end


"""
    find_maxiters(
        fname::String,
        nparticles::Int;
    )
Find the `nparticles` particles with most timesteps from the test particle
ensemble stored as HDF5 files with the filename `fname`.
"""
function find_maxiters(
    fname::String,
    nparticles::Int;
)
    x0, y0, z0, vparal0, mu0 = h5_getinitialstate(fname)
    nt = h5_getdataset(fname, "nt")
    sortby = nt
    idxs = sortperm(sortby)[end-nparticles+1:end]
    u0 = [[x0[i], y0[i], z0[i], vparal0[i]] for i in idxs]
    mu0 = mu0[idxs]

    return u0, mu0, idxs
end

function find_relativistic(
    fname::String,
)
    x0, y0, z0, vparal0, mu0 = h5_getinitialstate(fname)
    retmsg = h5_getdataset(fname, "retmsg")
    idxs = findall(x -> x == "Relativistic", retmsg)
    u0 = [[x0[i], y0[i], z0[i], vparal0[i]] for i in idxs]
    mu0 = mu0[idxs]

    return u0, mu0, idxs
end


function find_idxs(
    fname::String,
    idxs::Vector{Int}
)
    x0, y0, z0, vparal0, mu0 = h5_getinitialstate(fname)
    u0 = [[x0[i], y0[i], z0[i], vparal0[i]] for i in idxs]
    mu0 = mu0[idxs]

    return u0, mu0
end


"""
    get_u0(
        fname::String;
        mask=nothing
    )
Get the initial state of the test particles stored in the HDF5-file `fname`.
Optionally filter out particles using a mask.
"""
function get_u0(
    fname::String;
    mask=nothing
)
    expname, _ = splitext(fname)
    x0, y0, z0, vparal0, mu0 = h5_getinitialstate(expname)
    u0 = [[x0[i], y0[i], z0[i], vparal0[i]] for i in eachindex(x0)]

    if isnothing(mask)
        return u0, mu0
    else
        return u0[mask], mu0[mask]
    end
end


"""
    generate_probdistr(
        xvalues,
        yvalues,
        zvalues,
        ;
        kwargs...
        )
Generate a 2D probability distribution from discrete values of f_z(x,y).

The algorithm bins the `zvalues` on the xy-plane and uses the mean z-value
of each bin as the value of the distribution function at that point.
Finally it normalizes the distribution function to 1 and constructs a linear
interpolation object from the discrete probability distribution.

Returns the interpolation object, and the bin edges in x and y.
"""
function generate_probdistr(
    xvalues,
    yvalues,
    zvalues,
    args...
    ;
    kwargs...
)
    xx, yy, bm = binmap(xvalues, yvalues, zvalues, args...; kwargs...)
    normbm = copy(bm)
    # Empty bins have NaN-values. Set them to zero to reflect their values in
    # the probability distribution. Interpolations.jl doesn't like NaN values
    # either.
    normbm[ismissing.(normbm)] .= 0
    if minimum(normbm) < 0
        normbm .+= abs(minimum(normbm))
    end
    dxs = [diff(xx)[1], diff(yy)[1]]
    proddxs = abs(prod(dxs))
    normfactor = sum(normbm) * proddxs
    normbm /= normfactor
    itp = linear_interpolation((xx, yy), normbm, extrapolation_bc=Flat())
    return itp, xx, yy, normbm, bm
end


"""
    saturate_distribution(
        xx,
        yy,
        normbm,
        threshold
        )
Saturate the distribution over a threshold value.
"""
function saturate_distribution(
    xx,
    yy,
    gridded_dist,
    threshold;
    from_below=true
)
    saturated_dist = copy(gridded_dist)
    if from_below
        saturated_dist[saturated_dist.<threshold] .= threshold
    else
        saturated_dist[saturated_dist.>threshold] .= threshold
    end
    dxs = [diff(xx)[1], diff(yy)[1]]
    proddxs = abs(prod(dxs))
    normfactor = sum(saturated_dist) * proddxs
    saturated_dist /= normfactor
    itp = linear_interpolation((xx, yy), saturated_dist, extrapolation_bc=Flat())
    return itp, saturated_dist
end


"""
    energyfraction(
        energies::Vector{<:Real},
        fractionlimits::Tuple{Any,Any};
        weights=ones(length(energies)),
        logx=true,
        nbins=200,
    )
Calculate the fraction of the total energy in the specified range of energies.
"""
function energyfraction(
    energies::Vector{<:Real},
    fractionlimits::Tuple{Any,Any};
    weights=ones(length(energies)),
    logx=true,
    nbins=200,
)
    # Create the bin edges
    if logx
        binedges = 10.0 .^
                   range(
            log10(minimum(energies)),
            log10(maximum(energies)),
            length=nbins + 1
        )
    else
        binedges = range(
            minimum(energies),
            maximum(energies),
            length=nbins + 1
        )
    end

    # Fit the histogram
    hist = fit(
        Histogram,
        energies,
        Weights(weights),
        binedges
    )
    hist = normalize(hist, mode=:pdf)
    bincentres = midpoints(binedges)
    binwidths = diff(binedges)

    # Finda all indices of the bins that are within the limits
    lowerlimit = fractionlimits[1]
    upperlimit = fractionlimits[2]
    if isnothing(lowerlimit) && isnothing(upperlimit)
        fraction_idxs = findall(
            e -> true,
            bincentres
        )
    elseif isnothing(lowerlimit)
        fraction_idxs = findall(
            e -> e <= upperlimit,
            bincentres
        )
    elseif isnothing(upperlimit)
        fraction_idxs = findall(
            e -> lowerlimit < e,
            bincentres
        )
    else
        fraction_idxs = findall(
            e -> lowerlimit < e <= upperlimit,
            bincentres
        )
    end

    # Calculate the fraction of the total energy in the specified range
    nonthermalcounts = hist.weights[fraction_idxs]
    nothermalbinwidhts = binwidths[fraction_idxs]
    return sum(nonthermalcounts .* nothermalbinwidhts)
end


#____/\_____/\_________________________________________________________________
#
# Filtering/masking data
#


"""
    particlemaskfunction(
        ;
        x0lim=nothing,
        y0lim=nothing
        z0lim=nothing,
        xflim=nothing,
        yflim=nothing
        zflim=nothing,
        timelim=nothing,
        retmsg=nothing,
)
Create a function that filters out particles based on the specified limits.
The function accepts a dictionary which it filters based on the key-limit pairs
- `:x0` for `x0lim`
- `:y0` for `y0lim`
- `:z0` for `z0lim`
- `:xf` for `xflim`
- `:yf` for `yflim`
- `:zf` for `zflim`
- `:tf` for `timelim`
- `:retmsg` for `retmsg`

It assumes the dictionary has the key `:x0` to find out how many particles there
are in the dataset.
"""
function particlemaskfunction(
    ;
    x0lim=nothing,
    y0lim=nothing,
    z0lim=nothing,
    xflim=nothing,
    yflim=nothing,
    zflim=nothing,
    timelim=nothing,
    retmsg=nothing,
)

    conditions = Dict{Symbol,Function}()
    if !isnothing(retmsg)
        conditions[:retmsg] = x -> retmsg == x
    end
    for (sym, limit) in zip(
        (:x0, :y0, :z0, :xf, :yf, :zf, :tf),
        (x0lim, y0lim, z0lim, xflim, yflim, zflim, timelim)
    )
        if !isnothing(limit)
            conditions[sym] = x -> limit[1] <= x <= limit[2]
        end
    end
    return (data) -> begin
        findall(
            i -> begin
                all([func(data[sym][i]) for (sym, func) in conditions])
            end,
            eachindex(data[:x0])
        )
    end
end

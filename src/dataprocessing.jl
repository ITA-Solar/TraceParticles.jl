#-------------------------------------------------------------------------------
# Created 29.01.24
# Author: e.s.oyre@astro.uio.no
#
#                 dataprocessing.jl
#
# Contains functions for processing results.
#-------------------------------------------------------------------------------


#_______________________________________________________________________________
#
# Finding particular particles
#

"""
    initialstate_nonthermals(
        fname::String,
        nparticles::Int;
        relativegain::Bool=false,
        )
Find the initial state of the `nparticles` most energetic particles from the
test particle ensemble stored as HDF5 files with the filename `fname`.
"""
function initialstate_nonthermals(
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
    initialstate_maxiters(fname::String, nparticles::Int)
Find the initial state of the `nparticles` particles with most timesteps from
the test particle ensemble stored as HDF5 files with the filename `fname`.
"""
function initialstate_maxiters(fname::String, nparticles::Int)
    x0, y0, z0, vparal0, mu0 = h5_getinitialstate(fname)
    nt = h5_getdataset(fname, "nt")
    sortby = nt
    idxs = sortperm(sortby)[end-nparticles+1:end]
    u0 = [[x0[i], y0[i], z0[i], vparal0[i]] for i in idxs]
    mu0 = mu0[idxs]
    return u0, mu0, idxs
end


"""
    initialstate_terminationcode(fname, terminationcode)
Find the initial state of the particles with a certain `terminationcode` from
the test particle ensemble stored as HDF5 with the filename `fname`.
"""
function initialstate_terminationcode(fname::String, terminationcode::String)
    x0, y0, z0, vparal0, mu = h5_getinitialstate(fname)
    idxs = findall(x -> x == terminationcode, h5_getdataset(fname, "terminationcode"))
    u0 = [[x0[i], y0[i], z0[i], vparal0[i]] for i in idxs]
    mu = mu[idxs]
    return u0, mu, idxs
end


"""
    finalstate_terminationcode(fname, terminationcode)
Find the final state of the particles with a certain `terminationcode` from
the test particle ensemble stored as HDF5 with the filename `fname`.
"""
function finalstate_terminationcode(fname::String, terminationcode::String)
    xf, yf, zf, vparalf, mu = h5_getfinalstate(fname)
    idxs = findall(x -> x == terminationcode, h5_getdataset(fname, "terminationcode"))
    uf = [[xf[i], yf[i], zf[i], vparalf[i]] for i in idxs]
    mu = mu[idxs]
    return uf, mu, idxs
end


"""
    initialstate_idxs(fname, idxs)
Find the initial state of the particles with the indices `idxs` from the test
particle ensemble stored as HDF5 with the filename `fname`.
"""
function initialstate_idxs(fname::String, idxs::AbstractVector; mask=nothing)
    x0, y0, z0, vparal0, mu0 = h5_getinitialstate(fname)
    if !isnothing(mask)
        x0 = x0[mask]
        y0 = y0[mask]
        z0 = z0[mask]
        vparal0 = vparal0[mask]
        mu0 = mu0[mask]
    end
    u0 = [[x0[i], y0[i], z0[i], vparal0[i]] for i in idxs]
    mu0 = mu0[idxs]
    return u0, mu0
end


"""
    replace_finalstates(fname_original, fname_replacement, original_idxs)
    replace_finalstates(fname_original, fname_replacement, original_idxs;
        replace_keys=names)
    replace_finalstates(original_df, replacement_df, original_idxs)
Replace the particles in `fname_original` with the particles in
`fname_replacement`. The indices of the particles to be replaced in
`fname_original` are given by `idxs`. The length of `idxs` must be the same as
the number of particles in `fname_replacement`.
"""
function replaceparticles(
    fname_original::String,
    fname_replacement::String;
    original_idxs::String="original_idxs",
    kwargs...
)
    idxs = h5open(fname_replacement) do fid
        read(fid["metadata/$original_idxs"])
    end
    replaceparticles(fname_original, fname_replacement, idxs; kwargs...)
end
function replaceparticles(
    fname_original::String,
    fname_replacement::String,
    original_idxs::Vector{<:Int};
    replace_keys::Function=names
)
    original = h5_getall(fname_original)
    replacement = h5_getall(fname_replacement)

    replaceparticles(
        original,
        replacement,
        original_idxs;
        replace_keys=replace_keys
    )
end
function replaceparticles(
    original::DataFrame,
    replacement::DataFrame,
    original_idxs::Vector{<:Int};
    replace_keys::Function=names
)
    if length(original_idxs) != size(replacement, 1)
        error(
            "Length of `original_idxs` must be equal to the number of
            particles in the replacement"
        )
    end

    out = copy(original)
    for i in eachindex(original_idxs)
        for key in replace_keys(replacement)
            out[original_idxs[i], key] = replacement[i, key]
        end
    end
    return out
end

#_______________________________________________________________________________
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
        terminationcode=nothing,
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
- `:terminationcode` for `terminationcode`

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
    terminationcode=nothing,
)

    conditions = Dict{Symbol,Function}()
    if !isnothing(terminationcode)
        conditions[:terminationcode] = x -> terminationcode == x
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


"""
    particlemask(df, outerbounds, ignorezones, patternmatches, premask)
Create a mask for the DataFrame `df` based on the specified conditions.
- `outerbounds` is a dictionary with the keys being the column names and the
  values being the limits of the outer bounds of that column.
- `ignorezones` is a vector of dictionaries with the keys being the column names
  and the values being the limits of the "zone" to ignore particles from.
- `patternmatches` is a dictionary with the keys being the column names and the
  values being a pattern of the column values to match.
"""
function particlemask(
    df::DataFrame;
    outerbounds::Dict{Symbol,<:Any}=Dict{Symbol,Any}(),
    ignorezones::Vector{<:Dict{Symbol,<:Any}}=[Dict{Symbol,Any}()],
    patternmatches::Dict{Symbol,String}=Dict{Symbol,String}(),
    premask=nothing
)
    numparticles = size(df, 1)
    if isnothing(premask)
        mask = fill(true, numparticles)
    else
        mask = premask
    end

    # Check whether the particle is within the outer bounds
    for (sym, limit) in outerbounds
        mask = mask .& (limit[1] .<= df[!, sym] .<= limit[2])
    end

    # Check whether the particle is within certain zone
    for zone in ignorezones
        if isempty(zone)
            continue
        end
        insidezone = fill(true, numparticles)
        for (sym, limit) in zone
            insidezone = insidezone .& (limit[1] .<= df[!, sym] .<= limit[2])
        end
        # If the particle is inside the inner bounds, it should not be included
        # in the mask.
        mask = mask .& .!insidezone
    end


    # Check whether the particle matches any patterns.
    for (sym, pattern) in patternmatches
        mask = mask .& (df[!, sym] .== pattern)
    end
    return mask
end


#_______________________________________________________________________________
#
# Other
#

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
    args...
    ;
    kwargs...
)
    bm, edges = binmap(args...; kwargs...)
    edges = midpoints.(edges)
    normbm = copy(bm)
    # Empty bins have NaN-values. Set them to zero to reflect their values in
    # the probability distribution. Interpolations.jl doesn't like NaN values
    # either.
    normbm[ismissing.(normbm)] .= 0
    if minimum(normbm) < 0
        normbm .+= abs(minimum(normbm))
    end
    dxs = [diff(edge)[1] for edge in edges]
    proddxs = abs(prod(dxs))
    normfactor = sum(normbm) * proddxs
    normbm /= normfactor
    itp = linear_interpolation(edges, normbm, extrapolation_bc=Flat())
    return itp, normbm, bm, edges
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

"""
    localcollisiontime(particledata, densityitp, temperatureitp)
Calculate the Coulomb collision time for an ensemble of particles, where
the local density and temperature at the particle's initial position is used.
"""
function localcollisiontime(particledata, densityitp, temperatureitp)
    x0 = h5_getdataset(particledata, "x0")
    y0 = h5_getdataset(particledata, "y0")
    z0 = h5_getdataset(particledata, "z0")
    charge = h5_getdataset(particledata, "charge")
    mass = h5_getdataset(particledata, "mass")
    n_itp = load_object(densityitp)
    T_itp = load_object(temperatureitp)
    n = [n_itp(x0[i], y0[i], z0[i]) / TraceParticles.m_p for i in eachindex(x0)]
    T = [T_itp(x0[i], y0[i], z0[i]) for i in eachindex(x0)]
    # Assuming the like-particle collisions are the most frequent
    return spitzercollisionaltime.(charge, mass, T, n)
end

"""
    localcriticalvelocity_cgs(
        particledata,
        densityitp,
        temperatureitp,
        electricfielditp;
        coulomb_logarithm=20
    )
Calculate the critical velocity of an ensemble of particles in CGS units,
using the local density, temperature and parallel electric field at the
particle's initial position.
"""
function localcriticalvelocity_cgs(
    particledata,
    densityitp,
    temperatureitp,
    electricfielditp;
    coulomb_logarithm=20
)
    x0 = h5_getdataset(particledata, "x0")
    y0 = h5_getdataset(particledata, "y0")
    z0 = h5_getdataset(particledata, "z0")
    mass = h5_getdataset(particledata, "mass")
    n_itp = load_object(densityitp)
    T_itp = load_object(temperatureitp)
    E_itp = load_object(electricfielditp)
    n = [n_itp(x0[i], y0[i], z0[i]) * 1e-3 / TraceParticles.m_p_cgs for i in eachindex(x0)]
    T = [T_itp(x0[i], y0[i], z0[i]) for i in eachindex(x0)]
    E = [abs(E_itp(1e2x0[i], 1e2y0[i], 1e2z0[i])) for i in eachindex(x0)]
    return criticalvelocity_cgs.(
        n, T, E, mass;
        coulomb_logarithm=coulomb_logarithm
    )
end

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
    getinitialstate(df::DataFrame)
Return the initial state (or initial condition) of particle simulation stored
as a DataFrame.
"""
function getinitialstate(df::DataFrame)
    return df[:, [:x0, :y0, :z0, :vparal0]]
end
function getfinalstate(df::DataFrame)
    return df[:, [:xf, :yf, :zf, :vparalf]]
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
    threshold
)
    saturated_dist = copy(gridded_dist)
    saturated_dist[saturated_dist.<threshold] .= threshold
    dxs = [diff(xx)[1], diff(yy)[1]]
    proddxs = abs(prod(dxs))
    normfactor = sum(saturated_dist) * proddxs
    saturated_dist /= normfactor
    itp = linear_interpolation((xx, yy), saturated_dist, extrapolation_bc=Flat())
    return itp, saturated_dist
end




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
    return df[:,[:x0, :y0, :z0, :vparal0]]
end
function getfinalstate(df::DataFrame)
    return df[:,[:xf, :yf, :zf, :vparalf]]
end
function getmagneticmoments(df::DataFrame)
    return df.magneticmoment
end
function getinitialtimes(df::DataFrame)
    return df.t0
end
function getfinaltimes(df::DataFrame)
    return df.tf
end

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
    if minimum(bm) < 0
        bm .+= abs(minimum(bm))
    end
    # Empty bins have NaN-values. Set them to zero to reflect their values in
    # the probability distribution. Interpolations.jl doesn't like NaN values
    # either.
    normbm = copy(bm)
    normbm[findall(x -> isnan(x), bm)] .= 0
    dxs = [diff(xx)[1], diff(yy)[1]]
    proddxs = abs(prod(dxs))
    normfactor = sum(normbm)*proddxs
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
    saturated_dist[saturated_dist .< threshold] .= threshold
    dxs = [diff(xx)[1], diff(yy)[1]]
    proddxs = abs(prod(dxs))
    normfactor = sum(saturated_dist)*proddxs
    saturated_dist /= normfactor
    itp = linear_interpolation((xx, yy), saturated_dist, extrapolation_bc=Flat())
    return itp, saturated_dist
end


function collect_batches(filename::String, nbatches::Int)
    df = DataFrame(CSV.File(filename * "_1.csv"))
    for i in 2:nbatches
        df_tmp = DataFrame(CSV.File(filename * "_$i.csv"))
        df = vcat(df, df_tmp)
    end
    return df
end


function csv_to_h5(filename, batchnr)
    df = DataFrame(CSV.File(filename * "_$batchnr.csv"))
    filename_h5 = filename * ".h5"
    groupname = "batch_$batchnr"
    h5open(filename_h5, "cw") do h5_file
        h5group = create_group(h5_file, groupname)
        for key in names(df)
            if key == "retcode" || key == "retmsg"
                write_dataset(h5group, key, string.(df[!, key]))
            else
                write_dataset(h5group, key, df[!, key])
            end
        end
    end
end


function h5_getbatch(filename, batchnr)
    df = h5open(filename, "r") do h5_file
        DataFrame(read(h5_file["batch_$batchnr"]))
    end
end


function h5_getdataset(filename, dataset, batchnr)
    data = h5open(filename, "r") do h5_file
        read(h5_file["batch_$batchnr/" * dataset])
    end
end
function h5_getdataset(filename, dataset; batches=[1])
    nbatches = length(batches)
    groupnames = ["batch_$i/" for i in batches]

    data = reduce(vcat,
        h5open(filename, "r") do h5_file
            firstbatch = read(h5_file[groupnames[1] * dataset])
            npart = length(firstbatch)
            data = Vector{Vector{eltype(firstbatch)}}(undef, nbatches)
            data[1] = firstbatch
            for i in 2:nbatches
               data[i] = read(h5_file[groupnames[i] * dataset])
            end
            data
        end
    )
end


function h5_getall(filename; batches=[1])
    nbatches = length(batches)
    groupnames = ["batch_$i" for i in batches]

    df = reduce(vcat,
        h5open(filename, "r") do h5_file
            dfvec = Vector{DataFrame}(undef, nbatches)
            for i in 1:nbatches
                dfvec[i] = DataFrame(read(h5_file[groupnames[i]]))
            end
            dfvec
        end
    )
end


function h5_getenergies(filename; batches=[1], units="eV")
    fname0 = filename * "_gcastates0.h5"
    fnamef = filename * "_gcastatesf.h5"
    e0 = h5_getdataset(fname0, "energy"; batches=batches)
    ef = h5_getdataset(fnamef, "energy"; batches=batches)
    if units == "eV"
        e0 = e0 * J2eV
        ef = ef * J2eV
    end
    return e0, ef
end

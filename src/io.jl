#-------------------------------------------------------------------------------
# Created 11.02.25
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 io.jl
#
#-------------------------------------------------------------------------------
# Contains functions for saving and loading data, and for other I/O operations.
#-------------------------------------------------------------------------------



#____/\_____/\_________________________________________________________________
#
# Saving routines
#

"""
    save_gcastates(filename, fields_itp)
Calculates the GCAStates of the initial and final states of an ensamble of
test particles.
"""
function save_gcastates(
    filename::String,
    fields_itp::Vector{<:AbstractInterpolation};
    components="all",
    dimensionality="2Dxz",
)
    name, extension = splitext(filename)
    if extension == ".csv"
        df = DataFrame(CSV.File(filename))
        dfgca0, dfgcaf = GCAState(
            df,
            fields_itp;
            components=components,
            dimensionality=dimensionality
        )
        CSV.write(filename * "_gcastate0.csv", dfgca0)
        CSV.write(filename * "_gcastatef.csv", dfgcaf)
    elseif extension == ".h5"
        _, nbatches = h5_nbatches(filename)
        h5_sol = h5open(filename, "r")
        h5_gcastates0 = h5open(name * "_gcastates0.h5", "w")
        h5_gcastatesf = h5open(name * "_gcastatesf.h5", "w")
        for i in 1:nbatches
            batchname = "batch_$i"
            batch = read(h5_sol[batchname])
            df = DataFrame(batch)
            dfgca0, dfgcaf = GCAState(
                df,
                fields_itp;
                components=components,
                dimensionality=dimensionality,
            )
            group0 = create_group(h5_gcastates0, batchname)
            groupf = create_group(h5_gcastatesf, batchname)
            for key in names(dfgca0)
                write_dataset(group0, key, dfgca0[!, key])
                write_dataset(groupf, key, dfgcaf[!, key])
            end
        end
        close(h5_sol)
        close(h5_gcastates0)
        close(h5_gcastatesf)
    elseif extension == ""
        error("Filname must include extension.")
    else
        error("File extension not supported.")
    end
end


#____/\_____/\_________________________________________________________________
#
# HDF5 Loading routines
#

"""
    h5_nbatches(filename)
Returns the number of reduction batches in the HDF5 file.
"""
function h5_nbatches(filename)
    h5open(filename, "r") do h5_file
        groups = keys(h5_file)
        nbatches = length(
            findall(g -> occursin(r"batch_\d+", g), groups)
        )
        return 1:nbatches, nbatches
    end
end


"""
    h5_getbatch(filename, batchnr)
Returns batch number `batchnr` as a `DataFrame`.
"""
function h5_getbatch(filename, batchnr)
    h5open(filename, "r") do h5_file
        DataFrame(read(h5_file["batch_$batchnr"]))
    end
end

"""
    h5_getdataset(filename, dataset, batchnr)
Returns dataset `dataset` from batch number `batchnr` in a HDF5-file.

"""
function h5_getdataset(filename, dataset, batchnr)
    h5open(filename, "r") do h5_file
        read(h5_file["batch_$batchnr/"*dataset])
    end
end


"""
    h5_getdataset(filename, dataset; batches=nothing)
Returns the dataset `dataset` multiple batches in a HDF5-file. If `batches` is
not specified, all batches are returned.
"""
function h5_getdataset(filename, dataset; batches=nothing)
    if isnothing(batches)
        batches, nbatches = h5_nbatches(filename)
    else
        nbatches = length(batches)
    end
    groupnames = ["batch_$i/" for i in batches]
    data = reduce(vcat,
        h5open(filename, "r") do h5_file
            firstbatch = read(h5_file[groupnames[1]*dataset])
            npart = length(firstbatch)
            data = Vector{Vector{eltype(firstbatch)}}(undef, nbatches)
            data[1] = firstbatch
            for i in 2:nbatches
                data[i] = read(h5_file[groupnames[i]*dataset])
            end
            data
        end
    )
end


"""
    h5_getall(filename; batches=nothing)
Returns all datasats from a set of batches in a HDF5-file as a `DataFrame`.
If `batches` is not specified, the data from all batches are returned.
"""
function h5_getall(filename; batches=nothing)
    if isnothing(batches)
        batches, nbatches = h5_nbatches(filename)
    else
        nbatches = length(batches)
    end
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


"""
    h5_getenergies(filename; batches=nothing, units="eV")
Returns the initial and final energies of a set of batches in a HDF5-file.
If `batches` is not specified, the data from all batches are returned.
"""
function h5_getenergies(filename; batches=nothing, units="eV")
    expname, _ = splitext(filename)
    if isnothing(batches)
        batches, _ = h5_nbatches(expname * ".h5")
    end
    filename0 = expname * "_gcastates0.h5"
    filenamef = expname * "_gcastatesf.h5"
    e0 = h5_getdataset(filename0, "energy"; batches=batches)
    ef = h5_getdataset(filenamef, "energy"; batches=batches)
    if units == "eV"
        e0 *= J2eV
        ef *= J2eV
    end
    return e0, ef
end


"""
    h5_getinitialstate(filename)
Returns `x0`, `y0`, `z0`, `vparal0`, and `magneticmoment` from a test particle
ensemble stored in a HDF5-file.
"""
function h5_getinitialstate(filename; batches=nothing)
    _, ext = splitext(filename)
    if ext != ".h5"
        error("Filename must be a HDF5 file with extension .h5")
    end
    if isnothing(batches)
        batches, _ = h5_nbatches(filename)
    end
    x0 = h5_getdataset(filename, "x0"; batches=batches)
    y0 = h5_getdataset(filename, "y0"; batches=batches)
    z0 = h5_getdataset(filename, "z0"; batches=batches)
    vparal0 = h5_getdataset(filename, "vparal0"; batches=batches)
    mu0 = h5_getdataset(filename, "magneticmoment"; batches=batches)
    return x0, y0, z0, vparal0, mu0
end


"""
    h5_getfinalstate(filename)
Returns `xf`, `yf`, `zf`, `vparalf`, and `magneticmoment` from a test particle
ensemble stored in a HDF5-file.
"""
function h5_getfinalstate(filename; batches=nothing)
    _, ext = splitext(filename)
    if ext != ".h5"
        error("Filename must be a HDF5 file with extension .h5")
    end
    if isnothing(batches)
        batches, _ = h5_nbatches(filename)
    end
    xf = h5_getdataset(filename, "xf"; batches=batches)
    yf = h5_getdataset(filename, "yf"; batches=batches)
    zf = h5_getdataset(filename, "zf"; batches=batches)
    vparalf = h5_getdataset(filename, "vparalf"; batches=batches)
    mu = h5_getdataset(filename, "magneticmoment"; batches=batches)
    return xf, yf, zf, vparalf, mu
end


#____/\_____/\_________________________________________________________________
#
# Modifying data in HDF5-files
#

"""
    h5_mergebatches!(filename)
Merges all batches in a HDF5-file into a single group called `merged`.
"""
function h5_mergebatches!(filename)
    batches, nbatches = h5_nbatches(filename)
    fields, numparticles = h5open(filename, "r") do fid
        fields = keys(fid["batch_1"])
        numparticles = 0
        for b in batches
            h5group = fid["batch_$b"]
            numparticles += length(h5group[fields[1]])
        end
        return fields, numparticles
    end
    fid = h5open(filename, "cw")
    h5group = create_group(fid, "merged")
    try
        for f in fields
            if f == "timestamp"
                @info "In"
                data = Vector{String}(undef, nbatches)
            else
                data = Vector{eltype(read(fid["batch_1"][f]))}(undef, numparticles)
            end
            i = 1
            for b in batches
                n = length(fid["batch_$b"][f])
                data[i:i+n-1] .= read(fid["batch_$b"][f])
                i += n
            end
            @info "Writing $f"
            write(h5group, f, data)
        end
    catch
        delete_object(fid, "merged")
    end
    close(fid)
end


"""
    h5_replaceparticles!(fname_original, fname_replacement, original_idxs)
Replaces the particles in the HDF5-file `fname_original` with the particles in
`fname_replacement` at the indices `idxs`, which must either be a vector of
integers or a string pointing to a dataset in the replacement file under the
group "metadata".
"""
function h5_replaceparticles!(
    fname_original::String,
    fname_replacement::String;
    original_idxs::String="original_idxs"
)
    idxs = h5open(fname_replacement) do fid
        read(fid["metadata/$original_idxs"])
    end
    h5_replaceparticles!(fname_original, fname_replacement, idxs)
end
function h5_replaceparticles!(
    fname_original::String,
    fname_replacement::String,
    original_idxs::Vector{<:Int}
)
    original = replaceparticles(fname_original, fname_replacement, original_idxs)

    h5open(fname_original, "cw") do fid
        h5group = create_group(fid, "replaced")
        for key in names(original)
            write_dataset(h5group, key, original[!, key])
        end
    end
end


#____/\_____/\_________________________________________________________________
#
# CSV loading routines
#

"""
    csv_getbatches(filename, nbatches)
Reads data from multiple CSV-files and concatenates them into a single
`DataFrame`.
"""
function csv_getbatches(filename::String, nbatches::Int)
    df = DataFrame(CSV.File(filename * "_1.csv"))
    for i in 2:nbatches
        df_tmp = DataFrame(CSV.File(filename * "_$i.csv"))
        df = vcat(df, df_tmp)
    end
    return df
end


"""
    csv_to_h5(filename, batchnr)
Converts a CSV-file to a HDF5-file.
"""
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


#____/\_____/\_________________________________________________________________
#
# Creating interpolation objects from Bifrost input (MHD snapshots)
#

"""
    create_bifrost_interpolators(
        brxp::BifrostExperiment,
        snap::Integer,
        itp_type::Interpolations.InterpolationType,
        itp_bc::Interpolations.BoundaryCondition;
        filename::String
        ;
        units="si",
        auxvariables=["tg", "ne"],
    )
Create interpolation objects from the MHD fields in a Bifrost snapshot and save
them as JLD2-files.
"""
function create_bifrost_itps(
    brxp::BifrostExperiment,
    snap::Integer,
    itp_type::Interpolations.InterpolationType,
    itp_bc::Union{
        Interpolations.BoundaryCondition,
        NTuple{N, Interpolations.BoundaryCondition} where N
    };
    filename::String,
    units="si",
    auxvariables=[],
    normalise=fill(false, length(auxvariables)),
    destagger=fill(false, length(auxvariables)),
)
    if length(auxvariables) != length(normalise)
        throw(ArgumentError(
            "Length of `auxvariables` and `normalise` must be equal."
        ))
    end
    # Create interpolator for the electromagnetic field
    interpolator = get_br_emfield_vecof_interpolators(
        brxp,
        snap,
        itp_type=itp_type,
        itp_bc=itp_bc,
        units=units,
        destagger=true,
    )
    @save string(filename, "_BE.jld2") interpolator

    # Create interpolators for the auxiliary variables
    for (var, norm, dstgr) in zip(auxvariables, normalise, destagger)
        interpolator, _ = get_br_var_interpolator(
            brxp,
            snap,
            var,
            itp_type=itp_type,
            itp_bc=itp_bc,
            units=units,
            normalise=norm,
            destagger=dstgr,
        )
        @save string(filename, "_$(var)_norm$(norm).jld2") interpolator
    end

end

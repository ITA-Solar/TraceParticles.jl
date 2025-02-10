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
    if isnothing(batches)
        batches, _ = h5_nbatches(filename * ".h5")
    end
    filename0 = filename * "_gcastates0.h5"
    filenamef = filename * "_gcastatesf.h5"
    e0 = h5_getdataset(filename0, "energy"; batches=batches)
    ef = h5_getdataset(filenamef, "energy"; batches=batches)
    if units == "eV"
        e0 *= J2eV
        ef *= J2eV
    end
    return e0, ef
end


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
        return batches, fields, numparticles
    end
    fid = h5open(filename, "cw")
    h5group = create_group(fid, "merged")
    for f in fields[fields.!=="timestamp"]
        if f == "timestamp"
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
        write(h5group, f, data)
    end
    close(fid)
end


#____/\_____/\_________________________________________________________________
#
# CSV
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

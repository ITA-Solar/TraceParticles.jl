"""
    print_reduction_overview(batchnr, nbatches, expdir, filename, I)
Prints some meta-information about the current batch.
"""
function print_reduction_overview(batchnr, nbatches, expdir, filename, I)
    println("")
    println("=================================================================",
        "===============")
    println(pwd())
    println("Batch $(batchnr) of $(nbatches).")
    println("")
    println("Saving particles $I to")
    println("--> $expdir")
    println("--> $filename")
end


"""
    print_batch_statistics(batch)
Prints some statistics about the current batch.
"""
function print_batch_statistics(batch)
    println("")
    println("Timestamp:                   $(now())")
    println("Summary statistics:")
    println("Number of particles:         $(size(batch)[1])")
    @printf "Average number of timesteps: %.1f\n" mean(batch.nt)
    @printf "Average end time:            %.1f\n" mean(batch.tf)
    println("Number of maxiters:          $(sum(batch.retcode .== :maxiters))")
    println("Number of OutofDomain:       $(sum(batch.retmsg .== "OutofDomain"))")
    println("Number of GCABreakDown:      $(sum(batch.retmsg .== "GCABreakDown"))")
    println("Number of Relativistic:      $(sum(batch.retmsg .== "Relativistic"))")
    println("=================================================================",
        "===============")
end


"""
    SaveBatchAsCSV(expdir, expname, batchsize, nbatches)
Reduction functor that saves the batch as a DataFrame.
The struct instance is callable and its fields are
- `expdir::String`: The directory where the data is saved.
- `expname::String`: The name of the experiment.
- `batchnr::Int`: The current batch number.
- `batchsize::Int`: The size of the batch.
- `nbatches::Int`: The total number of batches.
"""
mutable struct SaveBatchAsCSV
    datadir::String
    expname::String
    batchnr::Int
    batchsize::Int
    nbatches::Int
    function SaveBatchAsDF(
        expdir::String,
        expname::String,
        batchsize::Int,
        nbatches::Int
    )
        batchnr = 0
        new(expdir, expname, batchnr, batchsize, nbatches)
    end
end
function (self::SaveBatchAsCSV)(u, batch, I)
    # Save the batch as a DataFrame
    self.batchnr += 1
    filename = self.expname * "_$(self.batchnr).csv"
    df = DataFrame(batch)
    expdir = joinpath(self.datadir, self.expname)
    filename = joinpath(expdir, filename)
    print_reduction_overview(self.batchnr, self.nbatches, expdir, filename, I)
    CSV.write(filename, df)
    println("success")
    print_batch_statistics(df)
    return u, false
end


"""
    SaveBatchAsHDF5(expdir, expname, batchsize, nbatches)
Reduction functor that saves the batch as an HDF5 file.
The struct instance is callable and its fields are
- `expdir::String`: The directory where the data is saved.
- `expname::String`: The name of the experiment.
- `suffix::String`: A suffix to the filename.
- `batchnr::Int`: The current batch number.
- `batchsize::Int`: The size of the batch.
- `nbatches::Int`: The total number of batches.
- `metadata::Dict{<:String, <:Any}`: A dictionary of user-specified metadata.
"""
mutable struct SaveBatchAsHDF5
    expdir::String
    expname::String
    suffix::String
    batchnr::Int
    batchsize::Int
    nbatches::Int
    metadata::Dict{<:String, <:Any}
    function SaveBatchAsHDF5(
        expdir::String,
        expname::String,
        batchsize::Int,
        nbatches::Int;
        suffix::String="",
        metadata::Dict{<:String, <:Any}=Dict{String, Any}()
    )
        batchnr = 0
        new(expdir, expname, suffix, batchnr, batchsize, nbatches, metadata)
    end
end
function (self::SaveBatchAsHDF5)(u, batch, I)
    # Save the batch as a DataFrame
    self.batchnr += 1
    filename = get_filename(self)

    print_reduction_overview(
        self.batchnr,
        self.nbatches,
        self.expdir,
        filename,
        I
    )

    df = DataFrame(batch)

    batchname = "batch_$(self.batchnr)"
    # Open data file
    h5open(filename, "cw") do fid
        # Create group for the batch and write the particle data
        h5group = create_group(fid, batchname)
        for key in names(df)
            write(h5group, key, df[!, key])
        end

        # If it is the first batch, write metadata to the file, and create the
        # group for batch specific metadata
        if self.batchnr == 1
            metagroup = create_group(fid, "metadata")
            for (key, value) in self.metadata
                write(metagroup, key, value)
            end
            batchmetagroup = create_group(fid, "batch specific metadata")
        else
            batchmetagroup = fid["batch specific metadata"]
        end
        # Write batch specific metadata
        write(batchmetagroup, batchname, string(now()))
    end
    # Write batch statistics to file?
    println("success")

    print_batch_statistics(df)

    return u, false
end


"""
    get_filename(reduction::SaveBatchAsHDF5)
Returns the filename for the current batch from the reduction functor.
"""
function get_filename(reduction::SaveBatchAsHDF5)
    joinpath(
        reduction.expdir,
        reduction.expname,
        reduction.expname
    ) * reduction.suffix * ".h5"
end

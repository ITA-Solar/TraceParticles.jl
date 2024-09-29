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
    println("=================================================================",
        "===============")
end


mutable struct SaveBatchAsDF
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
function (self::SaveBatchAsDF)(u, batch, I)
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


mutable struct SaveBatchAsHDF5
    datadir::String
    expname::String
    batchnr::Int
    batchsize::Int
    nbatches::Int
    function SaveBatchAsHDF5(
            expdir::String,
            expname::String,
            batchsize::Int,
            nbatches::Int
            )
        batchnr = 0
        new(expdir, expname, batchnr, batchsize, nbatches)
    end
end
function (self::SaveBatchAsHDF5)(u, batch, I)
    # Save the batch as a DataFrame
    self.batchnr += 1
    filename = self.expname * ".h5"
    expdir = joinpath(self.datadir, self.expname)
    filename = joinpath(expdir, filename)

    print_reduction_overview(self.batchnr, self.nbatches, expdir, filename, I)

    df = DataFrame(batch)

    h5open(filename, "cw") do fid
        h5group = create_group(fid, "batch_$(self.batchnr)")
        for key in names(df)
            write(h5group, key, df[!, key])
        end
        write(h5group, "timestamp", string(now()))
    end
    # Write batch statistics to file?
    println("success")

    print_batch_statistics(df)

    return u, false
end

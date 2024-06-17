mutable struct SaveBatchAsDF
    expdir::String
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
    filename = joinpath(self.expdir, filename)
    println("")
    println("Batch $(self.batchnr) of $(self.nbatches).")
    println("")
    println("Saving particles $I to")
    println("--> $(self.expdir)")
    println("--> $(filename)")
    CSV.write(filename, df)
    println("success")
    println("")
    println("Summary statistics:")
    println("Number of particles:         $(length(batch))")
    @printf "Average number of timesteps: %.1f\n" mean(df.nt)
    println("Number of maxiters:          $(sum(df.retcode .== :maxiters))")
    println("Number of OutofDomain:       $(sum(df.retmsg .== "OutofDomain"))")
    println("Number of GCABreakDown:      $(sum(df.retmsg .== "GCABreakDown"))")
    println("")
    println("=================================================================",
        "===============")
    return u, false
end

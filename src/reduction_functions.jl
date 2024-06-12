mutable struct SaveBatchAsDF
    expdir::String
    expname::String
    batchnr::Int
    batchsize::Int
    function SaveBatchAsDF(expdir::String, expname::String, batchsize::Int)
        batchnr = 0
        new(expdir, expname, batchnr, batchsize)
    end
end
function (self::SaveBatchAsDF)(u, batch, I)
    # Save the batch as a DataFrame
    self.batchnr += 1
    filename = self.expname * "_$(self.batchnr).csv"
    df = DataFrame(batch)
    filename = joinpath(self.expdir, filename)
    println("Saving batch $(self.batchnr) (particles $I) to")
    println("--> $(self.expdir)")
    println("--> $(filename)")
    CSV.write(filename, df)
    println("success")
    println("")
    println("Summary statistics:")
    println("Number of particles:         $(length(batch))")
    println("Average number of timesteps: $(mean(df.nt))")
    println("Number of maxiters:          $(sum(df.retcode .== :maxiters))")
    println("Number of OutofDomain:       $(sum(df.retmsg .== "OutofDomain"))")
    println("Number of GCABreakDown:      $(sum(df.retmsg .== "GCABreakDown"))")
    println("")
    println("=================================================================")
    println("")
    println("Running batch $(self.batchnr) of $(self.batchsize).")
    return u, false
end

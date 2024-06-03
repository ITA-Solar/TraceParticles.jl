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
    filename = self.expname * "_batch_$(self.batchnr)_size_$(self.batchsize).csv"
    df = DataFrame(batch)
    filename = joinpath(self.expdir, filename)
    println("=================================================================")
    println("Saving batch $(self.batchnr) (particles $I) to")
    println("--> $(self.expdir)")
    println("--> $(filename)")
    CSV.write(filename, df)
    println("success")
    println("=================================================================")
    return u, false
end

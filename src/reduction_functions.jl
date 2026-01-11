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
    println("--> $filename")
end


"""
    print_batch_statistics(batch)
Prints some statistics about the current batch.
"""
function print_batch_statistics(batch)
    npart = nrow(batch)
    c1 = 100 / npart
    nmaxiters = sum(batch.retcode .== 4)
    noob = sum(batch.terminationcode .== 1)
    ngradient = sum(batch.terminationcode .== 2)
    ncurvature = sum(batch.terminationcode .== 3)
    neratio = sum(batch.terminationcode .== 4)
    nrel = sum(batch.terminationcode .== 5)
    if hasproperty(batch, :nswitches)
        hybswitches = batch.nswitches[batch.nswitches.>0]
    end

    println("")
    println("Timestamp:                   $(now())\n")
    println("                  Summary statistics")
    @printf "Number of particles:         %.1E\n" npart
    @printf "Average number of timesteps: %.1f\n" mean(batch.nt)
    @printf "Median number of timesteps:  %.1f\n" median(batch.nt)
    @printf "Standard deviation of steps: %.1f\n" std(batch.nt)
    @printf "Fewest number of steps:      %i\n" minimum(batch.nt)
    @printf "Most number of steps:        %i\n" maximum(batch.nt)
    @printf "Total number of timesteps:   %.2E\n" sum(batch.nt)
    @printf "Initial time:                %.5f\n" batch.t0[1]
    @printf "Oldest particle:             %.5f\n" maximum(batch.tf)
    @printf "Average end time:            %.5f\n" mean(batch.tf)
    println("Termination statistics")
    @printf "  Nof. maxiters:             %i (%.2f%%)\n" nmaxiters nmaxiters * c1
    @printf "  Nof. OutofBounds:          %i (%.2f%%)\n" noob noob * c1
    @printf "  Nof. Relativistic:         %i (%.2f%%)\n" nrel nrel * c1
    if hasproperty(batch, :nswitches)
        nhybrid = sum(batch.nswitches .> 0)
        maxhyb = sum(@. (batch.nswitches > 0) & (batch.retcode == 4))
        println("Hybrid switch statistics")
        @printf "  Amount that switched:      %i (%.2f%%)\n" nhybrid nhybrid * c1
        if length(hybswitches) > 0
            @printf "  Average nof. switches:     %.4f\n" mean(hybswitches)
            @printf "  Median:                    %i\n" median(hybswitches)
            @printf "  Standard deviation:        %.4f\n" std(hybswitches)
            @printf "  Most nof. switches:        %i\n" maximum(hybswitches)
            @printf "  Least nof. switches:       %i\n" minimum(hybswitches)
            @printf(
                "  MaxIters & switched:       %i (%.0f%% of %i, %.0f%% of MaxIters)\n",
                maxhyb, maxhyb / nhybrid * 100, nhybrid, maxhyb / nmaxiters * 100
            )
        end
    else
        @printf "  Nof. MagneticGradient:     %i (%.0f%%)\n" ngradient ngradient * c1
        @printf "  Nof. MagneticCurvature:    %i (%.0f%%)\n" ncurvature ncurvature * c1
        @printf "  Nof. ParallelEfield:       %i (%.0f%%)\n" neratio neratio * c1
    end
    println("=================================================================",
        "===============")
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
    metadata::Dict{<:String,<:Any}
    verbose::Bool
    function SaveBatchAsHDF5(
        expdir::String,
        expname::String,
        batchsize::Int,
        nbatches::Int;
        suffix::String="",
        metadata::Dict{<:String,<:Any}=Dict{String,Any}(),
        verbose::Bool=true
    )
        batchnr = 0
        new(
            expdir,
            expname,
            suffix,
            batchnr,
            batchsize,
            nbatches,
            metadata,
            verbose
        )
    end
end
function (self::SaveBatchAsHDF5)(u, batch, I)
    # Save the batch as a DataFrame
    self.batchnr += 1
    filename = get_filename(self)

    if self.verbose
        print_reduction_overview(
            self.batchnr,
            self.nbatches,
            self.expdir,
            filename,
            I
        )
    end

    df = DataFrame(batch)

    batchname = "batch_$(self.batchnr)"
    # Open data file
    try
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
        if self.verbose
            println("success")
            print_batch_statistics(df)
        end
    catch e
        try
            CSV.write("./tp_panic.csv", df)
            @warn "Error writing batch to HDF5-file. The batch has been written " *
                  "as a `DataFrame` in 'tp_panic.csv'."
        catch e2
            @warn "Writing panic-file failed: $e2"
        finally
            @error e
        end
    end

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

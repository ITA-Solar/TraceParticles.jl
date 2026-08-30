"""
    reduction_overview(batchnr, nbatches, expdir, filename, I)
Returns some meta-information about the current batch.
"""
function reduction_overview(batchnr, nbatches, expdir, filename, I)
    return """
========================================================================
Timestamp: $(now())

Batch $(batchnr) of $(nbatches).

Saved particles $I to
--> $(joinpath(pwd(), filename))

"""
end


"""
    batch_statistics(batch)
Returns some statistics about the current batch.
"""
function batch_statistics(batch)
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

    msg = """
        Summary statistics
$(@sprintf "Number of particles:         %.1E" npart)
$(@sprintf "Average number of timesteps: %.1f" mean(batch.nt))
$(@sprintf "Median number of timesteps:  %.1f" median(batch.nt))
$(@sprintf "Standard deviation of steps: %.1f" std(batch.nt))
$(@sprintf "Fewest number of steps:      %i" minimum(batch.nt))
$(@sprintf "Most number of steps:        %i" maximum(batch.nt))
$(@sprintf "Total number of timesteps:   %.2E" sum(batch.nt))
$(@sprintf "Average start time           %.3E" mean(batch.t0))
$(@sprintf "Average end time:            %.3E" mean(batch.tf))

        Termination statistics
$(@sprintf "Nof. maxiters:               %i (%.2f%%)" nmaxiters nmaxiters * c1)
$(@sprintf "Nof. OutofBounds:            %i (%.2f%%)" noob noob * c1)
$(@sprintf "Nof. Relativistic:           %i (%.2f%%)" nrel nrel * c1)
"""

    if hasproperty(batch, :nswitches)
        nhybrid = sum(batch.nswitches .> 0)
        maxhyb = sum(@. (batch.nswitches > 0) & (batch.retcode == 4))
        msg *= "\n        Hybrid switch statistics\n"
        msg *= @sprintf "Amount that switched:        %i (%.2f%%)\n" nhybrid nhybrid * c1
        if length(hybswitches) > 0
            msg *= @sprintf "Average nof. switches:       %.4f\n" mean(hybswitches)
            msg *= @sprintf "Median:                      %i\n" median(hybswitches)
            msg *= @sprintf "Standard deviation:          %.4f\n" std(hybswitches)
            msg *= @sprintf "Most nof. switches:          %i\n" maximum(hybswitches)
            msg *= @sprintf "Least nof. switches:         %i\n" minimum(hybswitches)
            msg *= @sprintf(
                "MaxIters & switched:         %i (%.0f%% of %i, %.0f%% of MaxIters)\n",
                maxhyb, maxhyb / nhybrid * 100, nhybrid, maxhyb / nmaxiters * 100
            )
        end
    else
        msg *= @sprintf "Nof. MagneticGradient:       %i (%.0f%%)\n" ngradient ngradient * c1
        msg *= @sprintf "Nof. MagneticCurvature:      %i (%.0f%%)\n" ncurvature ncurvature * c1
        msg *= @sprintf "Nof. ParallelEfield:         %i (%.0f%%)\n" neratio neratio * c1
    end
    msg *= "=================================================================" *
        "============="
    return msg
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
@kwdef mutable struct SaveBatchAsHDF5{
    S<:String,
    I<:Integer,
    B<:Bool,
    T1
}
    expdir::S
    expname::S
    batchsize::I
    nbatches::I
    batchnr::I = 0
    suffix::S = ""
    verbose::B = true
    metadata::T1 = Dict{String,String}()
end
function (self::SaveBatchAsHDF5)(u, batch, I)
    # Save the batch as a DataFrame
    self.batchnr += 1
    filename = get_filename(self)

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
        if self.verbose
            @info reduction_overview(
                self.batchnr,
                self.nbatches,
                self.expdir,
                filename,
                I
            ) * batch_statistics(df)
        end
    catch e
        try
            JLD2.save("./tp_panic.jld2", "batch", df)
            @warn "Error writing batch to HDF5-file. The batch has been written " *
                "as a `DataFrame` in 'tp_panic.jld2'."
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

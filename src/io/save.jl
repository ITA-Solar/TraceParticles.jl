#-------------------------------------------------------------------------------
# Created 25.01.24
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 save.jl
#
#-------------------------------------------------------------------------------
# Contains functions for saving data.
#-------------------------------------------------------------------------------


"""
    save_gcastates(filename, nbatches, fields_itp)
Calculates the GCAStates of the initial and final states of an ensamble of
test particles.
"""
function save_gcastates(
    filename::String,
    fields_itp::Vector{<:AbstractInterpolation}
)
    name, extension = splitext(filename)
    if extension == ".csv"
        for i in 1:nbatches
            df = DataFrame(CSV.File(filename * "_$i.csv"))
            dfgca0, dfgcaf = GCAState(df, fields_itp, components="all")
            CSV.write(filename * "_gcastate0_$i.csv", dfgca0)
            CSV.write(filename * "_gcastatef_$i.csv", dfgcaf)
        end
    elseif extension == ".h5"
        _, nbatches = tp.h5_nbatches(filename)
        h5_sol = h5open(filename, "r")
        h5_gcastates0 = h5open(name * "_gcastates0.h5", "w")
        h5_gcastatesf = h5open(name * "_gcastatesf.h5", "w")
        for i in 1:nbatches
            batchname = "batch_$i"
            batch = read(h5_sol[batchname])
            df = DataFrame(batch)
            dfgca0, dfgcaf = GCAState(df, fields_itp, components="all")
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

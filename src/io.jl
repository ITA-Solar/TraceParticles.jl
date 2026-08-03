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
    write_dataset(group, dict)
Write a dictionary to a HDF5 group.
"""
function write_dict(h5group::HDF5.Group, dict::Dict{String,Any})
    for (key, value) in dict
        write_dataset(h5group, key, value)
    end
end

"""
    save_gcastates(filename, fields_itp)
Calculates the GCAStates of the initial and final states of an ensamble of
test particles.
"""
function save_gcastates(
    filename::String,
    electromagneticfield::Any;
    components="all",
    verbose=false,
)
    name, extension = splitext(filename)
    if extension == ".h5"
        nbatches = length(h5_getbatchnames(filename))
        h5_sol = h5open(filename, "r")
        name0 = name * "_gcastates0.h5"
        namef = name * "_gcastatesf.h5"
        h5_gcastates0 = h5open(name0, "w")
        h5_gcastatesf = h5open(namef, "w")
        try
            for i in 1:nbatches
                batchname = "batch_$i"
                batch = read(h5_sol[batchname])
                df = DataFrame(batch)
                if verbose
                    @info "Calculating GCA states for batch $i"
                end
                dfgca0, dfgcaf = GCAState(
                    df,
                    electromagneticfield;
                    components=components,
                )
                group0 = create_group(h5_gcastates0, batchname)
                groupf = create_group(h5_gcastatesf, batchname)
                for key in names(dfgca0)
                    write_dataset(group0, key, dfgca0[!, key])
                    write_dataset(groupf, key, dfgcaf[!, key])
                end
            end
        catch e
            rm(name0)
            rm(namef)
            @error "While saving GCA states:" exception = (e, catch_backtrace())
        end
        try
            groupname = "metadata"
            metadata = read(h5_sol[groupname])
            metagroup0 = create_group(h5_gcastates0, groupname)
            metagroupf = create_group(h5_gcastatesf, groupname)
            write_dict(metagroup0, metadata)
            write_dict(metagroupf, metadata)
        catch e
            @warn "While copying metadata:" exception = (e, catch_backtrace())
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

"""
    save_energy(filename, fields_itp)
Calculates the kinetic energy of the initial and final states of an ensamble of
test particles.
"""
function save_energy(
    filename::String,
    electromagneticfield::Any;
    eom::Function,
    components="all",
    verbose=false,
)
    name, extension = splitext(filename)
    if extension == ".h5"
        nbatches = length(h5_getbatchnames(filename))
        h5_sol = h5open(filename, "r")
        name0 = name * "_e0.h5"
        namef = name * "_ef.h5"
        h5_e0 = h5open(name0, "w")
        h5_ef = h5open(namef, "w")
        try
            for i in 1:nbatches
                batchname = "batch_$i"
                batch = read(h5_sol[batchname])
                df = DataFrame(batch)
                if verbose
                    @info "Calculating energy states for batch $i"
                end
                npart = nrow(df)
                e0 = Vector{Float64}(undef, npart)
                ef = Vector{Float64}(undef, npart)
                compute_energies!(e0, ef, eom, df, electromagneticfield)
                group0 = create_group(h5_e0, batchname)
                groupf = create_group(h5_ef, batchname)
                write_dataset(group0, "energy", e0)
                write_dataset(groupf, "energy", ef)
            end
        catch e
            rm(name0)
            rm(namef)
            @error "While saving energies:" exception = (e)#, catch_backtrace())
        end
        try
            groupname = "metadata"
            metadata = read(h5_sol[groupname])
            metagroup0 = create_group(h5_e0, groupname)
            metagroupf = create_group(h5_ef, groupname)
            write_dict(metagroup0, metadata)
            write_dict(metagroupf, metadata)
        catch e
            @warn "While copying metadata:" exception = (e, catch_backtrace())
        end
        close(h5_sol)
        close(h5_e0)
        close(h5_ef)
    elseif extension == ""
        error("Filname must include extension.")
    else
        error("File extension not supported.")
    end
end

function compute_energies!(
    e0::AbstractVector,
    ef::AbstractVector,
    eom::typeof(guidingcentreapproximation!),
    df::DataFrame,
    electromagneticfield::Any
)
    for i in eachindex(e0)
        e0[i] = kineticenergy(
            [
                df[i, :x0],
                df[i, :y0],
                df[i, :z0],
                df[i, :vparal0],
            ],
            df[i, :charge],
            df[i, :mass],
            df[i, :magneticmoment],
            electromagneticfield,
            df[i, :t0],
        )
        ef[i] = kineticenergy(
            [
                df[i, :xf],
                df[i, :yf],
                df[i, :zf],
                df[i, :vparalf],
            ],
            df[i, :charge],
            df[i, :mass],
            df[i, :magneticmoment],
            electromagneticfield,
            df[i, :tf],
        )
    end
end

function compute_energies!(
    e0::AbstractVector,
    ef::AbstractVector,
    eom::typeof(hybridgcafo!),
    df::DataFrame,
    electromagneticfield::Any
)
    for i in eachindex(e0)
        if df[i, :initialeomid] == 1
            e0[i] = kineticenergy(
                [
                    df[i, :x0],
                    df[i, :y0],
                    df[i, :z0],
                    df[i, :vx0],
                ],
                df[i, :charge],
                df[i, :mass],
                df[i, :initialmagneticmoment],
                electromagneticfield,
                df[i, :t0],
            )
        elseif df[i, :initialeomid] == 2
            e0[i] = kineticenergy(
                [
                    df[i, :vx0],
                    df[i, :vy0],
                    df[i, :vz0],
                ],
                df[i, :mass],
            )
        end
        if df[i, :eomid] == 1
            ef[i] = kineticenergy(
                [
                    df[i, :xf],
                    df[i, :yf],
                    df[i, :zf],
                    df[i, :vxf],
                ],
                df[i, :charge],
                df[i, :mass],
                df[i, :magneticmoment],
                electromagneticfield,
                df[i, :tf],
            )
        elseif df[i, :eomid] == 2
            ef[i] = kineticenergy(
                [
                    df[i, :vxf],
                    df[i, :vyf],
                    df[i, :vzf],
                ],
                df[i, :mass],
            )
        end
    end
end
#____/\_____/\_________________________________________________________________
#
# HDF5 Loading routines
#

"""
    h5_getbatchnames(filename)
Returns the HDF5 group names of the reduction batches.
"""
function h5_getbatchnames(filename)
    h5open(filename, "r") do h5_file
        groups = keys(h5_file)
        batchnumbers = findall(g -> occursin(r"batch_\d+", g), groups)
        return groups[batchnumbers]
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
    h5_getdataset(filename, dataset)
    h5_getdataset(filename, dataset, h5group::String)
    h5_getdataset(filename, dataset, batchnr::Int)
    h5_getdataset(filename, dataset, batches::Vector{<:Int}
    h5_getdataset(filename, dataset, batches::Vector{String}
Returns dataset `dataset` from all batches, the HDF5-group `h5group`,
`batchnr` or range of `batches`.
"""
function h5_getdataset(
    filename::String,
    dataset::Union{Symbol, String},
    h5_group::String
)
    dataset = string(dataset)
    h5open(filename, "r") do h5_file
        read(h5_file[h5_group][dataset])
    end
end
function h5_getdataset(
    filename::String,
    dataset::Union{Symbol, String},
    batchnr::Int
)
    h5_getdataset(filename, dataset, "batch_$batchnr")
end
function h5_getdataset(
    filename::String,
    dataset::Union{Symbol, String}
)
    batches = h5_getbatchnames(filename)
    h5_getdataset(filename, dataset, batches)
end
function h5_getdataset(
    filename::String,
    dataset::Union{Symbol, String},
    batches::Vector{<:Int}
)
    return h5_getdataset(filename, dataset, ["batch_$i/" for i in batches])
end
function h5_getdataset(
    filename::String,
    dataset::Union{Symbol, String},
    batches::Vector{String}
)
    dataset = string(dataset)
    nbatches = length(batches)
    data = reduce(vcat,
        h5open(filename, "r") do h5_file
            firstbatch = read(h5_file[batches[1]*dataset])
            npart = length(firstbatch)
            data = Vector{Vector{eltype(firstbatch)}}(undef, nbatches)
            data[1] = firstbatch
            for i in 2:nbatches
                data[i] = read(h5_file[batches[i]*dataset])
            end
            data
        end
    )
end
function h5_getdataset(filename::String, datasets::AbstractVector, args...)
    df = DataFrame()
    for var in datasets
        df[!, var] = h5_getdataset(filename, var, args...)
    end
    return df
end

"""
    h5_getall(filename)
    h5_getall(filename, group::String)
    h5_getall(filename, batchnr::Int)
    h5_getall(filename, batches::AbstractVector)
Return all datasets in the HDF5-file `filename` as a `DataFrame`.
"""
function h5_getall(filename::String, group::String)
    return h5open(filename, "r") do h5_file
        DataFrame(read(h5_file[group]))
    end
end
function h5_getall(filename::String, batchnr::Int)
    return h5_getall(filename, "batch_$batchnr")
end
function h5_getall(filename)
    batches = h5_getbatchnames(filename)
    return h5_getall(filename, batches)
end
function h5_getall(filename::String, batches::Vector{<:Int})
    nbatches = length(batches)
    groupnames = ["batch_$i" for i in batches]
    return h5_getall(filename, groupnames)
end
function h5_getall(filename::String, batches::Vector{String})
    nbatches = length(batches)
    return reduce(vcat,
        h5open(filename, "r") do h5_file
            dfvec = Vector{DataFrame}(undef, nbatches)
            for i in 1:nbatches
                dfvec[i] = DataFrame(read(h5_file[batches[i]]))
            end
            dfvec
        end
    )
end



"""
    h5_getenergies_gca(filename, args...; units="eV")
Returns the initial and final energies of a set of batches in a HDF5-file.
If `batches` is not specified, the data from all batches are returned.
"""
function h5_getenergies_gca(filename, args...; units="eV")
    expname, _ = splitext(filename)
    filename0 = expname * "_gcastates0.h5"
    filenamef = expname * "_gcastatesf.h5"
    e0 = h5_getdataset(filename0, "energy", args...)
    ef = h5_getdataset(filenamef, "energy", args...)
    if units == "eV"
        e0 *= J2eV
        ef *= J2eV
    end
    return e0, ef
end

"""
    h5_getenergies(filename, args...; units="eV")
Returns the initial and final energies of a set of batches in a HDF5-file.
If `batches` is not specified, the data from all batches are returned.
"""
function h5_getenergies(filename, args...; units="eV")
    expname, _ = splitext(filename)
    filename0 = expname * "_e0.h5"
    filenamef = expname * "_ef.h5"
    e0 = h5_getdataset(filename0, "energy", args...)
    ef = h5_getdataset(filenamef, "energy", args...)
    if units == "eV"
        e0 *= J2eV
        ef *= J2eV
    end
    return e0, ef
end

"""
    h5_getinitialstate(filename, args...)
Returns `x0`, `y0`, `z0`, `vx0`, `vy0`, `vz0`, and `magneticmoment` from a test
particle ensemble stored in a HDF5-file.
"""
function h5_getinitialstate(filename, args...)
    x0 = h5_getdataset(filename, "x0", args...)
    y0 = h5_getdataset(filename, "y0", args...)
    z0 = h5_getdataset(filename, "z0", args...)
    vx0 = h5_getdataset(filename, "vx0", args...)
    vy0 = h5_getdataset(filename, "vy0", args...)
    vz0 = h5_getdataset(filename, "vz0", args...)
    mu0 = h5_getdataset(filename, "magneticmoment", args...)
    return x0, y0, z0, vx0, vy0, vz0, mu0
end

"""
    h5_getinitialstate_gca(filename, args...)
Returns `x0`, `y0`, `z0`, `vparal0`, and `magneticmoment` from a test particle
ensemble using GCA, stored in a HDF5-file.
"""
function h5_getinitialstate_gca(filename, args...)
    x0 = h5_getdataset(filename, "x0", args...)
    y0 = h5_getdataset(filename, "y0", args...)
    z0 = h5_getdataset(filename, "z0", args...)
    vparal0 = h5_getdataset(filename, "vparal0", args...)
    mu0 = h5_getdataset(filename, "magneticmoment", args...)
    return x0, y0, z0, vparal0, mu0
end

"""
    h5_getfinalstate(filename, args...)
Returns `xf`, `yf`, `zf`, `vxf`, `vyf`, `vzf`, and `magneticmoment` from a test
particle ensemble stored in a HDF5-file.
"""
function h5_getfinalstate(filename, args...)
    xf = h5_getdataset(filename, "xf", args...)
    yf = h5_getdataset(filename, "yf", args...)
    zf = h5_getdataset(filename, "zf", args...)
    vxf = h5_getdataset(filename, "vxf", args...)
    vyf = h5_getdataset(filename, "vyf", args...)
    vzf = h5_getdataset(filename, "vzf", args...)
    mu = h5_getdataset(filename, "magneticmoment", args...)
    return xf, yf, zf, vxf, vyf, vzf, mu
end

"""
    h5_getfinalstate_gca(filename, args...)
Returns `xf`, `yf`, `zf`, `vparalf`, and `magneticmoment` from a test particle
ensemble using GCA, stored in a HDF5-file.
"""
function h5_getfinalstate_gca(filename, args...)
    xf = h5_getdataset(filename, "xf", args...)
    yf = h5_getdataset(filename, "yf", args...)
    zf = h5_getdataset(filename, "zf", args...)
    vparalf = h5_getdataset(filename, "vparalf", args...)
    mu = h5_getdataset(filename, "magneticmoment", args...)
    return xf, yf, zf, vparalf, mu
end

"""
    h5_getobservables(
        fname::String;
        obs_fname = joinpath(splitdir(fname, "observables.h5"))
    )
Load and return observables as a `DataFrame`.
"""
function h5_getobservables(
    fname::String;
    obs_fname=joinpath(splitdir(fname)[1], "observables.h5")
)
    h5open(obs_fname) do fid
        return DataFrame(read(fid))
    end
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
    batches = h5_getbatchnames(filename)
    nbatches = length(batches)
    fields, numparticles = h5open(filename, "r") do fid
        fields = keys(fid[batches[1]])
        numparticles = 0
        for groupid in batches
            h5group = fid[groupid]
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
                data = Vector{eltype(read(fid[batches[1]][f]))}(undef, numparticles)
            end
            i = 1
            for groupid in batches
                n = length(fid[groupid][f])
                data[i:i+n-1] .= read(fid[groupid][f])
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

"""
    h5_adddataset!(fname, dataset, key)
Adds a dataset to all batches in a HDF5-file.
"""
function h5_adddataset!(
    fname::String,
    dataset::AbstractVector,
    key::String,
)
    h5open(fname, "r+") do fid
        batches = h5_getbatchnames(fname)
        i = 1
        for b in batches
            group = fid[b]
            n = length(first(group))
            write(group, key,
                dataset[i:i+n-1])
            i += n
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
        filename::String,
        units::String,
        variables::Vector{String}
    )
Create interpolation objects from the MHD variables in a Bifrost snapshot
and save them as JLD2-files.

### Keyword arguments
- `normalise::Vector{Bool}`
- `destagger::Vector{Bool}`
- `xzinterpolation::Bool`: If true, wrap 2D interpolators in a struct that
    allows it to be called with three coordinates, using only the first and
    third coordinate. That is; itp(x,y,z) = itp(x,z).
"""
function create_bifrost_itps(
    brxp::BifrostExperiment,
    snaps::Union{Integer,AbstractVector},
    itp_type::Interpolations.InterpolationType,
    itp_bc::Union{
        Interpolations.BoundaryCondition,
        NTuple{N,Interpolations.BoundaryCondition} where N
    },
    filename::String,
    units,
    variables;
    normalise=fill(false, length(variables)),
    destagger=fill(false, length(variables)),
    xzinterpolation=false,
    deallocateafterwritten=false,
    kwargs...
)
    if length(variables) != length(normalise) != length(destagger)
        throw(ArgumentError(
            "Length of `variables` and `normalise` and `destagger` must be equal."
        ))
    end
    tstamp = length(snaps) == 1 ?
             "_t$snaps" : "_t$(first(snaps))-$(last(snaps))"
    if length(snaps) == 1
        if xzinterpolation
            wrapper = StaticXZInterpolation
        else
            wrapper = StaticInterpolation
        end
    elseif xzinterpolation
        wrapper = XZInterpolation
    else
        wrapper = nothing
    end
    # Create interpolator for the electromagnetic field
    if "BE" ∈ variables
        interpolator = get_br_emfield_vecof_interpolators(
            brxp,
            snaps;
            itp_type=itp_type,
            itp_bc=itp_bc,
            units=units,
            destagger=true,
            kwargs...
        )
        if !isnothing(wrapper)
            interpolator = [wrapper(itp) for itp in interpolator]
        end
        @save string(filename, tstamp, "_BE.jld2") interpolator
        # Remove "BE" from the variables
        mask = variables .!= "BE"
        variables = variables[mask]
        normalise = normalise[mask]
        destagger = destagger[mask]
        # For very large snapshots, it is a good idea to deallocate the
        # `interpolator` variable after loading the electromangetic field.
        # In this way, we do not eat up more memory when loading more
        # variables.
        if deallocateafterwritten
            interpolator = nothing
            GC.gc()
        end
    end

    # Create interpolators for the variables
    for (var, norm, dstgr) in zip(variables, normalise, destagger)
        interpolator, _ = get_br_var_interpolator(
            brxp,
            snaps,
            var;
            itp_type=itp_type,
            itp_bc=itp_bc,
            units=units,
            normalise=norm,
            destagger=dstgr,
            kwargs...
        )
        if !isnothing(wrapper)
            interpolator = wrapper(interpolator)
        end
        @save string(
            filename, tstamp, "_$(var)_norm$(norm).jld2"
        ) interpolator
        if deallocateafterwritten
            interpolator = nothing
            GC.gc()
        end
    end
end

#____/\_____/\_________________________________________________________________
#
# Creating `ODEProblem` arguments
#
"""
    create_diffeq_ic(
        ic_file::string,
        tf::Number,
        emfield::Any;
        f=hybridgcafo!
    )
Construct `u0`, `tspan` and `p` (parameters) for necessary for creating an
`ODEProblem`. The type of `p` depends on the value of `f`.

## Arguments
-`ic_file` is assumed to be a JLD2 file containing
the initial positions, initial velocities, initial times, charge, mass, magnetic
moment, statistical weight, and number of rejections (from sampling the condition)
for the particles to simulate.

-`f` is a function representing the equations of motion for the particles. If the
function is `hybridgcafo!`, then it is assumed that `ic_file` contains a field
`eomid` too, that reflects which state of the hybrid scheme the particles are
initialised in.

-`tf` is a vector of the final time of the particles.

-`emfield` is the callable object that returns the electric and magnetic field at
an arbitrary position. It is used to create the parameters for the `ODEProblem`
but should not be duplicated per particle as it may be large.
"""
function create_diffeq_ic(
    icfile::String,
    tf::Number,
    emfield::Any;
    f::Function=hybridgcafo!
)
    return h5open(icfile) do fid
        x0 = read(fid, "x0")
        y0 = read(fid, "y0")
        z0 = read(fid, "z0")
        vx0 = read(fid, "vx0")
        vy0 = read(fid, "vy0")
        vz0 = read(fid, "vz0")
        t0 = read(fid, "t0")
        charge = read(fid, "charge")
        mass = read(fid, "mass")
        magneticmoment = read(fid, "magneticmoment")
        weight = read(fid, "weight")
        nrejections = read(fid, "nrejections")
        npart = length(x0)
        t_precision = eltype(tf)
        u0 = Vector{Vector{Float64}}(undef, npart)
        tspan = Vector{Tuple{t_precision,t_precision}}(undef, npart)
        for i in 1:npart
            u0[i] = [
                x0[i],
                y0[i],
                z0[i],
                vx0[i],
                vy0[i],
                vz0[i]
            ]
            tspan[i] = (t0[i], tf)
        end
        if f == hybridgcafo!
            eomid = read(fid, "eomid")
            params = Vector{HybridParams}(undef, npart)
            for i in 1:npart
                params[i] = HybridParams(
                    charge=charge[i],
                    mass=mass[i],
                    electromagneticfield=emfield,
                    magneticmoment=magneticmoment[i],
                    weight=weight[i],
                    nrejections=nrejections[i],
                    eomid=eomid[i]
                )
            end
        elseif f == guidingcentreapproximation!
            params = Vector{GCAParams}(undef, npart)
            for i in 1:npart
                params[i] = GCAParams(
                    charge=charge[i],
                    mass=mass[i],
                    electromagneticfield=emfield,
                    magneticmoment=magneticmoment[i],
                    weight=weight[i],
                    nrejections=nrejections[i],
                )
            end
        else
            throw(ArgumentError("No implementation for `f`=$f"))
        end
        return u0, tspan, params
    end
end

function create_expdir(datadir::String, expname::String)
    # Create data-directory if not already present
    try
        mkdir(datadir)
    catch
    end

    # Create experiment-directory. Prompt for deletion if already present.
    expdir = joinpath(datadir, expname)
    try
        mkdir(expdir)
    catch
        println("The directory $expdir already exists.
         Do you want to delete it and re-run? y/n")
        if readline() == "y"
            rm(expdir, recursive=true)
            mkdir(expdir)
        else
            error("Exiting.")
        end
    end
end

function create_paramsbackup(datadir::String, expname::String)
    # Copy parameter-file to experiment-directory, giving it a name equal
    # to the experiment-name. This only work if the parameter file is given
    # as a command line argument.
    expdir = joinpath(datadir, expname)
    try
        paramsbackup = expdir * "/" * expname * ".jl"
        cp(PROGRAM_FILE, paramsbackup)
    catch
        @warn "Unable to copy parameter file to experiment directory."
    end
end

function logger2fileandstderr(
    loglevel::Logging.LogLevel,
    filepath::String,
)
    return TeeLogger(
        ConsoleLogger(stderr, loglevel),
        MinLevelLogger(
            FileLogger(filepath),
            loglevel
        )
    )
end

#-------------------------------------------------------------------------------
# Created 03.01.24 
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 io/bifrost_input.jl
#
#-------------------------------------------------------------------------------
# Contains functions for reading Bifrost snapshots and creating interpoaltion
# objects from the snapshot-fields. Depends on BifrostTools.jl.
#-------------------------------------------------------------------------------
function get_br_var_interpolator(
    brxp::BifrostExperiment,
    snaps::Union{Integer,AbstractVector},
    variable::String,
    ;
    itp_type=BSpline(Linear()),
    itp_bc=Flat(),
    units="si",
    normalise=false,
    destagger=true,
    kwargs...
)
    if variable == "ne"
        if destagger
            error("Electron density is not staggered. Set `destagger=false`")
        else
            var = get_electron_density(brxp, snaps; units=units, kwargs...)
        end
    elseif variable == "eparal"
        var = get_eparal(brxp, snaps; units=units, destagger=destagger, kwargs...)
    else
        var = get_var(
            brxp,
            snaps,
            variable;
            units=units,
            destagger=destagger,
            kwargs...
        )
    end
    var = stack(var)
    var = squeeze(var)
    xaxis, yaxis, zaxis = get_axes(brxp, units=units)
    if haskey(kwargs, :slicex)
        xaxis = xaxis[kwargs[:slicex]]
    end
    if haskey(kwargs, :slicey)
        yaxis = yaxis[kwargs[:slicey]]
    end
    if haskey(kwargs, :slicez)
        zaxis = zaxis[kwargs[:slicez]]
    end
    br_axes = squeeze( (xaxis, yaxis, zaxis) )
    if snaps isa AbstractVector

        times = get_var(brxp, snaps, "t"; units=units)
        br_axes = (br_axes..., times)
    end
    if normalise
        if !(variable in ("ne", "r", "tg"))
            error("Normalisation of vector components not implemented")
        end
        @warn "Normalising routine assumes uniform axes"
        dxvec = [diff(axes)[1] for axes in br_axes]
        dv = prod(dxvec)
        normfactor = sum(var) * dv
        var = var / normfactor
    end
    if itp_type == BSpline(Linear())
        itp_func = linear_interpolation
    else
        range_axes = Vector{StepRangeLen}(undef, length(br_axes))
        for i in eachindex(br_axes)
            gridspacing = diff(br_axes[i])
            if all(gridspacing[1] .== gridspacing)
                range_axes[i] = range(br_axes[i][1], br_axes[i][end], length(br_axes[i]))
            end
        end
        mask = [isassigned(range_axes, i) for i in eachindex(range_axes)]
        if all(mask)
            br_axes = Tuple(range_axes)
        else
            @warn "Non-uniform grid. Creating uniform axes " *
                  "using range."
            idx = findall(x -> !x, mask)
            for i in idx
                range_axes[i] = range(start=br_axes[i][1], stop=br_axes[i][end], length=length(br_axes[i]))
            end
            br_axes = Tuple(range_axes)
        end
        itp_func = cubic_spline_interpolation
        #var_itp = interpolate(var, itp_type)
        #var_itp = extrapolate(scale(var_itp, br_axes...), itp_bc)
    end
    var_itp = itp_func(br_axes, var, extrapolation_bc=itp_bc)
    return var_itp, br_axes
end


function get_br_emfield_vecof_interpolators(args...; kwargs...)
    vars = ("bx", "by", "bz", "ex", "ey", "ez")
    return [get_br_var_interpolator(args..., var; kwargs...)[1] for var in vars]
end

function get_eparal(
    brxp::BifrostExperiment,
    snap::Integer,
    ;
    units=units,
    destagger=true,
    kwargs...
)
    bx = get_var(brxp, snap, "bx"; units=units, destagger=destagger, kwargs...)
    by = get_var(brxp, snap, "by"; units=units, destagger=destagger, kwargs...)
    bz = get_var(brxp, snap, "bz"; units=units, destagger=destagger, kwargs...)
    ex = get_var(brxp, snap, "ex"; units=units, destagger=destagger, kwargs...)
    ey = get_var(brxp, snap, "ey"; units=units, destagger=destagger, kwargs...)
    ez = get_var(brxp, snap, "ez"; units=units, destagger=destagger, kwargs...)
    B = @. sqrt(bx^2 + by^2 + bz^2)
    eparal = @. (ex * bx + ey * by + ez * bz) / B
    return eparal
end

"""
    normfactor_2Duniformmesh(
        brxp, snap, var; xlim=[-Inf, Inf], zlim=[-Inf, Inf], units="si"
    )
Calculate the normalisation factor of a Bifrost scalar field.
"""
function normfactor_2Duniformmesh(
    brxp::BifrostExperiment,
    snap::Integer,
    var::String,
    ;
    xlim=[-Inf, Inf],
    zlim=[-Inf, Inf],
    units="si",
)
    if var == "ne"
        var = get_electron_density(brxp, snap; units=units)
    else
        var = get_var(brxp, snap, var; units=units, destagger=true)
    end

    # Normalise the variable over the full domain
    axes = get_axes(brxp, units=units)
    axes = squeeze(axes)
    var /= normfactor_2Duniformmesh(var[:, 1, :], axes)

    # Normalise the variable over the limits.
    if units == "si"
        axes = (1e6brxp.mesh.x, 1e6brxp.mesh.z)
    else
        axes = (brxp.mesh.x, brxp.mesh.z)
    end
    var = var[:, 1, :]
    return normfactor_2Duniformmesh(var, axes; xlim=xlim, ylim=zlim)
end


#____/\_____/\_________________________________________________________________
#
# Creating interpolation objects from Bifrost input (MHD snapshots)
#
"""
    create_bifrost_itps(
        brxp::BifrostExperiment,
        snaps::Union{Integer,AbstractVector},
        itp_type::Interpolations.InterpolationType,
        itp_bc::Interpolations.BoundaryCondition,
        variables::Vector{String};
        units::String,
        name=missing,
        kwargs...
    )
Create interpolation objects from the MHD variables in a Bifrost snapshot
and save them as JLD2-files.

Pass `"BE"` in `variables` to get the 6-element electromagnetic field vector
`[bx, by, bz, ex, ey, ez]` that `emfield_file` expects.

Each file is named `<name>_t<snap0>-t<snapf>_<var>_<units>_<itp>.itp.jld2`,
where `<itp>` is `linear` or `cubicspline`.

### Keyword arguments
- `units::String`: unit system the variables are converted to, e.g. `"si"`.
- `name`: prefix of the generated filenames, including any directory. Omitted
    from the name when left `missing`.
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
    variables;
    units::String,
    name=missing,
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
             "t$snaps" : "t$(first(snaps))-t$(last(snaps))"
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

    itpstr =
        if itp_type == BSpline(Linear())
            "linear"
        elseif itp_type == BSpline(Cubic())
            "cubicspline"
        else
            missing
        end
    postfix = join(skipmissing([units, itpstr]), "_")*".itp.jld2"
    prefix = join(skipmissing([name, tstamp]), "_")

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
        @save join([prefix, "BE", postfix], "_") interpolator
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
        @save join(skipmissing(
            [prefix, var, norm ? "normalised" : missing, postfix]
           ), "_"
        ) interpolator
        if deallocateafterwritten
            interpolator = nothing
            GC.gc()
        end
    end
end


function get_itp(
    brxp::BifrostExperiment,
    snaps::Union{Integer,AbstractVector},
    variables::Union{String, Vector{String}},
    itp_type,
    itp_bc,
    outputname=missing,
    deallocateafterwritten=false,
    destagger=true,
    units::String="code",
    slicex::AbstractVector{<:Integer}=Int[],
    slicey::AbstractVector{<:Integer}=Int[],
    slicez::AbstractVector{<:Integer}=Int[],
    precision::DataType=Float32,
)
    # Get spatial axes of the Bifrost experiment and squeeze away dimensions of
    # length 1.
    xaxis, yaxis, zaxis = get_axes(brxp, units=units)
    isempty(slicex) && ( slicex = 1:snapsize[1] )
    isempty(slicey) && ( slicey = 1:snapsize[2] )
    isempty(slicez) && ( slicez = 1:snapsize[3] )
    br_axes = squeeze( (xaxis[slicex], yaxis[slicey], zaxis[slicez]) )
    # Get temporal axis, or the time if only one snapshot.
    times = get_var(brxp, snaps, "t"; units=units)

    # Resolve the the interpolation type and boundary conditions from the
    # arguments. May throw `ArgumentError`.
    itp_type = resolve_itp_type(itp_type)
    n_spatial_dims = length(br_axes)
    static = length(times) == 1
    itp_bc = resolve_itp_bc(itp_bc, n_spatial_dims, static)

    # Put the snapshot(s) of each variable in a vector.
    # If `variable` is a String, not a vector of strings, `vectorfield` will
    # not be a vector, just the snapshot(s) itself.
    fetched_vars = map(variables) do var
        data = get_var(
            brxp, snaps, var;
            units=units,
            destagger=destagger,
            slicex=slicex,
            slicey=slicey,
            slicez=slicez,
            precision=precision,
            squeeze=true # dimensions of length 1 is not supported
        )
        # `data` is a vector if we fetched multiple snapshots. We want time to be
        # a dimension in addition to the spatial ones, hence the stacking.
        if data isa Vector
            return stack(data)
        else
            return data
        end
    end

    # From the fetched variables, axes and itp-type, crate the interpolation
    # object using interpolations.jl
    itp = create_interpolator(fetched_vars, itp_type, br_axes, times)

    # Create a filename for the stored interpolation object
    itp_string = resolve_itp_string(itp_type)
    tstamp = static == 1 ?
             "t$snaps" : "t$(first(snaps))-t$(last(snaps))"
    postfix = join(skipmissing([units, itp_string]), "_")*".itp.jld2"
    prefix = join(skipmissing([outputname, tstamp]), "_")
end

function create_interpolator(fetched_vars, itp_type, br_axes, times)
    n_components = length(fetched_vars)

    # Are we interpolating a vector field in time?
    static = length(times) > 1
    vectorfield = n_components > 1

    # If there is only one component, `fetched_vars` may be a
    # single-element vector. Discard the Vector-wrapping in this case.
    if fetched_vars isa Vector && n_components == 1
        fetched_vars = first(fetched_vars)
    end

    # Add the temporal axis if we are interpolating in time
    axes = static ? br_axes : (br_axes..., times)

    # Use an array of static vectors if we are interpolating a vector field
    field = vectorfield ? 
        to_array_of_svectors(n_components, fetched_vars) :
        fetched_vars

    return extrapolate(
        scale(interpolate(field, itp_type), axes...),
        itp_bc
    )
end

function to_array_of_svectors(::Val{N}, vec) where {N}
    A = first(vec)
    out = similar(A,SVector{N,(eltype(A))})
    @inbounds for i in eachindex(A)
        out[i] = SVector{N}(ntuple(k -> vec[k][i], Val(N)))
    end
    return out
end

function create_itp(
    vectorfield::Matrix,
    br_axes::Tuple{AbstractVector, AbstractVector},
    _::Real,
    itp_func,
    itp_bc
)
    itp_func(br_axes, vectorfield, extrapolation_bc=itp_bc)
end

function create_itp(
    vectorfield::Array{<:Real, 3},
    br_axes::Tuple{AbstractVector, AbstractVector},
    times::AbstractVector,
    itp_func,
    itp_bc
)
    itp_func((br_axes..., times), vectorfield, extrapolation_bc=itp_bc)
end

function create_itp(
    vectorfield::Array{<:Real, 3},
    br_axes::Tuple{AbstractVector, AbstractVector},
    _::Real,
    itp_func,
    itp_bc
)
    itp_func(br_axes, vectorfield, extrapolation_bc=itp_bc)
end

function create_itp(
    vectorfield::Array{<:Real, 4},
    br_axes::Tuple{AbstractVector, AbstractVector},
    times::AbstractVector,
    itp_func,
    itp_bc
)
    itp_func((br_axes..., times), vectorfield, extrapolation_bc=itp_bc)
end

function resolve_itp_type(itp_type)
    if itp_type isa String
        if itp_type == "linear"
            BSpline(Linear())
        elseif itp_type == "cubic-splines"
            BSpline(Cubic())
        else
            unknown_itp_type_error()
        end
    elseif itp_type isa ITP.InterpolationType
        itp_type
    else
        uknown_itp_type_error()
    end
end

function resolve_itp_bc(itp_bc, num_spatial_axes, static)
    if itp_bc isa Vector
        num_dims = num_spatial_axes + static ? 0 : 1
        if length(itp_bc) != num_dims
            throw(ArgumentError(
                "itp_bc is the wrong length. Should be 1 or equal to the"*
                "number of dimensions (including time)"
            ))
        end
    end
    if itp_bc isa Union{String, Vector{String}}
        map(string2itp, itp_bc)
    elseif itp_bc isa Union{
        ITP.BoundaryCondition,
        Vector{ITP.BoundaryCondition}
    }
        itp_bc
    else
        throw(ArgumentError("Uknown itp_bc. Try 'flat' or 'periodic'."))
    end
end

function create_itp(field, itp_type, axes, itp_bc)
end

function string2itpbc(s)
    if s == "flat"
        Flat()
    elseif s == "periodic"
        Periodic()
    else
        unknown_itp_bc_error()
    end
end

function resolve_itp_string(itp_func)
    if itp_func == linear_interpolation
        "linear"
    elseif itp_func == cubic_spline_interpolation
        "cubicspline"
    else
        missing
    end
end

function uknown_itp_type_error()
    throw(ArgumentError("Uknown itp_type. Try 'linear' or 'cubic-splines'."))
end

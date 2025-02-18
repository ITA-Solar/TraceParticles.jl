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
    snap::Integer,
    variable::String,
    ;
    itp_type=BSpline(Linear()),
    itp_bc=Flat(),
    units="si",
    normalise=false,
    kwargs...
)
    if variable == "ne"
        var = get_electron_density(brxp, snap; units=units, kwargs...)
    else
        var = get_var(brxp, snap, variable; units=units, kwargs...)
    end
    var = dropdims(var)
    br_axes = get_axes(brxp, units=units)
    br_axes = dropdims(br_axes)
    #var_itp = linear_interpolation(
    #    br_axes, var, extrapolation_bc=itp_bc
    #    )
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


function get_br_emfield_numdensity_gastemp_interpolator(
    expname::String,
    snap::Integer,
    expdir::String,
    ;
    itp_type=Gridded(Linear()),
    itp_bc=Flat(),
    units=units,
)
    #
    bx = get_var(expname, snap, expdir, "bx"; units="SI", destagger=true)
    by = get_var(expname, snap, expdir, "by"; units="SI", destagger=true)
    bz = get_var(expname, snap, expdir, "bz"; units="SI", destagger=true)
    ex = get_var(expname, snap, expdir, "ex"; units="SI", destagger=true)
    ey = get_var(expname, snap, expdir, "ey"; units="SI", destagger=true)
    ez = get_var(expname, snap, expdir, "ez"; units="SI", destagger=true)
    rho = get_var(expname, snap, expdir, "r"; units="SI")
    tg = get_var(expname, snap, expdir, "tg"; units="SI", destagger=false)
    #
    # Convert mass density to number density
    rho ./= tp.m_p
    fields = eachslice(
        [bx;;;; by;;;; bz;;;; ex;;;; ey;;;; ez;;;; rho;;;; tg],
        dims=(1, 2, 3)
    )
    #
    # Make interpolation-object
    fields = dropdims(fields)
    br_axes = get_axes(brxp, units=units)
    br_axes = dropdims(br_axes)
    fields_itp = interpolate(br_axes, fields, itp_type)
    fields_itp = extrapolate(fields_itp, itp_bc)
    #
    return fields_itp
end


"""
    findslice(x, y, z, xlim, ylim, zlim)
    findslice(x, z, xlim, zlim)
Find the ranges which slices the `BifrostExperiment` according to the given
limits.
"""
function findslice(
    x::AbstractVector,
    y::AbstractVector,
    z::AbstractVector,
    xlim::Any,
    ylim::Any,
    zlim::Any
    ;
    verbose=true
)
    if isnothing(xlim)
        slicex = 1:length(x)
    else
        slicex = findall(x -> xlim[1] <= x <= xlim[2], x)
        if length(slicex) == 0
            error("xlim did not match any points in the mesh. Wrong units?")
        end
    end
    if isnothing(zlim)
        slicez = 1:length(z)
    else
        slicez = findall(z -> zlim[1] <= z <= zlim[2], z)
        if length(slicez) == 0
            error("zlim did not match any points in the mesh. Wrong units?")
        end
    end
    if isnothing(ylim)
        slicey = 1:length(y)
    else
        slicey = findall(y -> y > ylim, y)[1]
        if length(slicez) == 0
            error("zlim did not match any points in the mesh. Wrong units?")
        end
    end
    if verbose
        @info "Range of x-axis: $(slicex[1]:slicex[end])"
        @info "Range of y-axis: $(slicey[1]:slicey[end])"
        @info "Range of z-axis: $(slicez[1]:slicez[end])"
    end
    return slicex, slicey, slicez
end
function findslice(
    x::AbstractVector,
    z::AbstractVector,
    xlim::Any,
    zlim::Any
    ;
    verbose=true
)
    if isnothing(xlim)
        slicex = 1:length(x)
    else
        slicex = findall(x -> xlim[1] <= x <= xlim[2], x)
        if length(slicex) == 0
            error("xlim did not match any points in the mesh. Wrong units?")
        end
    end
    if isnothing(zlim)
        slicez = 1:length(z)
    else
        slicez = findall(z -> zlim[1] <= z <= zlim[2], z)
        if length(slicez) == 0
            error("zlim did not match any points in the mesh. Wrong units?")
        end
    end
    if verbose
        @info "Range of x-axis: $(slicex[1]:slicex[end])"
        @info "Range of z-axis: $(slicez[1]:slicez[end])"
    end
    return slicex, slicez
end

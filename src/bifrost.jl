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
    var = dropdims(var)
    br_axes = get_axes(brxp, units=units)
    br_axes = dropdims(br_axes)
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
    axes = dropdims(axes)
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

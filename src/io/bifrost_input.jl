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
        expname ::String,
        snap    ::Integer,
        expdir  ::String,
        variable::String,
        ;
        itp_type=BSpline(Linear()),
        itp_bc=Flat(),
        kwargs...
    )
    var = get_var(expname, snap, expdir, variable; kwargs...)
    var = dropdims(var)
    br_axes = load_br_axes(expname, snap, expdir)
    br_axes = dropdims(br_axes)
        #var_itp = linear_interpolation(
    #    br_axes, var, extrapolation_bc=itp_bc
    #    )
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
        mask = [isassigned(range_axes,i) for i in eachindex(range_axes)]
        if all(mask)
            br_axes = Tuple(range_axes)
        else
            @warn "Non-uniform grid. Creating uniform axes " *
                "using range."
            idx = findall(x->!x, mask)
            for i in idx
                range_axes[i] = range(start=br_axes[i][1], stop=br_axes[i][end], length=length(br_axes[i]))
            end
            br_axes = Tuple(range_axes)
        end
        itp_func = cubic_spline_interpolation
        #var_itp = interpolate(var, itp_type)
        #var_itp = extrapolate(scale(var_itp, br_axes...), itp_bc)
    end

    var_itp = itp_func(br_axes, var, extrapolation_bc = itp_bc)
    return var_itp, br_axes
end


function get_br_emfield_vecof_interpolators(args...; kwargs...)
    vars = ("bx", "by", "bz", "ex", "ey", "ez")
    return [get_br_var_interpolator(args..., var; kwargs...)[1] for var in vars]
end


function get_br_gcafields_vecof_interpolators(
        expname ::String,
        snap    ::Integer,
        expdir  ::String,
        ;
        itp_bc=Flat(),
        itp_type=Gridded(Linear()),
        destagger=true,
        destaggergradients=true,
        diffscheme=bifrostscheme,
        kwargs...
        )
    vars = ("bx", "by", "bz", "ex", "ey", "ez")
    emfields = [get_var(expname, snap, expdir, var;
        destagger=destagger, kwargs...)
        for var in vars] # Vector{Array{3, <:Real}}
    B_vec = stack(emfields[1:3]) # Array{4, <:Real}
    E_vec = stack(emfields[4:6]) # Array{4, <:Real}

    println("Avg E: $(sum(E_vec)/prod(size(E_vec)))")
    println("Scaling:")
    #B_vec *= 1.121e3*1e-4 # u_B*cgs2SI_b
    #E *= 1e6*1.121e3*1e-6 # u_u*u_B*cgs2SI_E
    println("Avg E after scaling: $(sum(E_vec)/prod(size(E_vec)))")


    # br_axes = Tuple{Vector{<:Real}, Vector{<:Real}, Vector{<:Real}}
    br_axes = load_br_axes(expname, snap, expdir)

    println("Computing gradients...")
    # ∇B = Vector{Array{3, <:Real}}(3)
    # ∇b = ∇ExB = Vector{Array{3, <:Real}}(9)
    ∇B, ∇b̂, ∇ExBdrift = tp.gcagradients_from_emfield(
        B_vec,
        E_vec,
        br_axes...,
        diffscheme,
        )

    # Drop extra dimensions
    B_vec = dropdims(B_vec) # Array{2, <:Real}
    E_vec = dropdims(E_vec)         # Array{2, <:Real}
    ∇B = dropdims.(∇B)      # Vector{Array{2, <:Real}}(3)
    ∇b̂ = dropdims.(∇b̂)      # Vector{Array{2, <:Real}}(9)
    ∇ExBdrift = dropdims.(∇ExBdrift) # Vector{Array{2, <:Real}}(9)
    br_axes = dropdims(br_axes) # Tuple{Vector{<:Real}, Vector{<:Real}

    if itp_type == BSpline(Linear())
        itp_func = linear_interpolation
    elseif itp_type == BSpline(Cubic())
        range_axes = Vector{StepRangeLen}(undef, length(br_axes))
        for i in eachindex(br_axes)
            gridspacing = diff(br_axes[i])
            if all(gridspacing[1] .== gridspacing)
                range_axes[i] = range(br_axes[i][1], br_axes[i][end], length(br_axes[i]))
            end
        end
        mask = [isassigned(range_axes,i) for i in eachindex(range_axes)]
        if all(mask)
            br_axes = Tuple(range_axes)
        else
            @warn "Non-uniform grid. Creating uniform axes " *
                "using range."
            idx = findall(x->!x, mask)
            for i in idx
                range_axes[i] = range(start=br_axes[i][1], stop=br_axes[i][end], length=length(br_axes[i]))
            end
            br_axes = Tuple(range_axes)
        end
        itp_func = cubic_spline_interpolation
    elseif itp_type == Gridded(Linear())
        nothing
    else
        error("Unsupported interpolation type.")
    end

    itps = Vector{AbstractInterpolation}(undef, 27)
    if destaggergradients
        # E⃗ and B⃗ are cell centred. The gradients have been calculated with
        # an upwind scheme, so they are on the cell faces, needing different
        # interpolation axes than E⃗ and B⃗.
        if length(br_axes) == 2
            # Ill defined upper boundary
            xcentre_axis = br_axes[1]
            zcentre_axis = br_axes[2]
            if typeof(br_axes[1]) <: StepRangeLen
                dx = step(br_axes[1])
                dz = step(br_axes[2])
                xface_axis = br_axes[1][1] + dx/2:dx:br_axes[1][end] + dx/2
                zface_axis = br_axes[2][1] + dz/2:dz:br_axes[2][end] + dz/2
            else
                xface_axis = br_axes[1] .+ [diff(xcentre_axis)/2; 0.0]
                zface_axis = br_axes[2] .+ [diff(zcentre_axis)/2; 0.0]
            end
        else
            error("Cell centred gradient interpolation objects unavailable in 3D")
        end
        # Create interpolation objects of cell centred components.
        itps[1:6] = [itp_func(br_axes, var, extrapolation_bc=itp_bc)
            for var in ([B_vec[:,:,i] for i in 1:3]...,[E_vec[:,:,i] for i in 1:3]...)]
        # Create interpolation objects of upper x-face components, i.e ∂ₓfᵢ components.
        itps[[7, 10:12..., 19:21...]] = [itp_func(
            (xface_axis, zcentre_axis), var, extrapolation_bc=itp_bc)
            for var in (∇B[1], ∇b̂[1:3]...,∇ExBdrift[1:3]...)]
        # Create interpolation objects of upper z-face components, i.e ∂_zf_z components.
        itps[[9, 16:18..., 25:27...]] = [itp_func(
                (xcentre_axis, zface_axis), var, extrapolation_bc=itp_bc)
                for var in (∇B[3], ∇b̂[7:9]...,∇ExBdrift[7:9]...)]
        # zero itp struct for the ∂_zf_z components...
        # A Functor returning 0.0 that is a subtype og AbstractInterpolation is 3
        # times faster. A function returning 0.0 is 30 times faster. However, then
        # itps have to redefined, and maybe become type unstable, which may slow
        # down the code.
        itps[[8, 13:15..., 22:24...]] = [itp_func(
            (0:1, 0:1), zeros(2,2), extrapolation_bc=itp_bc) for i in 1:7]
        #
    else
        itps[:] = [itp_func(br_axes, var, extrapolation_bc=itp_bc)
            for var in (
                [B_vec[:,:,i] for i in 1:3]...,
                [E_vec[:,:,i] for i in 1:3]...,
                ∇B...,
                ∇b̂...,
                ∇ExBdrift...,
                )
            ]
    end
    return itps
end



function get_br_emfield_interpolator(
        expname::String,
        snap   ::Integer,
        expdir ::String,
        ;
        itp_type=Gridded(Linear()),
        itp_bc=Flat(),
        )
    bx = get_var(expname, snap, expdir, "bx"; units="SI", destagger=true)
    by = get_var(expname, snap, expdir, "by"; units="SI", destagger=true)
    bz = get_var(expname, snap, expdir, "bz"; units="SI", destagger=true)
    ex = get_var(expname, snap, expdir, "ex"; units="SI", destagger=true)
    ey = get_var(expname, snap, expdir, "ey"; units="SI", destagger=true)
    ez = get_var(expname, snap, expdir, "ez"; units="SI", destagger=true)
    #
    emfield = eachslice(stack([bx, by, bz, ex, ey, ez]), dims=(1,2,3))
   # nx,ny,nz = size(bx)
   # emfield = [[bx[i,j,k],by[i,j,k],bz[i,j,k],ex[i,j,k],ey[i,j,k],ez[i,j,k]] for i=1:nx, j = 1:ny, k = 1:nz]
    #
    # Make interpolation-object
    emfield = dropdims(emfield)
    br_axes = load_br_axes(expname, snap, expdir)
    br_axes = dropdims(br_axes)
    emfields_itp = linear_interpolation(
        br_axes, emfield, extrapolation_bc=itp_bc
        )
    #emfields_itp = interpolate(br_axes, emfield, itp_type)
    #emfields_itp = extrapolate(emfields_itp, itp_bc)
    #
    return emfields_itp
end


function get_br_emfield_numdensity_gastemp_interpolator(
        expname::String,
        snap   ::Integer,
        expdir ::String,
        ;
        itp_type=Gridded(Linear()),
        itp_bc=Flat(),
        )
    #
    bx = get_var(expname, snap, expdir, "bx"; units="SI", destagger=true)
    by = get_var(expname, snap, expdir, "by"; units="SI", destagger=true)
    bz = get_var(expname, snap, expdir, "bz"; units="SI", destagger=true)
    ex = get_var(expname, snap, expdir, "ex"; units="SI", destagger=true)
    ey = get_var(expname, snap, expdir, "ey"; units="SI", destagger=true)
    ez = get_var(expname, snap, expdir, "ez"; units="SI", destagger=true)
    rho= get_var(expname, snap, expdir, "r"; units="SI")
    tg = get_var(expname, snap, expdir, "tg"; units="SI", destagger=false)
    #
    # Convert mass density to number density
    rho ./= tp.m_p
    fields = eachslice(
        [bx;;;; by;;;; bz;;;; ex;;;; ey;;;; ez;;;; rho;;;; tg],
        dims=(1,2,3)
        )       
    #
    # Make interpolation-object
    fields = dropdims(fields)
    br_axes = load_br_axes(expname, snap, expdir)
    br_axes = dropdims(br_axes)
    fields_itp = interpolate(br_axes, fields, itp_type)
    fields_itp = extrapolate(fields_itp, itp_bc)
    #
    return fields_itp
end


function load_br_axes(
        expname::String,
        snap   ::Integer,
        expdir ::String,
        )
    basename = string(
        expdir, "/", expname
        )
    idl_filename = string(basename, "_", snap, ".idl")
    mesh_filename = string(basename, ".mesh")
    mesh = BifrostMesh(mesh_filename)
    params = read_params(idl_filename)
    # Scale br_axes
    code2cgs_l = parse(Float32, params["u_l"]) # exploit promotion
    x = code2cgs_l * tp.cgs2SI_l * mesh.x
    y = code2cgs_l * tp.cgs2SI_l * mesh.y
    z = code2cgs_l * tp.cgs2SI_l * mesh.z
    return (x, y, z)
end

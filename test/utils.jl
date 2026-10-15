"""
    create3Daxes(
        (x0, y0, z0)::Tuple{Real,Real,Real},
        (xf, yf, zf)::Tuple{Real,Real,Real},
        (nx, ny, nz)::Tuple{Integer,Integer,Integer}
    )
Create axes from axis limits and number of grid points.
"""
function create3Daxes(
    (x0, y0, z0)::Tuple{Real,Real,Real},
    (xf, yf, zf)::Tuple{Real,Real,Real},
    (nx, ny, nz)::Tuple{Integer,Integer,Integer}
)
    xx = range(x0, xf, length=nx)
    dx = xx[2] - xx[1]
    yy = range(y0, yf, length=ny)
    dy = yy[2] - yy[1]
    # Account for single point in the z-axis
    if nz == 1
        zz = [z0]
        dz = 0.0
    else
        zz = range(z0, zf, length=nz)
        dz = zz[2] - zz[1]
    end
    return xx, yy, zz, dx, dy, dz
end # function createaxes


function discretise3D(
    func::Function,
    xx::AbstractVector,
    yy::AbstractVector,
    zz::AbstractVector,
)
    f111 = func(xx[1], yy[1], zz[1])
    field = Array{typeof(f111)}(
        undef, length(xx), length(yy), length(zz)
    )
    for i in eachindex(xx)
        for j in eachindex(yy)
            for k in eachindex(zz)
                field[i, j, k] = func(xx[i], yy[j], zz[k])
            end
        end
    end
    return field
end

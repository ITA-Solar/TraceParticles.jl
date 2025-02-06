"""
    magneticmirrorfield(
        x::Real,
        y::Real,
        z::Real,
        B0::Real,
        L ::Real
        )
Return the magnetic field at (`x`, `y`, `z`) in a magnetic mirror/bottle with 
length `L` and strength `B0`.

B⃗(x, y, z) = -xzB₀/L² x̂  - yzB₀/L² ŷ  + (B₀ + zB₀/L²)ẑ
"""
function magneticmirrorfield(
    x::Real,
    y::Real,
    z::Real,
    B0::Real,
    L::Real,
)
    a = B0 * z / L^2
    return [-x * a, -y * a, B0 + z * a]
end # mirroringfield

"""
    magneticdipolefield(
        x::Real,
        y::Real,
        z::Real,
        M ::Real
        )
Return the magnetic field at (`x`, `y`, `z`) in a magnetic dipole field with 
mangetic moment `M`.

B⃗(x, y, z) = 3αzx x̂  + 3αzy ŷ  + α(2z² - x² - y²)ẑ, 
where 
α = M/(x² + y² + z²)^(2.5).
"""
function magneticdipolefield(
    x::Real,
    y::Real,
    z::Real,
    M::Real,
)
    a = M / (x^2 + y^2 + z^2)^(5 / 2)
    return [3a * z * x, 3a * z * y, a * (2z^2 - x^2 - y^2)]
end


"""
    magneticfield_2Dplasmoid(
        x::Real,
        y::Real,
        z::Real,
        μ_x::Real,
        μ_y::Real,
        σ_x::Real,
        σ_y::Real,
        normfactor::Real,
        B0::Vector{<:Real},
        )
Return the magnetic field at (`x`, `y`, `z`) in a plasmoid configuration with
background field `B0`. The field is analytical with a vector potential

A⃗(x,y,z) = normfactor * N_x(μ_x, σ_x)*N_y(μ_y, σ_y) * ẑ - B0_x * zŷ + B0_y * zx̂

where N_x is normal distribution of x around mean `μ_x` with standard deviation
`σ_x`.
"""
function magneticfield_2Dplasmoid(
    x::Real,
    y::Real,
    z::Real,
    μ_x::Real,
    μ_y::Real,
    σ_x::Real,
    σ_y::Real,
    normfactor::Real,
    B0::Vector{<:Real},
)
    A_z = normfactor / (σ_x * σ_y * 2π) * exp(-0.5 * (((x - μ_x) / σ_x)^2 + ((y - μ_y) / σ_y)^2))
    return [-A_z * (y - μ_y) / σ_y^2 + B0[1], A_z * (x - μ_x) / σ_x^2 + B0[2], B0[3]]
end


"""
    magneticfieldstrength_2Dplasmoid(args...)
Return the magnetic field strength of the 2D plasmoid configuration defined in
@magneticfield_2Dplasmoid.
"""
function magneticfieldstrength_2Dplasmoid(args...)
    B = magneticfield_2Dplasmoid(args...)
    return norm(B)
end

"""
    magneticfieldstrengthgradient_2Dplasmoid(
        x::Real,
        y::Real,
        z::Real,
        μ_x::Real,
        μ_y::Real,
        σ_x::Real,
        σ_y::Real,
        normfactor::Real,
        B0::Vector{<:Real},
        )
Return the gradient magnetic field strength of the 2D plasmoid configuration
defined in @magneticfield_2Dplasmoid.
"""
function magneticfieldstrengthgradient_2Dplasmoid(
    x::Real,
    y::Real,
    z::Real,
    μ_x::Real,
    μ_y::Real,
    σ_x::Real,
    σ_y::Real,
    normfactor::Real,
    B0::Vector{<:Real},
)
    frac_x = (x - μ_x) / σ_x
    frac_y = (y - μ_y) / σ_y
    A_z = normfactor / (2π * σ_x * σ_y) * exp(-0.5 * (frac_x^2 + frac_y^2))
    frac_3 = 2A_z * frac_x * frac_y / (σ_x * σ_y)
    α = -A_z * frac_y / σ_y + B0[1]
    β = A_z * frac_x / σ_x + B0[2]
    γ = sqrt(α^2 + β^2 + B0[3]^2)
    dβ2dx = 2β * A_z / σ_x^2 * (1 - frac_x^2)
    dβ2dy = -β * frac_3
    dα2dx = α * frac_3
    dα2dy = -2α * A_z / σ_y^2 * (1 - frac_y^2)

    ∇B_x = dα2dx + dβ2dx
    ∇B_y = dα2dy + dβ2dy
    #∇B_x = (x - μ_x)/σ_x^2 * (1/(γ*σ_x^2) - γ)
    #∇B_y = (y - μ_y)/σ_y^2 * (1/(γ*σ_y^2) - γ)
    return 1 / (2γ) * [∇B_x, ∇B_y, 0]
end



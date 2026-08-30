using OrdinaryDiffEq: ODEProblem, solve
using LinearAlgebra: norm
using TraceParticles
using CairoMakie

#-------------------------------------------------------------------------------
# Electromagnetic field
function magneticdipole(x, y, z; M)
    a = M / (x^2 + y^2 + z^2)^(5 / 2)
    return [3a * z * x, 3a * z * y, a * (2z^2 - x^2 - y^2)]
end

function T_dipole(; q, m, M, R0, v0, α)
    return @. 2π * q * M / (m * v0^2 * R0) * (1 - 1 / 3 * sin(α)^0.62)
end

#-------------------------------------------------------------------------------
# Parameters
tf = 80.0
tspan = (0.0, tf)
mass = 1.0
charge = 1.0
# Initial position is such that vperp is 1 and the guiding centre is
# at x⃗ = [1, 0, 0]. It is also assumed that the GCA is valid such that B is equal
# to B(x⃗) = 2M, directed in -ẑ. We may then compare with the Walt (1994)
# approximation of the azimuthal drift-period.
R0 = [1.0, 0, 0]
vparal = 0.5
vperp = 1.0
# Magnetic dipole strength
qMm = 5:5:60
dipolestrength = qMm * mass / charge

"""
    dipolesimulation(dipolestrength)
Given a `dipolestrength`, solve the Lorentz equation and the guiding centre
approximation. Also return the analytical drift
"""
function dipolesimulation(dipolestrength)
    # Dipole-specific field and initial conditions
    emfield(x, y, z, t) = zeros(3), magneticdipole(x, y, z; M=dipolestrength)
    E, B = emfield(R0..., 0)
    μ = TraceParticles.magneticmoment(vperp, mass, norm(B))
    rL = TraceParticles.larmorradius(mass, vperp, charge, norm(B))
    u0 = get_fullorbit(B, E, R0, vparal, μ, charge, mass, pi / 2)

    # Create problem
    fo_prob = ODEProblem(
        lorentzforce!,
        u0,
        tspan,
        (charge=charge, mass=mass, electromagneticfield=emfield)
    )
    gca_prob = ODEProblem(
        guidingcentreapproximation!,
        [R0...; vparal],
        tspan,
        (
            charge=charge,
            mass=mass,
            magneticmoment=μ,
            electromagneticfield=emfield
        )
    )

    # Run simulation
    fo_sim = solve(fo_prob; reltol=1e-7, abstol=1e-9)
    gca_sim = solve(gca_prob)


    # Analytical result
    v0 = norm(u0[4:6])
    α = atan(vperp / vparal)
    R0x = R0[1]
    T = T_dipole(; q=charge, m=mass, M=dipolestrength, R0=R0x, v0=v0, α=α)
    angularfreq = 2π / T
    phi = 2π * tf / T

    # Return positions and analytical drift
    times0 = range(0.0, tf, length=100)
    times1 = range(0.0, tf, length=1_000)
    times2 = range(0.0, tf, length=10_000)
    return (
        fo_sim(times2),
        gca_sim(times1),
        R0x * cos.(angularfreq * times0),
        R0x * sin.(angularfreq * times0),
        phi,
        rL
    )
end

pos_fo, pos_gca, driftx, drifty, _, _ = dipolesimulation(40 * mass / charge)

# Plot
fig = Figure()
ax = Axis3(fig[1:2,1:2]; aspect=:data)
lines!(ax, pos_fo[1, :], pos_fo[2, :], pos_fo[3, :]; label="Full orbit")
lines!(
    ax, pos_gca[1, :], pos_gca[2, :], pos_gca[3, :];
    label="Full orbit", linewidth=1
)
lines!(
    ax, driftx, drifty, [0.0 for _ in driftx];
    label="Analytical drift", linewidth=1
)
#axislegend(ax; position=:lt)
Legend(fig[3,2], ax)
ax.xlabel = "x"
ax.ylabel = "y"

ax2 = Axis(fig[3,1]; aspect=DataAspect())
#lines!(ax2, pos_fo[1, :], pos_fo[2, :])
lines!(ax2, pos_gca[1, :], pos_gca[2, :]; color=Makie.wong_colors()[2])
lines!(ax2, driftx, drifty; color=Makie.wong_colors()[3])
ax2.xlabel = "x"
ax2.ylabel = "y"
fig

N = length(dipolestrength)
phis_fo = Vector{Float64}(undef, N)
phis_gca = Vector{Float64}(undef, N)
phis_analytical = Vector{Float64}(undef, N)
rLs = Vector{Float64}(undef, N)
for i in eachindex(dipolestrength)
    global pos_fo
    global pos_gca
    global driftx
    global drifty
    M = dipolestrength[i]
    pos_fo, pos_gca, driftx, drifty, phi, rL = dipolesimulation(M)
    # The calculation of the angle phi is by evaluating atan(y/x).
    # Rotations more than π/2 degrees needs to be adjusted because atan(x/y)
    # jumps to -π/2 in the second circle quadrant.
    # Angles in (0, π/2) needs no addition
    # Angles in (π/2, 3π/2) needs +π 
    # Angles in (3/2, 5π/2) needs +2π
    # etc...
    # I use the analytical angle as reference for how much to add to the
    # numerical angles.
    n = div(phi, pi/2)
    extra_rads = (div(n, 2) + n % 2) * pi
    phis_fo[i] = atan(pos_fo[2, end] / pos_fo[1, end]) + extra_rads
    phis_gca[i] = atan(pos_gca[2, end] / pos_gca[1, end]) + extra_rads
    phis_analytical[i] = phi
    rLs[i] = rL
end


fig2 = Figure()
ax3 = Axis(fig2[1,1])
#scatterlines!(ax3, dipolestrength, phis_fo, label="Full orbit")
scatterlines!(ax3, dipolestrength, rLs, label="Larmor radius")
scatterlines!(ax3, dipolestrength, abs.(phis_gca .- phis_fo);
    label="GCA", linestyle=:dash, marker=:diamond
)
scatterlines!(ax3, dipolestrength, abs.(phis_analytical .- phis_fo);
    label="Analytical drift", linestyle=:dot, marker=:utriangle
)
#ax3.ylabel="Final azimuthal angle"
ax3.ylabel="Distance [m]"
ax3.xlabel="Magnitisation of particle, qM/m"
ax3.title = "Difference in azimuthal drift compared to fullorbit"
ax3.yscale=log10
ax3.xticks = (qMm)
ax3.yticks = (10.0 .^ (-3:0))
ax3.yminorticks=IntervalsBetween(10)
ax3.yminorticksvisible=true
ax3.yminorgridvisible=true
axislegend(ax3)
fig2


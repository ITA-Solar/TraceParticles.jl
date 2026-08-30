```@meta
CurrentModule = TraceParticles
```

# Gallery

## [Magnetic mirroring](@id mirror-verifaction)
When a charge particle experiences increased magnetic field strength, the gyration radius decreases.
Field strength increases in the plane of the gyration give rise to a drift along $\nabla B/\times\mathbf{B}$.
If the field strength increases perpendicular to the gyration, parallel to the magnetic field direction, the paricle might be mirrored if the magnetic moment of the particle is conserved.
As the the Larmor radius decreases, conservation of the moment requires a larger perpendicular velocity, and if the total energy is conserved (for example if there are no electric field), the parallel velocity must decrease.
If the magnetic field becomes strong enough, all the energy will be transfered from the parallel component into the perpendicular, and the particle will be mirrored.
In other words, in a field with a magnet gradient that i along its direction, a charged particle will experience a force in the opposite direction, $F_\text{mirror} \propto \nabla B \cdot \hat{\mathbf{b}}$.
```@example mirror-verification
using OrdinaryDiffEq
using LinearAlgebra
using TraceParticles
using CairoMakie

#...............................................................................
# Electromagnetic field
function magneticbottle(x, y, z; B0, L)
    a = B0 * z / L^2
    return [-x * a, -y * a, B0 + z * a]
end

B0 = 10.0
L = 0.1
emfield(x, y, z, t) = zeros(3), magneticbottle(x, y, z; B0=B0, L=L)

#...............................................................................
# Particle parameters
tf = 1.5
tspan = (0, tf)
mass = 1
charge = 10

#...............................................................................
# Initial conditions
vel0 = [0.0, 0.4, 0.4]
rL = mass * √(vel0[1]^2 + vel0[2]^2) / (charge * B0)
pos0 = [-1rL, 0.0, 0.0]
E, B = emfield(pos0..., 0)
R0, vparal, μ = get_guidingcentre(pos0, vel0, B, E, charge, mass)

#...............................................................................
# Create problem
prob_FO = ODEProblem(
    lorentzforce!,
    [pos0; vel0],
    tspan,
    (charge=charge, mass=mass, electromagneticfield=emfield),
)
prob_GCA = ODEProblem(
    guidingcentreapproximation!,
    [R0; vparal],
    tspan,
    (
        charge=charge,
        mass=mass,
        magneticmoment=μ,
        electromagneticfield=emfield
    )
)

#...............................................................................
# Run simulation
sol_FO = solve(prob_FO)
sol_GCA = solve(prob_GCA)

#-------------------------------------------------------------------------------
# Analytical solution

function z_analytical(t; μ, m, B0, L, A, ϕ)
    ω = √(2μ * B0 / (m* L^2))
    return @. A * sin(ω * t + ϕ)
end
Bmax = B0 * norm(vel0)^2 / (vel0[1]^2 + vel0[2]^2)
zmax = L * √(Bmax / B0 - 1)
ϕ = 0
times = range(0.0, tf, length=1000)
analytical = z_analytical(times; μ=μ, m=mass, B0=B0, L=L, A=zmax, ϕ=ϕ)

#-------------------------------------------------------------------------------
# Plot
fig = Figure()
ax = Axis(fig[1,1])
lines!(ax, [sol_FO(t)[3] for t in times], [sol_FO(t)[1] for t in times],
    label="Full orbit", linewidth=0.5
)
lines!(ax, sol_GCA[3,:], sol_GCA[1,:];
    label="GCA", linewidth=2.0
)
lines!(ax, analytical, [0.0 for _ in analytical];
    linestyle=:dash, label="Analytical solution", linewidth=1.0
)
ax.aspect = DataAspect()
ax.limits = ((nothing, nothing), (-0.01,0.03))
ax.yticks = ([-0.01,0.00,0.01])
axislegend(ax, position=:ct, orientation=:horizontal, framevisible=false)
```

## [Magnetic dipole](@id dipole-verification)
Charged particles in magnetic dipole can be trapped mirroring from pole to pole while drifting azimuthally due to the gradients and curvature of the magnetic field.
In this example we compare the results of a direct solution of the [`lorentzforce!`](@ref eom-api) with the solution of the `guidingcentreapproximation!` and an analytical approximation to the drift.
```@example dipole-verification
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
M = 40 * mass / charge

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

pos_fo, pos_gca, driftx, drifty, _, _ = dipolesimulation(M)

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
```

```@example dipole-verification
# Varying magnetisation
qMm = 10:5:60
dipolestrength = qMm * mass / charge
# (... see source file all code ...)
```
```@setup dipole-verification
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
    global M
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
scatterlines!(ax3, dipolestrength, rLs, label="Larmor radius")
scatterlines!(ax3, dipolestrength, abs.(phis_gca .- phis_fo);
    label="GCA", linestyle=:dash, marker=:diamond
)
scatterlines!(ax3, dipolestrength, abs.(phis_analytical .- phis_fo);
    label="Analytical drift", linestyle=:dot, marker=:utriangle
)
ax3.ylabel="Distance [m]" # hide
ax3.xlabel="Magnitisation of particle, qM/m" # hide
ax3.title = "Difference in azimuthal drift compared to fullorbit"
ax3.yscale=log10
ax3.xticks = (qMm)
ax3.yticks = (10.0 .^ (-3:0))
ax3.yminorticks=IntervalsBetween(10)
ax3.yminorticksvisible=true
ax3.yminorgridvisible=true
axislegend(ax3)
```
```@example dipole-verification
fig2
```

## [Dynamics at a current sheet (Speiser orbit)](@id speiser-verification)

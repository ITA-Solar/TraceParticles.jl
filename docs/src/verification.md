```@meta
CurrentModule = TraceParticles
```

# Verification gallery

Every figure on this page is generated when the documentation is built, by
running `TraceParticles` on a problem whose answer is known in closed form and
plotting the numerical result on top of the analytical one. Nothing here is
pre-rendered, so if a solver regresses the plots stop agreeing and you can see it
at a glance.

All of these use analytically prescribed fields, so each example is
self-contained and can be pasted into the REPL.

```@setup gallery
using TraceParticles
# Physics helpers, callback affects and samplers that the package defines but
# does not export, since the high-level interface reaches them for you.
using TraceParticles: get_guidingcentre, get_fullorbit, kineticenergy,
    larmorradius, magneticmoment, J2eV, hybridswitchaffect!,
    maxwellianvelocitysample, rejectionsample
using OrdinaryDiffEq
using DiffEqCallbacks
using StaticArrays
using LinearAlgebra
using Statistics
using Random
using Printf
using CairoMakie
CairoMakie.activate!(type="png")
```

## E×B drift

A uniform magnetic field `B = B₀ẑ` with a perpendicular electric field
`E = E₀ŷ`. The guiding centre drifts at

```math
\mathbf{v}_E = \frac{\mathbf{E}\times\mathbf{B}}{B^2},
```

independently of charge, mass and speed. The full-orbit particle gyrates about
that drifting centre; the guiding-centre solver should follow the straight drift
line exactly.

```@example gallery
q = -TraceParticles.e     # electron
m = TraceParticles.m_e
B0 = 3.37213e-4           # T
vperp0 = 1.0e5            # m/s
v_E = 0.25vperp0          # choose E so the drift is a quarter of v⊥
E0 = v_E * B0

emfield(x, y, z, t) = SVector(0.0, E0, 0.0), SVector(0.0, 0.0, B0)

gyroperiod = 2π * m / abs(q * B0)
tspan = (0.0, 8gyroperiod)
u0 = [0.0, 0.0, 0.0, vperp0, 0.0, 0.0]

sol_fo = solve(
    ODEProblem(lorentzforce!, u0, tspan,
        FullOrbitParams(charge=q, mass=m, electromagneticfield=emfield)),
    Tsit5(); reltol=1e-10, abstol=1e-14,
)

E, B = emfield(0.0, 0.0, 0.0, 0.0)
R0, vparal0, μ = get_guidingcentre(
    SVector(u0[1], u0[2], u0[3]), SVector(u0[4], u0[5], u0[6]), B, E, q, m)

sol_gca = solve(
    ODEProblem(guidingcentreapproximation!, [R0..., vparal0], tspan,
        GCAParams(charge=q, mass=m, electromagneticfield=emfield, magneticmoment=μ)),
    Tsit5(); reltol=1e-10, abstol=1e-14,
)
nothing # hide
```

The dashed line is `R₀ + v_E t`. The right-hand panel reconstructs the guiding
centre from the full-orbit solution with [`get_guidingcentre`](@ref) and compares
it with the guiding-centre solver directly — this checks the coordinate
transformation as well as the drift.

```@example gallery
t = range(tspan..., length=2000)
mm = 1e3   # metres -> millimetres

# Guiding centre reconstructed from the full-orbit trajectory
gc_from_fo = map(t) do tᵢ
    u = sol_fo(tᵢ)
    Eᵢ, Bᵢ = emfield(u[1], u[2], u[3], tᵢ)
    R, _, _ = get_guidingcentre(
        SVector(u[1], u[2], u[3]), SVector(u[4], u[5], u[6]), Bᵢ, Eᵢ, q, m)
    R
end

fig = Figure(size=(900, 360))
ax = Axis(fig[1, 1], xlabel="x [mm]", ylabel="y [mm]", title="Trajectories")
lines!(ax, mm .* [sol_fo(tᵢ)[1] for tᵢ in t], mm .* [sol_fo(tᵢ)[2] for tᵢ in t],
    label="full orbit")
lines!(ax, mm .* [sol_gca(tᵢ)[1] for tᵢ in t], mm .* [sol_gca(tᵢ)[2] for tᵢ in t],
    linewidth=3, label="GCA")
lines!(ax, mm .* (v_E .* t), mm .* fill(R0[2], length(t)),
    linestyle=:dash, color=:black, label="R₀ + v_E t")
axislegend(ax, position=:lt)

ax2 = Axis(fig[1, 2], xlabel="t / T_c", ylabel="x [mm]",
    title="Guiding centre from the full orbit")
lines!(ax2, t ./ gyroperiod, mm .* [g[1] for g in gc_from_fo],
    label="get_guidingcentre(full orbit)")
lines!(ax2, t ./ gyroperiod, mm .* [sol_gca(tᵢ)[1] for tᵢ in t],
    linestyle=:dash, linewidth=3, label="GCA solver")
axislegend(ax2, position=:lt)
fig
```

```@example gallery
err = maximum(norm(sol_gca(tᵢ)[1:3] .- (R0 .+ SVector(v_E * tᵢ, 0.0, 0.0))) for tᵢ in t)
@printf("Larmor radius            : %.3e m\n", larmorradius(m, vperp0, q, B0))
@printf("max |R_GCA − R_analytic| : %.3e m\n", err)
```

## ∇B drift

A straight magnetic field whose strength increases linearly across the field
lines, `B = (B₀ + a y) ẑ`. The guiding centre drifts perpendicular to both `B`
and `∇B` at

```math
\mathbf{v}_{\nabla B} = \frac{\mu}{qB}\,\hat{\mathbf{b}}\times\nabla B
= -\frac{m v_\perp^2 a}{2qB^2}\,\hat{\mathbf{x}}.
```

```@example gallery
m = 1.0; q = 1.0
B0 = 6.0; a = 10.0
gradbfield(x, y, z, t) = SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, B0 + a * y)

pos0 = SVector(0.5, 0.5, 0.5)
vel0 = SVector(0.0, 0.1, 0.0)
tspan = (0.0, 10 / 1.7507)   # ten gyrations

E, B = gradbfield(pos0..., 0.0)
R0, vparal0, μ = get_guidingcentre(pos0, vel0, B, E, q, m)

sol_fo = solve(
    ODEProblem(lorentzforce!, [pos0..., vel0...], tspan,
        FullOrbitParams(charge=q, mass=m, electromagneticfield=gradbfield)),
    Tsit5(); reltol=1e-10, abstol=1e-12)
sol_gca = solve(
    ODEProblem(guidingcentreapproximation!, [R0..., vparal0], tspan,
        GCAParams(charge=q, mass=m, electromagneticfield=gradbfield, magneticmoment=μ)),
    Tsit5(); reltol=1e-10, abstol=1e-12)

v_gradB = -μ * a / (q * norm(B))    # analytical drift speed, along -x

t = range(tspan..., length=2000)
fig = Figure(size=(900, 360))
ax = Axis(fig[1, 1], xlabel="x", ylabel="y", title="∇B drift")
lines!(ax, [sol_fo(tᵢ)[1] for tᵢ in t], [sol_fo(tᵢ)[2] for tᵢ in t], label="full orbit")
lines!(ax, [sol_gca(tᵢ)[1] for tᵢ in t], [sol_gca(tᵢ)[2] for tᵢ in t],
    linewidth=3, label="GCA")
lines!(ax, R0[1] .+ v_gradB .* t, fill(R0[2], length(t)),
    linestyle=:dash, color=:black, label="analytic")
axislegend(ax)

ax2 = Axis(fig[1, 2], xlabel="t", ylabel="x")
lines!(ax2, t, [sol_gca(tᵢ)[1] for tᵢ in t], label="GCA")
lines!(ax2, t, R0[1] .+ v_gradB .* t, linestyle=:dash, color=:black, label="analytic")
axislegend(ax2)
fig
```

```@example gallery
@printf("analytic v_∇B  : %.8e\n", v_gradB)
@printf("final x, GCA   : %.8f\n", sol_gca(last(tspan))[1])
@printf("final x, exact : %.8f\n", R0[1] + v_gradB * last(tspan))
```

## Magnetic mirror: bounce motion and the adiabatic invariant

The magnetic bottle of Ripperda et al. (2018),

```math
\mathbf{B}(x,y,z) = -\frac{xzB_0}{L^2}\hat{\mathbf{x}}
                    -\frac{yzB_0}{L^2}\hat{\mathbf{y}}
                    + \Big(B_0 + \frac{z^2B_0}{L^2}\Big)\hat{\mathbf{z}}.
```

A particle launched on the axis at the centre of the bottle bounces harmonically
between the mirror points,

```math
z(t) = z_\mathrm{max}\sin(\omega t),\qquad
\omega = \sqrt{\frac{2\mu B_0}{mL^2}},\qquad
z_\mathrm{max} = L\sqrt{\frac{B_\mathrm{max}}{B_0}-1}.
```

This exercises the mirror force in both solvers. It also checks the two things
the guiding-centre approximation *assumes*: that the magnetic moment `μ` and the
kinetic energy are conserved along the full-orbit trajectory.

```@example gallery
m = 1.0; q = 1.0
B0 = 10.0; L = 0.4
function bottlefield(x, y, z, t)
    α = B0 * z / L^2
    return SVector(0.0, 0.0, 0.0), SVector(-x * α, -y * α, B0 + z * α)
end

vel0 = SVector(0.0, 0.1, 0.1)
rL = m * sqrt(vel0[1]^2 + vel0[2]^2) / (q * B0)
pos0 = SVector(-rL, 0.0, 0.0)
tspan = (0.0, 300.0)

E, B = bottlefield(pos0..., 0.0)
R0, vparal0, μ = get_guidingcentre(pos0, vel0, B, E, q, m)

sol_fo = solve(
    ODEProblem(lorentzforce!, [pos0..., vel0...], tspan,
        FullOrbitParams(charge=q, mass=m, electromagneticfield=bottlefield)),
    Tsit5(); reltol=1e-9, abstol=1e-12)
sol_gca = solve(
    ODEProblem(guidingcentreapproximation!, [R0..., vparal0], tspan,
        GCAParams(charge=q, mass=m, electromagneticfield=bottlefield, magneticmoment=μ)),
    Tsit5(); reltol=1e-9, abstol=1e-12)

Bmax = B0 * norm(vel0)^2 / (vel0[1]^2 + vel0[2]^2)
zmax = L * sqrt(Bmax / B0 - 1)
ω = sqrt(2μ * B0 / (m * L^2))
zanalytic(t) = zmax * sin(ω * t)

t = range(tspan..., length=3000)
μ_fo = [magneticmoment(sol_fo(tᵢ, idxs=1:3), sol_fo(tᵢ, idxs=4:6), tᵢ, m, bottlefield)
        for tᵢ in t]
Ek_fo = [kineticenergy(collect(sol_fo(tᵢ, idxs=4:6)), m) for tᵢ in t]

fig = Figure(size=(950, 620))
ax = Axis(fig[1, 1:2], xlabel="t [s]", ylabel="z [m]",
    title="Bounce motion in a magnetic bottle (first 60 s)")
lines!(ax, t, [sol_fo(tᵢ)[3] for tᵢ in t], label="full orbit")
lines!(ax, t, [sol_gca(tᵢ)[3] for tᵢ in t], linewidth=2.5, label="GCA")
lines!(ax, t, zanalytic.(t), linestyle=:dash, color=:black, label="analytic")
xlims!(ax, 0, 60)
axislegend(ax)

ax2 = Axis(fig[2, 1], xlabel="t [s]", ylabel="μ/μ₀ − 1",
    title="Magnetic moment (full orbit)")
lines!(ax2, t, μ_fo ./ μ .- 1)
ax3 = Axis(fig[2, 2], xlabel="t [s]", ylabel="E_k/E_k(0) − 1",
    title="Kinetic energy (full orbit)")
lines!(ax3, t, Ek_fo ./ first(Ek_fo) .- 1)
fig
```

```@example gallery
@printf("bounce amplitude, analytic : %.6f m\n", zmax)
@printf("bounce frequency, analytic : %.6f rad/s\n", ω)
@printf("max relative drift in μ    : %.3e\n", maximum(abs, μ_fo ./ μ .- 1))
@printf("max relative drift in E_k  : %.3e\n", maximum(abs, Ek_fo ./ first(Ek_fo) .- 1))
```

## Magnetic dipole: azimuthal drift period

A particle trapped in a dipole field drifts in azimuth. Walt (1994) gives the
drift period as

```math
T = \frac{2\pi q M}{m v^2 R_0}\left(1 - \tfrac{1}{3}\sin^{0.62}\alpha\right),
```

with `α` the equatorial pitch angle. Both solvers should complete the same
fraction of an orbit in a given time. This is the setup of Fig. 19 in Ripperda
et al. (2018).

```@example gallery
m = 1.0; q = 1.0
M = 60.0 * m / q                       # dipole moment, qM/m = 60
function dipolefield(x, y, z, t)
    α = M / (x^2 + y^2 + z^2)^(5 / 2)
    return SVector(0.0, 0.0, 0.0),
           SVector(3α * z * x, 3α * z * y, α * (2z^2 - x^2 - y^2))
end

R0 = SVector(1.0, 0.0, 0.0)
vparal = 0.5; vperp = 1.0
tspan = (0.0, 80.0)

E, B = dipolefield(R0..., 0.0)
μ = magneticmoment(vperp, m, norm(B))
# Start the full-orbit particle on a gyro-orbit about R0, at gyrophase π/2
u0 = collect(get_fullorbit(B, E, R0, vparal, μ, q, m, π / 2))

sol_fo = solve(
    ODEProblem(lorentzforce!, u0, tspan,
        FullOrbitParams(charge=q, mass=m, electromagneticfield=dipolefield)),
    Tsit5(); reltol=1e-8, abstol=1e-10)
sol_gca = solve(
    ODEProblem(guidingcentreapproximation!, [R0..., vparal], tspan,
        GCAParams(charge=q, mass=m, electromagneticfield=dipolefield, magneticmoment=μ)),
    Tsit5(); reltol=1e-6, abstol=1e-10)

v = sqrt(vparal^2 + vperp^2)
pitch = atan(vperp, vparal)
T_walt = 2π * q * M / (m * v^2 * R0[1]) * (1 - sin(pitch)^0.62 / 3)

function unwrapangle(φ)
    out = copy(φ); offset = 0.0
    for i in 2:length(φ)
        d = φ[i] - φ[i-1]
        d > π ? (offset -= 2π) : d < -π && (offset += 2π)
        out[i] = φ[i] + offset
    end
    return out
end

t = range(tspan..., length=4000)
fig = Figure(size=(950, 400))
ax = Axis(fig[1, 1], xlabel="x", ylabel="y", aspect=DataAspect(),
    title="Drift in the equatorial plane")
lines!(ax, [sol_fo(tᵢ)[1] for tᵢ in t], [sol_fo(tᵢ)[2] for tᵢ in t],
    linewidth=0.7, label="full orbit")
lines!(ax, [sol_gca(tᵢ)[1] for tᵢ in t], [sol_gca(tᵢ)[2] for tᵢ in t],
    linewidth=2, label="GCA")
lines!(ax, R0[1] .* cos.(2π .* t ./ T_walt), R0[1] .* sin.(2π .* t ./ T_walt),
    linestyle=:dash, color=:black, label="Walt (1994)")
axislegend(ax)

ax2 = Axis(fig[1, 2], xlabel="t [s]", ylabel="azimuth φ [rad]")
lines!(ax2, t, unwrapangle([atan(sol_fo(tᵢ)[2], sol_fo(tᵢ)[1]) for tᵢ in t]),
    label="full orbit")
lines!(ax2, t, unwrapangle([atan(sol_gca(tᵢ)[2], sol_gca(tᵢ)[1]) for tᵢ in t]),
    label="GCA")
lines!(ax2, t, 2π .* t ./ T_walt, linestyle=:dash, color=:black, label="Walt (1994)")
axislegend(ax2, position=:lt)
fig
```

```@example gallery
@printf("Walt drift period : %.4f s (%.3f of an orbit in %.0f s)\n",
        T_walt, last(tspan) / T_walt, last(tspan))
```

## When the guiding-centre approximation breaks down

The guiding-centre approximation is an expansion in the ratio of the Larmor
radius to the length scale over which the field varies. [`scalesratio`](@ref)
computes `r_L / L_B`; [`magneticcurvatureratio`](@ref) computes `r_L |κ|`. Both
are what the [`MagneticGradientCondition`](@ref),
[`MagneticCurvatureCondition`](@ref) and [`HybridSwitchCondition`](@ref)
callbacks watch.

Here the same magnetic bottle is run at four speeds. Faster particles have larger
Larmor radii and therefore worse scale separation, and the guiding centre
predicted by the approximation drifts away from the one reconstructed from the
exact orbit.

```@example gallery
tspan = (0.0, 25.0)
t = range(tspan..., length=1200)
ratios = Float64[]; errors = Float64[]; solutions = []

for vfac in (1.0, 4.0, 8.0, 16.0)
    vel0 = SVector(0.0, 0.1, 0.1) .* vfac
    rL = m * sqrt(vel0[1]^2 + vel0[2]^2) / (q * B0)
    pos0 = SVector(-rL, 0.0, 0.0)
    E, B = bottlefield(pos0..., 0.0)
    R0, vparal0, μ = get_guidingcentre(pos0, vel0, B, E, q, m)

    sol_fo = solve(
        ODEProblem(lorentzforce!, [pos0..., vel0...], tspan,
            FullOrbitParams(charge=q, mass=m, electromagneticfield=bottlefield)),
        Tsit5(); reltol=1e-10, abstol=1e-13)
    sol_gca = solve(
        ODEProblem(guidingcentreapproximation!, [R0..., vparal0], tspan,
            GCAParams(charge=q, mass=m,
                      electromagneticfield=bottlefield, magneticmoment=μ)),
        Tsit5(); reltol=1e-10, abstol=1e-13)

    gc = map(t) do tᵢ
        u = sol_fo(tᵢ)
        Eᵢ, Bᵢ = bottlefield(u[1], u[2], u[3], tᵢ)
        R, _, _ = get_guidingcentre(
            SVector(u[1], u[2], u[3]), SVector(u[4], u[5], u[6]), Bᵢ, Eᵢ, q, m)
        R
    end

    push!(ratios,
        maximum(abs, get_observable(sol_gca, :scalesratio; times=t, EoM="GCA")))
    push!(errors,
        maximum(norm(sol_gca(t[i])[1:3] .- gc[i]) for i in eachindex(t)) / L)
    push!(solutions, (vfac, sol_gca, gc))
end
nothing # hide
```

```@example gallery
fig = Figure(size=(950, 400))
for (col, k) in enumerate((1, 4))
    vfac, sol_gca, gc = solutions[k]
    ax = Axis(fig[1, col], xlabel="t [s]", ylabel="z [m]",
        title=@sprintf("v₀ × %d,  max r_L/L_B = %.2f", Int(vfac), ratios[k]))
    lines!(ax, t, [g[3] for g in gc], label="full orbit guiding centre")
    lines!(ax, t, [sol_gca(tᵢ)[3] for tᵢ in t], linestyle=:dash, linewidth=2,
        label="GCA solver")
    axislegend(ax, position=:lb)
end

ax = Axis(fig[1, 3], xscale=log10, yscale=log10,
    xlabel="max r_L / L_B", ylabel="max |R_GCA − R_FO| / L",
    title="Error vs scale separation")
scatterlines!(ax, ratios, errors, color=:black)
fig
```

The error grows steeply once `r_L/L_B` approaches unity — which is exactly what
the hybrid scheme below is for.

## Hybrid GCA ↔ full-orbit switching

[`hybridgcafo!`](@ref) integrates either equation and switches between them.
[`hybridswitchaffect!`](@ref) performs the switch: guiding centre → full orbit
requires inventing a gyrophase (drawn by `p.getphase`), and full orbit → guiding
centre discards it.

To verify the transformations in isolation, the switches are forced at fixed
times with a `PresetTimeCallback` rather than left to
[`HybridSwitchCondition`](@ref), and the result is compared against a pure
full-orbit reference in the same uniform E×B field.

```@example gallery
# Same uniform E x B field as the first example. Redefined here so the section
# does not depend on names bound further up the page.
q = -TraceParticles.e
m = TraceParticles.m_e
B0 = 3.37213e-4
vperp0 = 1.0e5
E0 = 0.25vperp0 * B0
hybridfield(x, y, z, t) = SVector(0.0, E0, 0.0), SVector(0.0, 0.0, B0)

gyroperiod = 2π * m / abs(q * B0)
tspan = (0.0, 12gyroperiod)
u0 = [0.0, 0.0, 0.0, vperp0, 0.0, 0.0]

switchtimes = [3gyroperiod, 6gyroperiod, 9gyroperiod]
hybridparams = HybridParams(
    charge=q, mass=m, electromagneticfield=hybridfield,
    initialeomid=EoMID.FullOrbit,
    getphase=integrator -> π / 2,     # deterministic, for reproducibility
)

sol_hyb = solve(
    ODEProblem(hybridgcafo!, copy(u0), tspan, hybridparams), Tsit5();
    reltol=1e-10, abstol=1e-14,
    callback=PresetTimeCallback(switchtimes, hybridswitchaffect!))

sol_ref = solve(
    ODEProblem(lorentzforce!, copy(u0), tspan,
        FullOrbitParams(charge=q, mass=m, electromagneticfield=hybridfield)),
    Tsit5(); reltol=1e-10, abstol=1e-14)
nothing # hide
```

Between `3T_c` and `6T_c`, and again after `9T_c`, the hybrid solution runs in
guiding-centre mode and the gyration disappears. The *guiding centre* tracks the
reference throughout — the full-orbit segments do not retrace the reference orbit
because the gyrophase is chosen arbitrarily on each switch back, which is
physically correct: the approximation discards that information.

The switch times are known, so which equation is active at any moment is known
too: full orbit on `[0, 3T_c)` and `[6T_c, 9T_c)`, guiding centre in between and
after. That lets the hybrid state be reduced to a guiding centre everywhere and
compared against the reference guiding centre directly.

```@example gallery
t = range(tspan..., length=3000)
mm = 1e3   # metres -> millimetres

guidingcentre(u, tᵢ) = begin
    Eᵢ, Bᵢ = hybridfield(u[1], u[2], u[3], tᵢ)
    R, _, _ = get_guidingcentre(
        SVector(u[1], u[2], u[3]), SVector(u[4], u[5], u[6]), Bᵢ, Eᵢ, q, m)
    R
end

infullorbit(tᵢ) = tᵢ < switchtimes[1] ||
                  (switchtimes[2] <= tᵢ < switchtimes[3])

gcref = [guidingcentre(sol_ref(tᵢ), tᵢ) for tᵢ in t]
# In guiding-centre mode the state already *is* the guiding centre.
gchyb = [infullorbit(tᵢ) ? guidingcentre(sol_hyb(tᵢ), tᵢ) :
         SVector(sol_hyb(tᵢ)[1], sol_hyb(tᵢ)[2], sol_hyb(tᵢ)[3]) for tᵢ in t]

fig = Figure(size=(950, 380))
ax = Axis(fig[1, 1], xlabel="x [mm]", ylabel="y [mm]", title="Hybrid vs full orbit")
lines!(ax, mm .* [sol_ref(tᵢ)[1] for tᵢ in t], mm .* [sol_ref(tᵢ)[2] for tᵢ in t],
    label="full orbit")
lines!(ax, mm .* [sol_hyb(tᵢ)[1] for tᵢ in t], mm .* [sol_hyb(tᵢ)[2] for tᵢ in t],
    linewidth=2, label="hybrid")
axislegend(ax, position=:lt)

ax2 = Axis(fig[1, 2], xlabel="t / T_c", ylabel="guiding centre x [mm]",
    title="Guiding centre through the switches")
lines!(ax2, t ./ gyroperiod, mm .* [g[1] for g in gcref], label="full orbit")
lines!(ax2, t ./ gyroperiod, mm .* [g[1] for g in gchyb],
    linestyle=:dash, linewidth=2, label="hybrid")
vlines!(ax2, switchtimes ./ gyroperiod, linestyle=:dash, color=:black)
axislegend(ax2, position=:lt)
fig
```

The gyration in the left panel disappears during the guiding-centre segments and
reappears at an arbitrary phase — the straight radial jumps — while the guiding
centre in the right panel is continuous across every switch.

```@example gallery
gcerr = maximum(norm(gchyb[i] .- gcref[i]) for i in eachindex(t))
@printf("switches performed               : %d\n", hybridparams.nswitches)
@printf("Larmor radius                    : %.3e m\n", larmorradius(m, vperp0, q, B0))
@printf("max |hybrid − FO guiding centre| : %.3e m\n", gcerr)
```

## Speiser orbits in a current sheet

The classic non-adiabatic problem: a Harris-type current sheet
`B = b(η, −x/d, 0)` with a reconnection electric field `E = E_z ẑ`. Particles
entering the sheet are trapped, oscillate across it with an amplitude damped as
`t^{-1/4}`, are accelerated along `E`, and are finally ejected by the guide field
`ηb` at

```math
t_\mathrm{eject} = \frac{\pi m}{q\eta b}.
```

The guiding-centre approximation has no chance here — the Larmor radius is
comparable to the sheet width — which is precisely why the full-orbit solver
exists.

```@example gallery
m = TraceParticles.m_p
q = TraceParticles.e
η = 0.025      # guide-field fraction
d = 1e-4       # current sheet half-width [m]
b = 1e-2       # characteristic field strength [T]
Ez = -1e-2 * 1e7 * b

speiserfield(x, y, z, t) =
    SVector(0.0, 0.0, Ez), SVector(η * b, -x * b / d, 0.0)

t_eject = π * m / (q * η * b)
tspan = (0.0, 1.15t_eject)

sol = solve(
    ODEProblem(lorentzforce!, [1e-2, 0.0, 0.0, 0.0, 0.0, 0.0], tspan,
        FullOrbitParams(charge=q, mass=m, electromagneticfield=speiserfield)),
    Tsit5(); reltol=1e-9, abstol=1e-14, maxiters=Int(1e8))
nothing # hide
```

The damping law is checked against the peaks of the oscillation rather than the
raw signal, since `|x|` touches zero at every sheet crossing. It describes the
*trapped* phase only: once the guide field has turned the particle far enough
that it leaves the sheet, the amplitude grows again. The comparison below is
therefore restricted to the first half of the trapping time, and the point where
it stops applying is marked.

```@example gallery
t = range(1e-8, last(tspan), length=8000)
x = [sol(tᵢ)[1] for tᵢ in t]
absx = abs.(x)

# local maxima of |x| = the oscillation envelope
allpeaks = [i for i in 2:length(t)-1 if absx[i] > absx[i-1] && absx[i] >= absx[i+1]]
allpeaks = filter(i -> absx[i] > 1e-6, allpeaks)
trapped = filter(i -> t[i] < 0.5t_eject, allpeaks)
tref = t[trapped[length(trapped)÷4]]; xref = absx[trapped[length(trapped)÷4]]

fig = Figure(size=(950, 620))
ax = Axis(fig[1, 1], xlabel="t [s]", ylabel="|x| [m]",
    xscale=log10, yscale=log10,
    title="Oscillation damped as t^(−1/4) while trapped")
scatter!(ax, t[trapped], absx[trapped], markersize=5, label="envelope of |x|")
lines!(ax, t[trapped], xref .* (t[trapped] ./ tref) .^ (-0.25),
    color=:black, linestyle=:dash, linewidth=2, label="∝ t^(−1/4)")
axislegend(ax, position=:lb)

ax2 = Axis(fig[1, 2], xlabel="t [s]", ylabel="v_y [m/s]", title="Ejection")
lines!(ax2, t, [sol(tᵢ)[5] for tᵢ in t])
vlines!(ax2, [t_eject], linestyle=:dash, color=:black, label="πm/(qηb)")
axislegend(ax2, position=:lb)

ax3 = Axis(fig[2, 1], xlabel="t [s]", ylabel="E_k [eV]",
    title="Energy gained from E∥")
lines!(ax3, t, [kineticenergy(collect(sol(tᵢ, idxs=4:6)), m) * J2eV for tᵢ in t])
vlines!(ax3, [t_eject], linestyle=:dash, color=:black)

# The particle is accelerated hundreds of metres along z while oscillating over
# a fraction of a millimetre in x, so the x-z trajectory is only legible once
# the particle has entered the sheet.
tenter = t[findfirst(<(2e-4), absx)]
tz = range(tenter, tenter + 0.15t_eject, length=6000)
ax4 = Axis(fig[2, 2], xlabel="x [mm]", ylabel="z [m]",
    title="Trajectory once inside the sheet")
lines!(ax4, [1e3sol(tᵢ)[1] for tᵢ in tz], [sol(tᵢ)[3] for tᵢ in tz])
fig
```

```@example gallery
@printf("predicted ejection time : %.4e s\n", t_eject)
@printf("final kinetic energy    : %.4e eV\n",
        kineticenergy(collect(sol(last(tspan), idxs=4:6)), m) * J2eV)
```

## Sampling initial conditions

The sampler behind [`SampleICsFromMHD`](@ref) relies on two primitives, both of
which can be checked against their target distributions directly.

[`maxwellianvelocitysample`](@ref) draws one velocity component from a
Maxwell–Boltzmann distribution at temperature `T`, i.e. a zero-mean Gaussian of
width `σ = √(k_B T/m)`. [`rejectionsample`](@ref) draws a position and time from
an arbitrary distribution over `(x, y, z, t)`, given an upper bound on its
value. Both take the floating-point type to sample in as their second argument,
which is how a single-precision run stays single precision.

```@example gallery
rng = Xoshiro(1234)

T = 1e6                       # K
m = TraceParticles.m_p
σ = sqrt(TraceParticles.k_B * T / m)
v = [maxwellianvelocitysample(rng, Float64, T, m) for _ in 1:200_000]

# The sampler is four-dimensional; holding z and t fixed leaves a 2D target.
target(x, y, z, t) = exp(-((x - 0.3)^2 + (y + 0.2)^2) / (2 * 0.15^2)) +
                     0.6exp(-((x + 0.5)^2 + (y - 0.4)^2) / (2 * 0.25^2))
maxvalue = 1.05               # must be ≥ max(target) over the domain

npoints = 60_000
xs = Vector{Float64}(undef, npoints); ys = similar(xs); nrejections = 0
for i in 1:npoints
    xs[i], ys[i], _, _, r = rejectionsample(
        rng, Float64, target, maxvalue,
        -1.0, 1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    )
    global nrejections += r
end

fig = Figure(size=(950, 380))
ax = Axis(fig[1, 1], xlabel="v [m/s]", ylabel="probability density",
    title="maxwellianvelocitysample")
hist!(ax, v, bins=120, normalization=:pdf, label="samples")
vv = range(-4σ, 4σ, length=400)
lines!(ax, vv, @.(exp(-0.5(vv / σ)^2) / (σ * sqrt(2π))),
    color=:black, linestyle=:dash, linewidth=2, label="Maxwell–Boltzmann")
axislegend(ax)

ax2 = Axis(fig[1, 2], xlabel="x", ylabel="y", aspect=DataAspect(),
    title="rejectionsample (contours = target)")
edges = range(-1, 1, length=81)
counts = zeros(80, 80)
for i in 1:npoints
    ix = clamp(searchsortedlast(edges, xs[i]), 1, 80)
    iy = clamp(searchsortedlast(edges, ys[i]), 1, 80)
    counts[ix, iy] += 1
end
heatmap!(ax2, edges, edges, counts)
grid = range(-1, 1, length=200)
contour!(ax2, grid, grid, [target(x, y, 0.0, 0.0) for x in grid, y in grid],
    color=:white, levels=6)
fig
```

```@example gallery
@printf("σ analytic       : %.5e m/s\n", σ)
@printf("σ sampled        : %.5e m/s\n", std(v))
@printf("sample mean      : %.5e m/s\n", mean(v))
@printf("acceptance rate  : %.3f\n", npoints / (npoints + nrejections))
```

## References

- Ripperda, B. et al. (2018), *A comprehensive comparison of relativistic
  particle integrators*, ApJS 235, 21 — magnetic bottle and dipole setups.
- Walt, M. (1994), *Introduction to Geomagnetically Trapped Radiation* —
  dipole drift period.
- Speiser, T. W. (1965), *Particle trajectories in model current sheets*,
  JGR 70, 4219.
- Chen, F. F. (2016), *Introduction to Plasma Physics and Controlled Fusion* —
  drift expressions.

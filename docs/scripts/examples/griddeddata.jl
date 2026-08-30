using TraceParticles: lorentzforce!, electron_mass, electron_charge
using OrdinaryDiffEq: ODEProblem, solve
using CairoMakie
CairoMakie.activate!(type="png") # hide

# The electromagnetic field should to be a callable with signature
# (x,y,z,t), returning two 3D vectors: (E, B).
B0 = 5e-4  # Tesla
E0 = 1e2  # V/m
emfield(x, y, z, t) = [0.0, E0, 0.0], [0.0, 0.0, B0]

charge = electron_charge # Coloumb
mass = electron_mass     # kg
gyroperiod = 2π * mass / abs(charge * B0) # Seconds

# The initial conditions is a 6D state vector: [x, y, z, vx, vy, vz]
u0 = [0.0, 0.0, 0.0, 1.0e6, 0.0, 0.0] # [m, m, m, m/s, m/s, m/s]

prob = ODEProblem(
    lorentzforce!,
    u0,
    (0.0, 8gyroperiod), # Time span of simulation [seconds]
    (; charge=charge, mass=mass, electromagneticfield=emfield),
)
sol = solve(prob)

# We can interpolate the dense solution to any tᵢ with sol(tᵢ)
t = range(0.0, 8gyroperiod, length=2000)
fig, ax, _ = lines(
    1e3 .* [sol(tᵢ)[1] for tᵢ in t], # x-position in mm
    1e3 .* [sol(tᵢ)[2] for tᵢ in t]; # y-position in mm
    axis=(
        xlabel="x [mm]",
        ylabel="y [mm]",
        title="An electron gyrating and E×B-drifting"
    ),
)
fig

using TraceParticles: guidingcentreapproximation!, get_guidingcentre

E, B = emfield(u0[1], u0[2], u0[3], 0.0)
R0, vparal0, μ = get_guidingcentre(
    [u0[1], u0[2], u0[3]],   # position
    [u0[4], u0[5], u0[6]],   # velocity
    B, E, charge, mass,
)

prob_gca = ODEProblem(
    guidingcentreapproximation!,
    [R0..., vparal0],
    (0.0, 8gyroperiod),
    (; 
        charge=charge,
        mass=mass,
        electromagneticfield=emfield,
        magneticmoment=μ
    ),
)
sol_gca = solve(prob_gca)

lines!(ax, 1e3 .* [sol_gca(tᵢ)[1] for tᵢ in t],
            1e3 .* [sol_gca(tᵢ)[2] for tᵢ in t],
       linewidth=3, label="GCA")
axislegend(ax, position=:lt)
fig


using Interpolations

xaxis = 0:0.01:1
yaxis = 0:0.01:1

# Gridded components
bz = [B0*(1 - 4.2*x) for x in xaxis, y in yaxis]
ey = [E0 for x in xaxis, y in yaxis]
bz_itp = linear_interpolation((xaxis, yaxis), bz, extrapolation_bc=Flat())
ey_itp = linear_interpolation((xaxis, yaxis), ey, extrapolation_bc=Flat())

emfield(x, y, _, _) = [0, ey_itp(x, y), 0], [0, 0, bz_itp(x, y)]

sol = solve(prob; reltol=1e-8)
sol_gca = solve(prob_gca)

fig, ax, _ = lines(
    1e3 .* [sol(tᵢ)[1] for tᵢ in t],
    1e3 .* [sol(tᵢ)[2] for tᵢ in t];
    axis=( #
        xlabel="x [mm]",
        ylabel="y [mm]",
        title="An electron gyrating and E×B-drifting"
    ), label="Full orbit"
) # hide
lines!(ax, 1e3 .* [sol_gca(tᵢ)[1] for tᵢ in t],
            1e3 .* [sol_gca(tᵢ)[2] for tᵢ in t],
       linewidth=3, label="GCA") # hide
axislegend(ax, position=:lt) # hide
# ...
fig

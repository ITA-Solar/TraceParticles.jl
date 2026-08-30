```@meta
CurrentModule = TraceParticles
```

# Examples

Complete, runnable recipes. Every code block on this page is executed when the
documentation is built, so the numbers and figures you see are the ones the
current version of the package produces.

```@setup examples
using TraceParticles
# The conditions, affects and physics helpers below are part of the package but
# deliberately not exported: the high-level interface builds callbacks from
# `callback_specs` instead. Import them explicitly to use them directly.
using TraceParticles: get_guidingcentre, get_fullorbit,
    OutOfBoundsCondition, outofboundsaffect!,
    MagneticGradientCondition, magneticgradientaffect!
using OrdinaryDiffEq
using SciMLBase: EnsembleSerial, EnsembleThreads
using DiffEqCallbacks
using Interpolations
using StaticArrays
using LinearAlgebra
using Statistics
using Random
using Printf
using DataFrames
using CairoMakie
CairoMakie.activate!(type="png")
```

## An ensemble of particles

[`PredefinedICs`](@ref) is the `prob_func` to use when you already know where
every particle starts. You give it three vectors — initial states, time spans and
parameter containers — and it hands out the `i`-th element of each. Everything
else is the standard `EnsembleProblem` interface, so `EnsembleThreads()` and
`EnsembleDistributed()` work as usual.

```@example examples
const B₀ = 3.37213e-4                     # T
emfield(x, y, z, t) = SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, B₀)

charge = -TraceParticles.e
mass = TraceParticles.m_e
gyroperiod = 2π * mass / abs(charge * B₀)

npart = 5
u0 = [[0.0, 0.0, 0.0, 1e5 * (1 + 0.2i), 0.0, 1e4 * i] for i in 1:npart]
tspans = [(0.0, 3gyroperiod) for _ in 1:npart]
params = [FullOrbitParams(charge=charge, mass=mass, electromagneticfield=emfield)
          for _ in 1:npart]

prob = ODEProblem(lorentzforce!, u0[1], tspans[1], params[1])
eprob = EnsembleProblem(prob; prob_func=PredefinedICs(u0, tspans, params))
sim = solve(eprob, Tsit5(), EnsembleThreads(); trajectories=npart, reltol=1e-9)

fig = Figure(size=(700, 500))
ax = Axis3(fig[1, 1], xlabel="x [mm]", ylabel="y [mm]", zlabel="z [mm]",
           title="Five electrons in a uniform field")
for (i, s) in enumerate(sim.u)
    t = range(s.t[1], s.t[end], length=600)
    lines!(ax, 1e3 .* [s(tᵢ)[1] for tᵢ in t],
               1e3 .* [s(tᵢ)[2] for tᵢ in t],
               1e3 .* [s(tᵢ)[3] for tᵢ in t], label="particle $i")
end
axislegend(ax)
fig
```

Different perpendicular speeds give different Larmor radii, different parallel
speeds give different pitch of the helix. The individual `ODESolution`s live in
`sim.u`:

```@example examples
[s.retcode for s in sim.u]
```

## A field on a grid

Fields do not have to be analytical. Anything callable as `(x, y, z, t)` that
returns `(E, B)` will do, and [`ElectromagneticFieldInterpolator`](@ref) turns a
set of interpolators into exactly that. It expects the six components in the
order `[bx, by, bz, ex, ey, ez]`; the wrappers
[`StaticInterpolation`](@ref), [`XZInterpolation`](@ref) and
[`StaticXZInterpolation`](@ref) adapt interpolators built over fewer dimensions
to the four-argument interface.

Here we tabulate `B = (B₀ + a y) ẑ` on a grid and let a cubic spline reconstruct
it, so the resulting ∇B drift can be compared with the analytical answer
`v_∇B = -μ a / (q B)`.

```@example examples
xx = range(0.0, 1.0, length=64)
yy = range(0.0, 1.0, length=64)
zz = range(0.0, 1.0, length=2)

Bg, a = 6.0, 10.0
components = [(x, y, z) -> 0.0,          # bx
              (x, y, z) -> 0.0,          # by
              (x, y, z) -> Bg + a * y,   # bz
              (x, y, z) -> 0.0,          # ex
              (x, y, z) -> 0.0,          # ey
              (x, y, z) -> 0.0]          # ez

itps = [cubic_spline_interpolation(
            (xx, yy, zz),
            [f(x, y, z) for x in xx, y in yy, z in zz],
            extrapolation_bc=Flat(),
        ) for f in components]

gridfield = TraceParticles.ElectromagneticFieldInterpolator([TraceParticles.StaticInterpolation(i) for i in itps])
gridfield(0.5, 0.5, 0.5, 0.0)   # (E, B) at the box centre
```

```@example examples
q, m = 1.0, 1.0
pos0 = SVector(0.5, 0.5, 0.5)
vel0 = SVector(0.0, 0.1, 0.0)
E, B = gridfield(pos0..., 0.0)
R₀, vparal₀, μ = get_guidingcentre(pos0, vel0, B, E, q, m)

sol = solve(
    ODEProblem(guidingcentreapproximation!, [R₀..., vparal₀], (0.0, 5.0),
               GCAParams(charge=q, mass=m,
                         electromagneticfield=gridfield, magneticmoment=μ)),
    Tsit5(); reltol=1e-8,
)

analytical = R₀[1] - μ * a / (q * norm(B)) * 5.0
@printf("numerical x(5) = %.8f\nanalytical x(5) = %.8f\n", sol.u[end][1], analytical)
```

The spline reproduces the drift to the accuracy of the tabulated field. Note
that `extrapolation_bc=Flat()` matters: the guiding-centre equations differentiate
the field with `ForwardDiff`, and a particle that steps a hair outside the grid
must still get a finite answer.

## Stopping particles with callbacks

Callbacks are how particles leave a simulation. Each condition is a callable
struct holding its tolerance; each affect writes a
[`TraceParticles.TerminationCode`](@ref) into
the parameters and calls `terminate!`. Because the parameter containers are
mutable, the reason a particle stopped is available afterwards on `sol.prob.p`.

```@example examples
gradfield(x, y, z, t) = SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, Bg + a * y)

E, B = gradfield(pos0..., 0.0)
R₀, vparal₀, μ = get_guidingcentre(pos0, vel0, B, E, q, m)
p = GCAParams(charge=q, mass=m, electromagneticfield=gradfield, magneticmoment=μ)

callbacks = CallbackSet(
    # stop when state variable 1 (x) leaves [0.4, 0.6]
    DiscreteCallback(OutOfBoundsCondition(((0.4, 0.6),), [1]), outofboundsaffect!),
    # stop when the Larmor radius exceeds 10 % of the field's scale length
    DiscreteCallback(MagneticGradientCondition(1e-1), magneticgradientaffect!),
)

sol = solve(
    ODEProblem(guidingcentreapproximation!, [R₀..., vparal₀], (0.0, 1e4), p),
    Tsit5(); reltol=1e-8, callback=callbacks,
)

@printf("retcode          = %s\nterminationcode  = %s\nstopped at t     = %.2f\nx at stop        = %.4f\n",
        sol.retcode, p.terminationcode, last(sol.t), sol.u[end][1])
```

The particle drifted out of the box long before its Larmor radius grew large, so
the termination code is `OutOfBounds`. The available conditions are
[`OutOfBoundsCondition`](@ref), [`MagneticGradientCondition`](@ref),
[`MagneticCurvatureCondition`](@ref), [`ParallelElectricFieldCondition`](@ref),
[`RelativisticCondition`](@ref) (and its GCA and hybrid variants) and
[`HybridSwitchCondition`](@ref).

## Post-processing a solution

[`get_observable`](@ref) evaluates derived quantities along a solution, using its
dense output, so you are not restricted to the time steps the solver happened to
take. Quantities that depend on the equation of motion need `EoM="GCA"` or
`EoM="FO"` unless the solution carries hybrid switch information, in which case
they are inferred.

```@example examples
B_bottle, L = 10.0, 0.4
function bottle(x, y, z, t)
    α = B_bottle * z / L^2
    SVector(0.0, 0.0, 0.0), SVector(-x * α, -y * α, B_bottle + z * α)
end

vel0 = SVector(0.0, 0.1, 0.1)
r_L = m * 0.1 / (q * B_bottle)
pos0 = SVector(-r_L, 0.0, 0.0)
E, B = bottle(pos0..., 0.0)
R₀, vparal₀, μ = get_guidingcentre(pos0, vel0, B, E, q, m)

solb = solve(
    ODEProblem(guidingcentreapproximation!, [R₀..., vparal₀], (0.0, 30.0),
               GCAParams(charge=q, mass=m,
                         electromagneticfield=bottle, magneticmoment=μ)),
    Tsit5(); reltol=1e-9,
)

times = range(0.0, 30.0, length=600)
fig = Figure(size=(760, 520))
for (k, (sym, label)) in enumerate([(:z, "z"), (:pitchangle, "pitch angle [rad]"),
                                    (:vperp, "v⊥"), (:larmorradius, "r_L")])
    ax = Axis(fig[fld(k - 1, 2)+1, mod(k - 1, 2)+1], xlabel="t", ylabel=label)
    lines!(ax, times, get_observable(solb, sym; times=times, EoM="GCA"))
end
fig
```

For a single time, [`GCAState`](@ref) computes everything derived from a
guiding-centre state at once — all four drifts, the accelerations, the energy,
the Larmor radius, the gyrofrequency, the pitch angle and the characteristic
field length:

```@example examples
gs = GCAState(collect(solb(10.0, idxs=1:4)), 10.0, q, m, μ, bottle)
(energy=gs.energy, r_L=gs.r_L, cospitchangle=gs.β, ω_c=gs.ω_c,
 ∇Bdrift=gs.∇Bdrift, mirroracc=gs.mirroracc, L_B=gs.L_B)
```

A `Vector{GCAState}` converts directly to a `DataFrame`, which is convenient for
turning a whole ensemble into a table.

## Sampling an ensemble and streaming it to disk

This is the production workflow: draw initial conditions from a distribution,
reduce each solution to a handful of numbers, and write the results to HDF5 batch
by batch so that memory usage stays flat no matter how many particles you run.

The sampler draws positions and times from a *proposal* distribution by
rejection sampling, draws a velocity from a Maxwell–Boltzmann distribution at the
local temperature, and assigns each particle an importance weight
`target(x,y,z,t) / proposal(x,y,z,t)`. Sampling uniformly and weighting by the
target is often the simplest choice; sampling from the target itself makes all
weights equal to one.

[`SampleGCA`](@ref) wraps it into a `prob_func`. It is the guiding-centre one of
three — [`SampleFullOrbit`](@ref) and [`SampleHybrid`](@ref) build the
6-component state and their own parameter containers instead. The high-level
interface picks between them for you from
[`SampleICsFromMHD`](@ref); here the distributions are closures rather than
files on disk, so the sampler is built directly.

```@example examples
outdir = mktempdir()
expname = "demo"
mkpath(joinpath(outdir, expname))

B_s = 1e-3
uniformfield(x, y, z, t) = SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, B_s)

proposal(x, y, z, t) = 1.0                                        # sample uniformly
density(x, y, z, t) = exp(-((x - 0.5)^2 + (y - 0.5)^2) / (2 * 0.2^2))
temperature(x, y, z, t) = 1e6                                     # K

npart, batchsize = 5_000, 1_000
T_c = 2π * mass / abs(charge * B_s)

sampler = TraceParticles.MHDSampler(
    density, proposal, temperature,       # target, proposal, temperature
    0.0, 1.0,
    0.0, 1.0,
    0.0, 1.0,
    0.0, 0.0,                             # t₀ bounds
    1.0,                                  # upper bound on the proposal
    Float64,                              # precision to sample in
)
probfunc = SampleGCA(sampler, 20T_c)      # 20T_c is the final time

reduction = SaveBatchAsHDF5(
    expdir=outdir, expname=expname, batchsize=batchsize,
    nbatches=cld(npart, batchsize), verbose=false,
)

baseprob = ODEProblem(
    guidingcentreapproximation!, zeros(4), (0.0, 20T_c),
    GCAParams(charge=charge, mass=mass,
              electromagneticfield=uniformfield, magneticmoment=0.0),
)
eprob = EnsembleProblem(baseprob; prob_func=probfunc,
                        output_func=output_func_lightweight,
                        reduction=reduction, safetycopy=false)

sim = solve(eprob, Tsit5(), EnsembleSerial(); trajectories=npart,
            batch_size=batchsize, reltol=1e-8, abstol=1e-10)

filename = get_filename(reduction)
h5_getbatchnames(filename)
```

The base problem's `u0` and `tspan` are placeholders — `prob_func` overwrites
both for every trajectory. What matters in it is the equation of motion and the
type of the parameter container.

Reading the results back gives a `DataFrame`:

```@example examples
df = h5_getall(filename)
first(df[:, [:x0, :y0, :vparal0, :magneticmoment, :weight, :nt, :terminationcode]], 5)
```

Because the sampling was uniform, the physical distribution lives entirely in the
weights. Histogramming `x0` with and without them shows the importance weighting
doing its job:

```@example examples
edges = range(0.0, 1.0, length=26)
centres = 0.5 .* (edges[1:end-1] .+ edges[2:end])
idx = [clamp(searchsortedlast(edges, x), 1, length(centres)) for x in df.x0]
raw = [count(==(k), idx) for k in eachindex(centres)]
weighted = [sum(df.weight[idx.==k]; init=0.0) for k in eachindex(centres)]

fig = Figure(size=(700, 380))
ax = Axis(fig[1, 1], xlabel="x₀", ylabel="normalised counts",
          title="Uniform sampling, Gaussian target")
barplot!(ax, centres, raw ./ sum(raw), color=(:grey, 0.6), label="unweighted")
barplot!(ax, centres, weighted ./ sum(weighted), color=(:dodgerblue, 0.6),
         label="importance weighted")
target = [exp(-(x - 0.5)^2 / (2 * 0.2^2)) for x in centres]
lines!(ax, centres, target ./ sum(target), color=:black, linewidth=2,
       label="target ∝ exp(-(x-0.5)²/2σ²)")
axislegend(ax, position=:lt)
fig
```

[`batch_statistics`](@ref) summarises a batch — this is what `SaveBatchAsHDF5`
logs as it writes:

```@example examples
print(TraceParticles.batch_statistics(df))
```

Other readers are available for parts of the file:
[`h5_getbatch`](@ref), [`h5_getdataset`](@ref), and [`h5_getenergies`](@ref).

## The high-level interface

For production runs you rarely assemble the `EnsembleProblem` by hand.
[`TraceParticlesParameters`](@ref) collects every knob into one
keyword-constructed struct — equation of motion, particle count, tolerances,
which field file to read, where initial conditions come from, what terminates a
particle, where to write output — and validates the combination on
construction. [`TraceParticlesProblem`](@ref) turns that into a fully assembled
problem, and `solve` runs it.

The field is read from a JLD2 file holding six interpolation objects in the order
`[bx, by, bz, ex, ey, ez]`; [`create_bifrost_itps`](@ref) produces one from a
Bifrost snapshot.

```julia
using TraceParticles
using OrdinaryDiffEq

params = TraceParticlesParameters(
    eom=guidingcentreapproximation!,
    npart=1_000_000,
    tf=1.0,
    mass=TraceParticles.m_e,
    charge=-TraceParticles.e,
    alg=Tsit5(),
    maxiters=Int(1e7),
    reltol=1e-8, abstol=1e-10,
    expname="myexperiment",
    datadir="/path/to/output",
    batchsize=10_000,
    emfield_file="/path/to/emfield.jld2",

    # Sample initial conditions from the snapshot during the run
    prob_func=SampleICsFromMHD(
        tg_file="/path/to/temperature.jld2",
        target_distr_file="/path/to/density.jld2",
        xbounds=(0.0, 1.0),
        ybounds=(0.0, 1.0),
        zbounds=(0.0, 1.0),
        t0bounds=(0.0, 0.0),
    ),

    # Terminate particles that leave the box in z, or become relativistic
    callback_specs=(
        OutOfBounds(z=(0.0, 1.0)),
        Relativistic(fraction=0.12),
    ),
)

prob = TraceParticlesProblem(params)
sim = solve(prob, params)
```

`TraceParticlesProblem` sets up logging to both `stderr` and
`datadir/expname/expname.log.txt` — set `log2file=false` in the parameters to
keep the logger you pass as the `logger` keyword instead. `solve` chooses
distributed, threaded or serial execution from the number of available processes
and threads, streams the batches to HDF5 in `datadir/expname`, and finally
computes the initial and final energies of every particle and appends them to the
file.

Once a run has finished, particles re-integrated with different settings usin
[`rerun`](@ref).

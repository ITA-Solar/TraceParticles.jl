```@meta
CurrentModule = TraceParticles
```

# Concepts

This page describes the pieces `TraceParticles` provides and how they plug into
a `DifferentialEquations.jl` problem. Everything here is ordinary SciML
machinery: if you know how to build an `ODEProblem` and an `EnsembleProblem`,
you already know how to use this package.

## The electromagnetic field

Every equation of motion, callback and physics function in `TraceParticles`
takes the external field as a **callable**

```julia
E, B = electromagneticfield(x, y, z, t)
```

that returns the **electric field first** and the magnetic field second, each as
a 3-component vector. That is the entire interface. It can be an anonymous
function over an analytical expression, a closure over interpolation objects, or
one of the wrappers below.

Returning `StaticArrays.SVector`s is strongly recommended: the guiding-centre
equations differentiate the field with `ForwardDiff`, and static vectors keep
that allocation-free.

```julia
# Analytical field
emfield(x, y, z, t) = SVector(0.0, E0, 0.0), SVector(0.0, 0.0, B0)
```

For gridded fields the package ships [`ElectromagneticFieldInterpolator`](@ref),
which wraps either a single interpolator returning all six components, or a
6-element vector of per-component interpolators ordered
`[bx, by, bz, ex, ey, ez]`, and presents them through the `(x, y, z, t) -> (E, B)`
interface:

```julia
emfield = ElectromagneticFieldInterpolator(itps)   # itps[1] -> bx, ..., itps[6] -> ez
```

Interpolators built over fewer dimensions are adapted with the small wrappers
[`StaticInterpolation`](@ref) (a time-independent 3D field, `(x,y,z,t) -> itp(x,y,z)`),
[`XZInterpolation`](@ref) (a time-dependent 2.5D field, `(x,y,z,t) -> itp(x,z,t)`)
and [`StaticXZInterpolation`](@ref) (`(x,y,z,t) -> itp(x,z)`).

## State vectors and equations of motion

Which equation of motion you pick determines the length and meaning of the state
vector `u`.

| Equation of motion | `u` | Length |
| :--- | :--- | :--- |
| [`lorentzforce!`](@ref) | `[x, y, z, vx, vy, vz]` | 6 |
| [`guidingcentreapproximation!`](@ref) | `[Rx, Ry, Rz, v∥]` | 4 |
| [`hybridgcafo!`](@ref) | `[·, ·, ·, ·, ·, ·]` | 6 (interpretation depends on `p.eomid`) |

[`lorentzforce!`](@ref) integrates the full Lorentz force, resolving the
gyration. It is exact but has to resolve the gyroperiod, which is expensive when
the region of interest is many gyroradii across.

[`guidingcentreapproximation!`](@ref) integrates the drift motion of the
*guiding centre* instead. The gyration is averaged out and replaced by the
adiabatic invariant `μ = m v⊥² / 2B`, the magnetic moment, which is carried in
the parameters rather than in the state vector. Included are the E×B, ∇B,
curvature and polarisation drifts, plus the mirror force, parallel electric
acceleration and the first-order Fermi term. The approximation is valid when the
Larmor radius is small compared with the scales over which the field varies —
see [`scalesratio`](@ref) and [`magneticcurvatureratio`](@ref), and the
[Verification gallery](@ref) for how the error actually behaves.

[`hybridgcafo!`](@ref) runs *either* of the two, switching at run time based on
`p.eomid`. Because the state vector must accommodate both, it is always 6
components long; in guiding-centre mode only the first four carry meaning.

Converting between the two representations is done by [`get_guidingcentre`](@ref)
and [`get_fullorbit`](@ref) (and their in-place `!` variants). Going from full
orbit to guiding centre discards the gyrophase; going back requires you to
supply one, since 5 quantities have to become 6.

## Parameter containers

The `ODEProblem` parameters `p` carry the particle's charge, mass, the field, and
any per-particle bookkeeping. `TraceParticles` provides **mutable** structs for
this, because callbacks need to write into them mid-integration — a `NamedTuple`
would not do.

- [`FullOrbitParams`](@ref) — charge, mass, field, statistical weight, number of
  rejections, termination code.
- [`GCAParams`](@ref) — the same plus the magnetic moment `magneticmoment`.
- [`HybridParams`](@ref) — the same plus `eomid` (which equation of motion is
  currently active), an RNG and a `getphase` function used to pick a gyrophase
  when switching back to full orbit, and a switch counter.
- [`HybridParamsWithDetection`](@ref) — as `HybridParams`, but additionally
  records the time, the new equation of motion and the new magnetic moment at
  every switch, so the trajectory can be reconstructed in post-processing.

All of them expose the fields the equations of motion read: `charge`, `mass`,
`electromagneticfield` and — where relevant — `magneticmoment`. Any object with
those fields works, including a plain `NamedTuple`, as long as nothing needs to
mutate it.

## Enumerations

Three [`EnumX`](https://github.com/fredrikekre/EnumX.jl) enums appear in the
parameters and in saved output:

- [`EoMID`](@ref) — `FullOrbit` or `GuidingCentreApproximation`. Selects the
  active equation of motion in the hybrid scheme.
- [`TraceParticles.TerminationCode`](@ref) — why a particle stopped:
  `NotTerminated`, `OutOfBounds`, `MagneticGradient`, `MagneticCurvature`,
  `ParallelElectricField`, `Relativistic` or `MaxSwitches`. Terminating callback
  affects set it, and the output functions save it as an `Int32`.
- [`HybridScheme`](@ref) — whether a sampling `prob_func` should build hybrid
  parameters (`No`, `Yes`, `YesSaveSwitches`).

## Callbacks

Callbacks are how a particle is stopped, or how the hybrid scheme decides to
change equations. They follow the SciML split into a *condition* (a callable
returning `true` when something should happen) and an *affect* (what to do).

The conditions are callable structs holding their tolerances:

- [`OutOfBoundsCondition`](@ref) — a state variable left a box.
- [`MagneticGradientCondition`](@ref) — `r_L / L_B` exceeded a tolerance.
- [`MagneticCurvatureCondition`](@ref) — `r_L |κ|` exceeded a tolerance.
- [`ParallelElectricFieldCondition`](@ref) — `E∥ / E⊥` exceeded a tolerance.
- [`RelativisticCondition`](@ref), [`RelativisticConditionGCA`](@ref),
  [`RelativisticConditionHybrid`](@ref) — kinetic energy exceeded a given
  fraction of the rest energy.
- [`HybridSwitchCondition`](@ref) — the guiding-centre approximation just became
  invalid (or valid again).

The terminating affects — [`outofboundsaffect!`](@ref),
[`magneticgradientaffect!`](@ref), [`magneticcurvatureaffect!`](@ref),
[`parallelelectricfieldaffect!`](@ref), [`relativisticaffect!`](@ref) — each set
`p.terminationcode` and call `terminate!`. The switching affects
[`hybridswitchaffect!`](@ref) and [`hybridswitchaffect_withdetection!`](@ref)
transform the state vector between representations instead of stopping.

```julia
using DiffEqCallbacks

cb = DiscreteCallback(
    OutOfBoundsCondition(((0.0, 1.0), (-1.0, 1.0)), [1, 3]),  # bound x and z
    outofboundsaffect!,
)
```

### Callback specifiers

Pairing a condition with the right affect by hand is error-prone, and which
relativistic condition is correct depends on the equations of motion. The
high-level interface therefore takes *specifiers* instead — one per criterion,
each naming only what it needs:

- [`OutOfBounds`](@ref) — bounds on any of `x`, `y`, `z`.
- [`Relativistic`](@ref) — the fraction of the rest energy at which to stop.
- [`HighMagneticGradient`](@ref) — the tolerance on `r_L / L_B`.

```julia
callback_specs = (
    OutOfBounds(x = (0.0, 1.0), z = (-1.0, 1.0)),
    Relativistic(fraction = 0.02),
)
```

A specifier validates its own parameters on construction — an `OutOfBounds`
cannot exist without bounds — and [`create_callback`](@ref) builds the actual
callback from it, choosing e.g. between [`RelativisticCondition`](@ref),
[`RelativisticConditionGCA`](@ref) and [`RelativisticConditionHybrid`](@ref)
according to the equations of motion. Adding a criterion means adding a
[`TraceParticles.CallbackSpec`](@ref) subtype and one `create_callback` method,
and touching nothing that already exists.

The hybrid GCA/full-orbit switch is not a specifier: it is part of the
equations of motion rather than an optional add-on, and is configured by the
`gradient_tolerance`, `curvature_tolerance`, `eparal_ratio` and
`switchback_tolerance` fields of [`TraceParticlesParameters`](@ref) whenever
`eom = hybridgcafo!`.

## The ensemble

A test-particle run is an
[`EnsembleProblem`](https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/)
built from three user-supplied functions. `TraceParticles` provides ready-made
implementations of all three.

### `prob_func` — where each particle starts

Called once per trajectory to `remake` the base problem with that particle's
initial condition, time span and parameters.

- [`PredefinedICs`](@ref) takes vectors of `u0`, `tspan` and `params` and hands
  out the `i`-th of each.
- [`SampleFullOrbit`](@ref), [`SampleGCA`](@ref) and [`SampleHybrid`](@ref) draw
  each particle from a sampler — one problem function per equation of motion,
  since they differ only in the state vector and parameter container they build.

The sampler behind the last three draws a position and time from a proposal
distribution by rejection sampling, draws a velocity from a Maxwell–Boltzmann
distribution at the local temperature, and assigns an importance weight from the
ratio of target to proposal distribution. This is the workhorse for sampling an
MHD snapshot.

At the high level you do not pick between them. `prob_func` in
[`TraceParticlesParameters`](@ref) takes a *specification* —
[`SampleICsFromMHD`](@ref) for sampling during the run, [`ICsFromFile`](@ref)
for reading initial conditions from disk — and
[`TraceParticlesProblem`](@ref) resolves it into the matching problem function
once it knows the equations of motion. Passing your own `(prob, ctx)` function
instead bypasses this entirely; it is used unchanged.

### `output_func` — what is kept from each solution

An `ODESolution` carries a lot of metadata. For a ten-million-particle ensemble
that overhead dominates memory, so the output function reduces each solution to a
`NamedTuple` before it is stored.

- [`output_func_lightweight`](@ref) — initial and final guiding-centre state,
  charge, mass, magnetic moment, weight, times, and return/termination codes.
- [`output_func_lightweight_hybrid`](@ref) — the same for a 6-component hybrid
  state, plus switch counts and equation-of-motion IDs.
- [`output_func_max_lightweight`](@ref) — as `output_func_lightweight`, plus the
  largest Larmor radius, smallest characteristic field length and largest scales
  ratio encountered along the trajectory.

### `reduction` — what happens to a batch

Called once per batch of trajectories.

- [`SaveBatchAsHDF5`](@ref) writes each batch into its own group of an HDF5 file
  and logs [`batch_statistics`](@ref) for it, so a long run streams to disk
  instead of accumulating in memory. [`get_filename`](@ref) returns the path it
  writes to.

## Post-processing

[`GCAState`](@ref) takes a guiding-centre state and the field and computes
everything derived from it: all four drifts, the accelerations, the kinetic
energy, the Larmor radius, gyrofrequency, pitch angle and the characteristic
field length. A `Vector{GCAState}` converts straight to a `DataFrame`.

[`get_observable`](@ref) is the general accessor for an `ODESolution`. It knows
about primary variables (`:x`, `:vparal`, …), derived quantities (`:energy`,
`:pitchangle`, `:larmorradius`, `:scalesratio`, `:magneticcurvatureratio`,
`:magneticmoment`, `:gyrofrequency`, …), guiding-centre-only quantities
(`:fermi`, `:betatron`, `:parallelenergy`, …) and field quantities (`:bfield`,
`:eparal`, `:exbdrift`, `:characteristicfieldlength`, …). It evaluates them at
whatever times you ask for, using the solution's dense interpolation.

Results written to HDF5 are read back with [`h5_getall`](@ref),
[`h5_getbatch`](@ref), [`h5_getdataset`](@ref), and [`h5_getenergies`](@ref).

## The high-level workflow

For production runs, [`TraceParticlesParameters`](@ref) collects every knob —
equation of motion, particle count, tolerances, field file, where initial
conditions come from, which criteria terminate a particle — into one
keyword-constructed struct. Most requirements are encoded in the types of the
specifications it holds; what is left spans several fields and is validated on
construction, so a misconfigured run fails there rather
than after an hour of integration. [`TraceParticlesProblem`](@ref) turns it into
a fully assembled `EnsembleProblem` (loading the field, resolving the
`prob_func`, building the callbacks, output function and reduction), and
`solve(prob, params)` runs it with logging and automatic selection of serial,
threaded or distributed execution.

See [Examples](@ref) for how this fits together in practice.

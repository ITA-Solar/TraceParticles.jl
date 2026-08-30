```@meta
CurrentModule = TraceParticles
```

# Concepts

## The electromagnetic field
Every equation of motion (and some callback and physics functions) in TraceParticles.jl take the external electromagnetic field as a **callable**

```julia
E, B = electromagneticfield(x, y, z, t)
```

that returns the **electric field first** and the magnetic field second, each as a 3-component vector.
It can be an anonymous function over an analytical expression, or a closure over interpolation objects.

Returning `StaticArrays.SVector`s is recommended for performance; the guiding-centre equations differentiate the field with [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl), and static vectors keep that allocation-free.

```julia
# Analytical field
emfield(x, y, z, t) = SVector(0.0, E0, 0.0), SVector(0.0, 0.0, B0)
```

## State vectors and equations of motion

Which equation of motion you pick determines the length and meaning of the state
vector `u`.

| Equation of motion | `u` | Expected fields of `p` |
| :--- | :--- | :--- |
| [`lorentzforce!`](@ref) | `[x, y, z, vx, vy, vz]` | `mass`, `charge`, `electromagneticfield` |
| [`guidingcentreapproximation!`](@ref) | `[Rx, Ry, Rz, v∥]` | `mass`, `charge`, `electromagneticfield`, `magneticmoment` |
| [`hybridgcafo!`](@ref) | `[x, y, z, vx, vy, vz]` or `[Rx, Ry, Rz, v∥, (unused), (unused)]` | `mass`, `charge`, `electromagneticfield`, `magneticmoment`, `eomid` |

[`lorentzforce!`](@ref) integrates the Lorentz force directly, resolving the
gyration.

[`guidingcentreapproximation!`](@ref) integrates the drift motion of the *guiding centre* instead.
The gyration is averaged out and replaced by the adiabatic invariant $\mu$, the magnetic moment, which is carried in the parameters rather than in the state vector.
The state variables are the guiding centre position and velocity parallel to the magnetic field.
The equations model the parallel acceleration of the guiding centre and its perpendicular drifts.
Included are the $\mathbf{E}\times\mathbf{B}$, $\nabla B$, curvature and polarisation drifts.
The approximation is valid when the Larmor radius is small compared with the scales over which the field varies —
see [`scalesratio`](@ref) and [`magneticcurvatureratio`](@ref).

[`hybridgcafo!`](@ref) runs *either* of the two, switching at run time based on `p.eomid`.
Because the state vector must accommodate both, it is always 6 components long; in guiding centre mode only the first four carry meaning.
Converting between the two representations is done by [`get_guidingcentre`](@ref) and [`get_fullorbit`](@ref) (and their in-place `!` variants).
Going from full orbit to guiding centre discards the gyrophase; going back requires you to supply one, since 5 variables have to become 6.
[`HybridSwitchCondition`](@ref) and [`hybridswitchaffect!`](@ref) can be used together as a [`callback`](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/) to detect the need to switch, to transform particle state, and to mutate `p.eomid`.
`p.eomid` must be of the [`EnumX`](https://github.com/fredrikekre/EnumX.jl)-type `EoMID.T`.
```@example eomid
using TraceParticles: EoMID
EoMID.T
```

## Parameter containers
The `ODEProblem`'s parameters `p` is expected to carry the particle's charge, mass, the electromagnetic field, and any per-particle bookkeeping.
TraceParticles.jl provides **mutable** structs for this.
- [`FullOrbitParams`](@ref) — charge, mass, electromagnetic field, statistical weight, number of rejections, termination code.
- [`GCAParams`](@ref) — the same plus the magnetic moment `magneticmoment`.
- [`HybridParams`](@ref) — the same plus `eomid` (which equation of motion is
  currently active), an `rng`, a `getphase` function used to pick a gyrophase
  when switching back to full orbit, and a switch counter `nswitches`.
- [`HybridParamsWithDetection`](@ref) — as `HybridParams`, but additionally
  holds the arrays `timeatswitch`, `eomidafterswitch`, and `magneticmomentafterswitch`, so the trajectory can be reconstructed in post-processing.
  It also holds `maxswitches`, so the particle can be terminated if it switches to often.

## [Callbacks](@id cancepts-callbacks)
TraceParticles.jl provides a set of conditions and affects that can be used in [callbacks](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/).
Callbacks are passed to `DifferentialEquations.solve` to handles event during the integration.

### Conditions
- [`OutOfBoundsCondition`](@ref) — a state variable left the domain
- [`RelativisticCondition`](@ref), [`RelativisticConditionGCA`](@ref),
  [`RelativisticConditionHybrid`](@ref) — kinetic energy exceeded a given
  fraction of the rest energy.
- [`MagneticGradientCondition`](@ref) — $r_L/L_B$ exceeded a tolerance.
- [`MagneticCurvatureCondition`](@ref) — $r_L/\kappa$ exceeded a tolerance.
- [`ParallelElectricFieldCondition`](@ref) — $E_\parallel / E_\perp$ exceeded a tolerance.
- [`HybridSwitchCondition`](@ref) — the guiding-centre approximation just became
  invalid (or valid again).

### Affects
#### Affects that terminate the particle
The following affects terminate the integration and sets its parameters
`p.terminationcode` to a [`TraceParticles.TerminationCode`](@ref) type:
- [`outofboundsaffect!`](@ref) -- sets to `TerminationCode.OutOfBounds`,
- [`magneticgradientaffect!`](@ref) -- sets to `TerminationCode.MagneticGradient`,
- [`magneticcurvatureaffect!`](@ref) -- sets to `TerminationCode.MagneticCurvature`,
- [`parallelelectricfieldaffect!`](@ref) -- sets to `TerminationCode.ParallelElectricField`,
- [`relativisticaffect!`](@ref) -- sets to `TerminationCode.Relativistic,

#### Hybrid-scheme affects
[`hybridswitchaffect!`](@ref) and [`hybridswitchaffect_withdetection!`](@ref) transform the state vector between full orbit and guiding centre representations and mutates `p.eomid`.
The latter also mutates the arrays `p.timeatswitch`, `p.eomidafterswitch`, and `p.magneticmomentafterswitch`, and terminates if `p.nswitches` exceeds `p.maxswitches`.

## Tools for running particle ensembles
[`EnsembleProblem`](https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/)s builds from three user-supplied functions.
`TraceParticles` provides implementations of all three.

### `prob_func`
Called once per trajectory to `remake` the base problem with that particle's initial condition, time span and parameters.
- [`PredefinedICs`](@ref) takes vectors of `u0`, `tspan` and `params` and hands
  out the `i`-th of each.
- [`SampleFullOrbit`](@ref), [`SampleGCA`](@ref) and [`SampleHybrid`](@ref) draw
  each particle from a sampler — one problem function per equation of motion,
  since they differ in the state vector and parameter container they build.

### `output_func`
Decides what to keep of each particle solution.
The following output function reduces each particle trajectory solution to a `NamedTuple`:
- [`output_func_lightweight_gca`](@ref) — initial and final guiding centre state, charge, mass, magnetic moment, weight, times, and return/termination codes.
- [`output_func_lightweight_hybrid`](@ref) — the same for a 6-component hybrid state, plus switch counts and equation-of-motion IDs.
- [`output_func_max_lightweight`](@ref) — as `output_func_lightweight`, plus the largest Larmor radius, smallest characteristic field length and largest scales ratio encountered along the trajectory.

### `reduction`
Called once per batch of trajectories.
- [`SaveBatchAsHDF5`](@ref) writes each batch into its own group of an HDF5 file and logs [`batch_statistics`](@ref) for it.
  [`get_filename`](@ref) returns the path it writes to.

## The high-level workflow
For convenience, [`TraceParticlesParameters`](@ref) collects every knob necessary to run an ensemble simulation with the TraceParticles.jl-tools.
A [`TraceParticlesProblem`](@ref) turns the `TraceParticlesParameters` into a fully assembled `EnsembleProblem` (loading the electromagnetic field, resolving the `prob_func`, building the callbacks, output function and reduction), and [`TraceParticles.solve`](@ref) runs it with logging and automatic selection of serial, threaded or distributed execution.
```julia
params = TraceParticlesParameters(;
    # ...keyword-defined experiment parameters...
)
tpprob = TraceParticlesProblem(params)
sol = solve(tpprob, params)
```
Alternative solver options can be used with
```julia
sol = solve(
    tpprob;
    alg=Tsit5(),
    trajectories=20_000_000,
    batch_size=1_000_000,
    solver_kwargs...
)
```
This high-level workflow is useful when producing multiple large ensemble simulations in a structured manner.
The `TraceParticlesParameters`-struct defines the experiment, is easily adjustable and promotes reproducability.
The possibly time consuming problem-building is outsourced to the `TraceParticlesProblem` constructor.
See below or the [high-level interface example](@ref interface-example).

### The electromagnetic-field in a file
`TraceParticlesProblem` expects the electromagnetic field to be stored in a [JLD2.jl](https://github.com/juliaio/jld2.jl)-file, defined by the parameter `emfield_file::String`.

### The `prob_func`-specifiers
[`TraceParticlesParameters`](@ref) accepts any type for the `prob_func`-field.
However, TraceParticles.jl provides `prob_func`-*specifiers* that `TraceParticlesProblem` automatically resolves.
- [`SampleICsFromMHD`](@ref) for sampling initial conditions during the run, using the [`MHDSampler`](@ref).
- [`ICsFromFile`](@ref) for reading initial conditions from disk.

### The `output_func` and `reduction` 
[`TraceParticlesParameters`](@ref) accepts any type for the `output_func`- and `reduction`-fields.
If nothing is given, `TraceParticlesProblem` uses [`output_func_lightweight`](@ref outputfunc-api) and [`SaveBatchAsHDF5`](@ref), respectively.

### Callback-specifiers
The `TraceParticlesParameters`-field `callbacks` accepts a tuple of any type, but only `SciML.DECallback` and [`TraceParticles.CallbackSpec`](@ref callbackspec-api) are accepted by `TraceParticlesProblem`.
`CallbackSpec`-types were made to standardise the pairing of conditions and affects, and to ease the definition of callbacks that depended on other experiment-parameters.
A `CallbackSpec` is defined as an element in `TraceParticlesProblem.callbacks`, but is resolved into a `SciML.DECallback` by `TraceParticlesProblem` via
```julia
create_callback(cs::CallbackSpec, params)
```

## Post-processing
[`get_observable`](@ref) is the general accessor for an `ODESolution`.
It knows about primary variables (`:x`, `:vparal`, …), derived quantities (`:energy`,`:pitchangle`, `:larmorradius`, `:scalesratio`, `:magneticcurvatureratio`, `:magneticmoment`, `:gyrofrequency`, …), guiding-centre-only quantities (`:fermi`, `:betatron`, `:parallelenergy`, …) and field quantities (`:bfield`, `:eparal`, `:exbdrift`, `:characteristicfieldlength`, …).
It evaluates them at whatever times you ask for, using the solution's dense interpolation.

Results written to HDF5 are read back with [`h5_getall`](@ref), [`h5_getbatch`](@ref), [`h5_getdataset`](@ref), and [`h5_getenergies`](@ref).

See [Examples](@ref) for how everything fits together in practice.

```@meta
CurrentModule = TraceParticles
```

# TraceParticles.jl

With `TraceParticles` you can *trace* charged *particles*  in electromagnetic fields (often called test particles).
The particles move independently of each other by solving the Lorentz force and do not feed back to the electromagnetic field.
To solve the equations the package builds on [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/), which
provides adaptive time-stepping, dense output, callback and parallelism of ensembles of particles.

## What the package provides

| Tool | Description |
| :--- | :--- |
| [Equations of motion](@ref eom-api) | `lorentzforce!`, `guidingcentreapproximation!` and the hybrid `hybridgcafo!` |
| [Parameter containers](@ref params-api) | `FullOrbitParams`, `GCAParams`, `HybridParams` — mutable `ODEProblem` parameters |
| [Problem functions](@ref probfunc-api) | `prob_func`s that draw initial conditions, including importance-weighted sampling of particles from a fluid |
| [Output functions](@ref outputfunc-api) | `output_func`s that keep the ensemble's memory footprint small |
| [Reduction functions](@ref reduction-api) | `reduction`s that store particle batches to HDF5 |
| [Callbacks](@ref callbacks-api) | Conditions and affects that terminate particles, or switch between equations of motion |
| [Physics](@ref physics-api) | Drifts, Larmor radius, magnetic moment, pitch angle, guiding-centre ↔ full-orbit transforms |
| [Post-processing](@ref io-api) | `get_observable`, and HDF5 readers |
| [High-level interface](@ref interface-api) | Container for paremeters and functions for easy assembly and running of a trace particle experiment |

Many of the tools are general and can be used to trace particles with custom
equations of motion.

```@contents
Pages = ["examples.md", "concepts.md"]
Depth = 5
```

### Usage examples
The [Examples](@ref)-page demonstrate the use of `TraceParticles.jl`-tools,
covering simple use of the provided equations of motion, to using the high-level interface for running an ensemble of particles with importance sampled initial condition, hybrid equations, and efficient data-storage.
* [Using the equations of motion](@ref equations-example)
* [Using callbacks](@ref callbacks-example)
* [Using statistical sampling of initial conditions](@ref sampling-example)
* [Using the hybrid equations scheme](@ref hybrid-example)
* [Efficient storage of particle ensembles](@ref ensembles-example)
* [Using the high-level interface](@ref interface-example))
* [Using post-processing tools](@ref post-processing-example))

### Design concepts
[Concepts](@ref) explains design concepts of `TraceParticles.jl`, like the electromagnetic field interface, the parameter containers, the state vectors of the equations, and how the ensemble pieces fit together.

### Package API
The [API](@ref) lists every exported function and type.

### Numerical verification gallery
The [Verification gallery](@ref) reproduces analytically solvable problems.

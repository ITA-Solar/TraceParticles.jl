# TestParticles.jl documentation

`TestParticles` was created to embed test particles in magnetohydrodynamic
environments, or any electromagnetic field. Test particles move independently
and their motion is governed by the Lorentz force. This code uses [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) to solve the particles' equations of motion.

## Basic usage
`TestParticles` provides tools to simulate test particle ensembles using the
[`EnsembleProblem`](https://docs.sciml.ai/DiffEqDocs/dev/features/ensemble/) from `DifferentialEquations`. The tools provided by `TestParticles.jl` are:
- [Equations of motion (EoMs)](http://localhost:8000/api/#Equations-of-motion)
- [`prob_func`s for draw initial conditions and importance weighing.](http://localhost:8000/api/#Problem-functions)
- [`output_func`s for managing the solutions.](http://localhost:8000/api/#Output-functions)
- [`reduction` functions for storing data.](http://localhost:8000/api/#Reduction-functions)
- [`callback`s for e.g. determine the validity of EoMs.](http://localhost:8000/api/#Callbacks)
- [Post-processing tools for fetching and analysing the data.](http://localhost:8000/api/#I/O)

See [examples](http://localhost:8000/examples/) for simulations in both
analytical and gridded electromagnetic fields, with and wihtout statistical
sampling of the particle initial conditions.


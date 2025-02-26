# TestParticles.jl documentation

`TestParticles` was created to embed test particles in magnetohydrodynamic
environments, or any electromagnetic field. Test particles move independently
and their motion is governed by the Lorentz force. This code uses [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) to solve the particles' equations of motion.

## Basic usage
`TestParticles` provides tools to simulate test particle ensembles using the
[`EnsembleProblem`](https://docs.sciml.ai/DiffEqDocs/dev/features/ensemble/) from DifferentialEquations.jl. The tools provided by `TestParticles.jl` are:
- [Equations of motion (EoMs)](https://ita-solar.github.io/TestParticles.jl/dev/api/#Equations-of-motion)
- [`prob_func`s for drawing initial conditions and importance weighing.](https://ita-solar.github.io/TestParticles.jl/dev/api/#Problem-functions)
- [`output_func`s for managing the solutions.](https://ita-solar.github.io/TestParticles.jl/dev/api/#Output-functions)
- [`reduction` functions for storing data.](https://ita-solar.github.io/TestParticles.jl/dev/api/#Reduction-functions)
- [`callback`s for e.g. determine the validity of EoMs.](https://ita-solar.github.io/TestParticles.jl/dev/api/#Callbacks)
- [Post-processing tools for fetching and analysing the data.](https://ita-solar.github.io/TestParticles.jl/dev/api/#I/O)

See [examples](https://ita-solar.github.io/TestParticles.jl/dev/examples/)
for simulations in both analytical and gridded electromagnetic fields, with and
wihtout statistical sampling of the particle initial conditions.
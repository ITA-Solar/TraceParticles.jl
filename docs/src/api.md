```@meta
CurrentModule = TraceParticles
CollapsedDocStrings = true
```

# API

## [Equations of motion](@id eom-api)
```@docs
lorentzforce!
guidingcentreapproximation!
hybridgcafo!
```

## [Parameter containers](@id params-api)
```@docs
FullOrbitParams
GCAParams
HybridParams
HybridParamsWithDetection
```

## Enumerations
```@docs
EoMID
TerminationCode
HybridScheme
```

## [Problem functions](@id probfunc-api)
The `prob_func` of an
[`EnsembleProblem`](https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/)
-- how each individual particle problem is created.
```@docs
PredefinedICs
SampleFullOrbit
SampleGCA
SampleHybrid
```

## [Output functions](@id outputfunc-api)
The `output_func` of an `EnsembleProblem` -- what is kept from each problem solution.
```@docs
output_func_lightweight_gca
output_func_lightweight_fo
output_func_lightweight_hybrid
output_func_max_lightweight
```

## [Reduction functions](@id reduction-api)
The `reduction` of an `EnsembleProblem` -- what happens to a batch of solved
problems.
```@docs
SaveBatchAsHDF5
get_filename
batch_statistics
```

## [Callbacks](@id callbacks-api)
Conditions and affects to use in [callbacks](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/), for handling events during the particle trajectories.
Conditions decide *when* something happens; affects decide *what*.

### Conditions
```@docs
OutOfBoundsCondition
MagneticGradientCondition
MagneticCurvatureCondition
ParallelElectricFieldCondition
RelativisticCondition
RelativisticConditionGCA
RelativisticConditionHybrid
HybridSwitchCondition
```

### Affects
```@docs
outofboundsaffect!
magneticgradientaffect!
magneticcurvatureaffect!
parallelelectricfieldaffect!
relativisticaffect!
hybridswitchaffect!
hybridswitchaffect_withdetection!
```

## [Physics](@id physics-api)
```@docs
get_guidingcentre
get_fullorbit
scalesratio
magneticcurvatureratio
```

## [Statistics](@id statistics-api)
```@docs
mhdsample
rejectionsample
maxwellianvelocitysample
MHDSampler
```

## [Post-processing and I/O](@id io-api)
```@docs
get_observable
h5_getall
h5_getbatch
h5_getdataset
h5_getenergies
h5_getbatchnames
save_energy
```

## Particle reruns
```@docs
rerun
```

## [The high-level interface](@id interface-api)
```@docs
TraceParticlesParameters
TraceParticlesProblem
solve
```

### `prob_func`-specifications
`prob_func`-specifications resolves into one of the [problem functions above](@ref probfunc-api) via `init_probfunc`.
```@docs
SampleICsFromMHD
ICsFromFile
init_probfunc
```

### [Callback specifiers](@id callbackspec-api)
```@docs
CallbackSpec
create_callback
OutOfBounds
Relativistic
HighMagneticGradient
```

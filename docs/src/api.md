```@meta
CurrentModule = TraceParticles
CollapsedDocStrings = true
```

# API

Everything the package exports, grouped by what it is for. See
[Concepts](@ref) for how the pieces fit together.

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

The `prob_func` of an `EnsembleProblem`: what each trajectory starts from.

[`TraceParticlesParameters`](@ref) takes a *specification* — `SampleICsFromMHD`
or `ICsFromFile` — which [`TraceParticlesProblem`](@ref) resolves into one of
the problem functions below. Any function callable as `(prob, ctx)` may be
passed instead, in which case it is used unchanged.

### Specifications

```@docs
SampleICsFromMHD
ICsFromFile
```

### Problem functions

```@docs
PredefinedICs
SampleFullOrbit
SampleGCA
SampleHybrid
```

## [Output functions](@id outputfunc-api)

The `output_func` of an `EnsembleProblem`: what is kept from each solution.

```@docs
output_func_lightweight_gca
output_func_lightweight_fo
output_func_lightweight_hybrid
output_func_max_lightweight
```

## [Reduction functions](@id reduction-api)

The `reduction` of an `EnsembleProblem`: what happens to a finished batch.

```@docs
SaveBatchAsHDF5
get_filename
batch_statistics
```

## [Callback specifiers](@id callbackspec-api)

What `callback_specs` in [`TraceParticlesParameters`](@ref) holds. Each
specifier validates its own parameters on construction and knows how to build
its callback, so an inconsistent configuration fails there rather than during
the run.

```@docs
CallbackSpec
create_callback
OutOfBounds
Relativistic
HighMagneticGradient
```

## [Callbacks](@id callbacks-api)

The conditions and affects the specifiers above are built from. Conditions
decide *when* something happens; affects decide *what*.

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

## Statistics

```@docs
mhdsample
rejectionsample
maxwellianvelocitysample
```

## [Post-processing and I/O](@id io-api)
```@docs
get_observable
h5_getall
h5_getbatch
h5_getdataset
h5_getenergies
h5_getbatchnames
create_diffeq_ic
ElectromagneticFieldInterpolator
StaticInterpolation
XZInterpolation
StaticXZInterpolation
GCAState
create_bifrost_itps
save_energy
```

### Particle reruns

```@docs
rerun
```

## [The high-level interface](@id interface-api)

```@docs
TraceParticlesParameters
TraceParticlesProblem
solve
```

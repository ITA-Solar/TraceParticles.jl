```@meta
CollapsedDocStrings = true
```

## Equations of motion
```@docs
guidingcentreapproximation!
lorentzforce!
hybridgcafo!
```

## Callbacks
```@docs
OutOfBoundsCondition
MagneticGradientCondition
MagneticCurvatureCondition
ParallelElectricFieldCondition
RelativisticConditionGCA
outofboundsaffect!
magneticgradientaffect!
magneticcurvatureaffect!
parallelelectricfieldaffect!
relativisticaffect!
HybridSwitchCondition
hybridswitchaffect!
TerminationCode
```

## Problem functions
```@docs
DposMBvel
PposPvel
```

## Output functions
```@docs
output_func_max_lightweight
output_func_lightweight
```

## Reduction functions
```@docs
SaveBatchAsHDF5
get_filename
```

## Parameter structs
```@docs
GCAParams
HybridParams
```

## I/O
```@docs
h5_getall
h5_getbatch
h5_getdataset
h5_getenergies
h5_getinitialstate
h5_getfinalstate
save_gcastates
create_bifrost_itps
```

## Dataprocessing
```@docs
GCAState
get_observable
initialstate_idxs
initialstate_terminationcode
finalstate_terminationcode
initialstate_maxiters
particlemaskfunction
particlemask
generate_probdistr
replaceparticles
rerun
```

## Physics
```@docs
get_fullorbit
get_fullorbit!
get_guidingcentre
get_guidingcentre!
kineticenergy
perpendicular_velocity
larmorradius
characteristicfieldlength
scalesratio
magneticmoment
```

## Statistics
```@docs
maxwellianvelocitysample
rejectionsample
binmap
```

## Utilities
```@docs
ElectromagneticFieldInterpolator
StaticInterpolation
StaticXZInterpolation
XZInterpolation
```

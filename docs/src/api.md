```@meta
CollapsedDocStrings = true
```

# Functions
```@docs
rerun_nonthermals
rerun_maxiters
rerun_idxs
```

## Equations of motion
```@docs
guidingcentreapproximation!
lorentzforce!
```

## Statistics
```@docs
maxwellianvelocitysample
rejectionsample
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

# Reduction functions
```@docs
SaveBatchAsHDF5
get_filename
```

## I/O
```@docs
h5_getall
h5_getbatch
h5_getdataset
h5_getenergies
h5_getinitialstate
save_gcastates
```

## Dataprocessing
```@docs
GCAState
get_observable
```

# Physics
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

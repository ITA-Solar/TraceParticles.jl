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
gca_2Dxz!
lorentzforce!
```

## Statistics
```@docs
maxwellianvelocitysample
```

## Callbacks
```@docs
OutOfDomainCondition_2Dxz
RelativisticConditionGCA_2Dxz
GCABreakDownCondition_2Dxz
outofdomainaffect!
relativisticaffect!
gcabreakdownaffect!
```

## Problem functions
```@docs
UposMBvel
DposMBvel
PposPvel
```

## Output functions
```@docs
output_func_max_lightweight
```

# Reduction functions
```@docs
SaveBatchAsCSV
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
get_stateidx
get_exbdrift
get_fermi
get_betatron
get_polarisationacc
get_driftenergy
get_perpenergy
get_parallelenergy
binmap
```
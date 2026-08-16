"""
An EnumX module containing termination codes. The codes are set by callback
affects that terminate the solver. E.g. the callback affect `outofboundsaffect!`
will terminate the integration and set the parameter `terminationcode` to
`TerminationCode.OutOfBounds`.
"""

"""
    TerminationCode
NotTerminated          = 0
OutOfBounds            = 1
MagneticGradient       = 2
MagneticCurvature      = 3
ParallelElectricField  = 4
Relativistic           = 5
MaxSwitches            = 6
"""
@enumx TerminationCode begin
    NotTerminated         # = 0
    OutOfBounds           # = 1
    MagneticGradient      # = 2
    MagneticCurvature     # = 3
    ParallelElectricField # = 4
    Relativistic          # = 5
    MaxSwitches           # = 6
end

"""
    EoMID
FullOrbit                   = 0
GuidingCentreApproximation  = 1
"""
@enumx EoMID begin
    FullOrbit                  # = 0
    GuidingCentreApproximation # = 1
end

"""
    HybridScheme
Default      = 0
SaveSwitches = 1

Whether a hybrid GCA/full-orbit run records each switch between the two sets of
equations. `SaveSwitches` selects [`HybridParamsWithDetection`](@ref) over
[`HybridParams`](@ref); it is chosen by `save_switchinfo` in the parameters.
"""
@enumx HybridScheme begin
    Default      # = 0
    SaveSwitches # = 1
end

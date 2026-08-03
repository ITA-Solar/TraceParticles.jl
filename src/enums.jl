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
No              = 0
Yes             = 1
YesSaveSwitches = 2
"""
@enumx HybridScheme begin
    No              # = 0
    Yes             # = 1
    YesSaveSwitches # = 2
end

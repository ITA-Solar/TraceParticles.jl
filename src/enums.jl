"""
An EnumX module containing termination codes. The codes are set by callback
affects that terminate the solver. E.g. the callback affect `outofboundsaffect!`
will terminate the integration and set the parameter `terminationcode` to
`TerminationCode.OutOfBounds`.
"""
@enumx TerminationCode begin
    NotTerminated
    OutOfBounds
    MagneticGradient
    MagneticCurvature
    ParallelElectricField
    Relativistic
end



#-------------------------------------------------------------------------------
#           callback_specifiers.jl
#-------------------------------------------------------------------------------
# This file defines structs that contain and validate the parameters used for
# constructing specific callbacks. The struct can be passed to
# `TraceParticlesParameters` and `TraceParticlesProblem` will construct the
# actual callback with its condition and affect.
#
# A callback specifier should be a subtype of the abstract type `CallbackSpec`
# and have a corresponding method
#
#   create_callback(::<:CallbackSpec, ::TraceParticlesParameters)
#
# Adding a callback means adding a struct and that one method *here*, and
# touching nothing that already exists.
#
# Encoding the requirements in the types means an inconsistent configuration is
# largely unrepresentable rather than something `validate` has to catch: an
# `OutOfBounds` cannot exist without bounds.
#-------------------------------------------------------------------------------

"""
    CallbackSpec
Supertype of the structs specifying callbacks. See e.g. [`OutOfBounds`](@ref),
[`Relativistic`](@ref) and [`HighMagneticGradient`](@ref).
"""
abstract type CallbackSpec end

"""
    create_callback(spec::CallbackSpec, p::TraceParticlesParameters)
Build the callback that enforces `spec`, matching callback conditions and
affects. `p` supplies the run-wide quantities a criterion may need, such as the
particle mass and the equations of motion.
"""
function create_callback end


#-------------------------------------------------------------------------------
# Termination callbacks
#-------------------------------------------------------------------------------

"""
    OutOfBounds(; x=nothing, y=nothing, z=nothing, save_positions=(false, true))

Terminate a particle that leaves the given spatial bounds. Each of `x`, `y` and
`z` is either `nothing` (that axis is not checked) or a `(lower, upper)` tuple.
At least one axis must be bounded.

The bounds are only ever compared against the state vector, never multiplied
with it, so they are promoted to a common float type among themselves and do
not affect the precision the particles are integrated in.

```julia
OutOfBounds(x = (1.58e7, 1.63e7), z = (-6.6e6, -6.1e6))
```
"""
struct OutOfBounds{N,T<:AbstractFloat} <: CallbackSpec
    u_idxs::NTuple{N,Int}
    bounds::NTuple{N,Tuple{T,T}}
    save_positions::Tuple{Bool,Bool}
end
function OutOfBounds(
    ;
    x=nothing,
    y=nothing,
    z=nothing,
    save_positions::Tuple{Bool,Bool}=(false, true)
)
    given = ((1, x), (2, y), (3, z))
    selected = Tuple(pair for pair in given if !isnothing(pair[2]))
    if isempty(selected)
        throw(ArgumentError(
            "`OutOfBounds` needs bounds on at least one of `x`, `y` or `z`."
        ))
    end
    T = float(promote_type((typeof(v) for pair in selected for v in pair[2])...))
    for (axis, bound) in selected
        if !(bound[1] < bound[2])
            throw(ArgumentError(
                "`OutOfBounds` axis $axis has lower bound $(bound[1]) ≥ upper "*
                "bound $(bound[2])."
            ))
        end
    end
    return OutOfBounds(
        map(first, selected),
        map(pair -> (T(pair[2][1]), T(pair[2][2])), selected),
        save_positions,
    )
end


"""
    Relativistic(; fraction=0.12, save_positions=(true, false))
Terminate a particle whose kinetic energy reaches `fraction` of its rest
energy.
"""
Base.@kwdef struct Relativistic{T<:AbstractFloat} <: CallbackSpec
    fraction::T = 0.12
    save_positions::Tuple{Bool,Bool} = (true, false)
end


"""
    HighMagneticGradient(; tolerance=1e-4, save_positions=(true, false))
Terminate a particle in a region where the magnetic field varies too rapidly
for the guiding-centre approximation, i.e. where the ratio of the Larmor radius
to the characteristic field length exceeds `tolerance`.
"""
Base.@kwdef struct HighMagneticGradient{T<:AbstractFloat} <: CallbackSpec
    tolerance::T = 1e-4
    save_positions::Tuple{Bool,Bool} = (true, false)
end


# `OutOfBoundsCondition` wants the state-vector indices as a vector.
function create_callback(s::OutOfBounds, p)
    return DiscreteCallback(
        OutOfBoundsCondition(s.bounds, collect(s.u_idxs)),
        outofboundsaffect!,
        save_positions=s.save_positions,
    )
end

function create_callback(s::Relativistic, p)
    condition = if p.eom == hybridgcafo!
        RelativisticConditionHybrid
    elseif p.eom == guidingcentreapproximation!
        RelativisticConditionGCA
    elseif p.eom == lorentzforce!
        RelativisticCondition
    else
        throw(ArgumentError(
            "No relativistic termination condition for eom=$(p.eom)."
        ))
    end
    return DiscreteCallback(
        condition(; mass=p.mass, fraction=s.fraction),
        relativisticaffect!,
        save_positions=s.save_positions,
    )
end

function create_callback(s::HighMagneticGradient, p)
    return DiscreteCallback(
        MagneticGradientCondition(s.tolerance),
        magneticgradientaffect!,
        save_positions=s.save_positions,
    )
end

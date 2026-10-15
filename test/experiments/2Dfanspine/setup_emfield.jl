using Interpolations
using JLD2
using Printf
using StaticArrays

# The interpolator has a call signature `(x, z)`, returning a 6D Static vector:
# `[Ex, Ey, Ez, Bx, By, Bz]`
# What TraceParticles.jl needs is a function with call signature `(x, y, z, t)`
# that returns two 3D static vectors:
# `([Ex, Ey, Ez], [Bx, By, Bz])`
#
# Construct a wrapper that satisfies this.
#
# The parametric wrapper below also keeps the interpolator's concrete type so
# the field evaluations are type-stable.
struct EMFieldWrapper{T}
    itp::T
end
function (self::EMFieldWrapper)(x, y, z, t)
    v = self.itp(x, z)
    E = SVector(v[1], v[2], v[3])
    B = SVector(v[4], v[5], v[6])
    return E, B
end

function setup_emfield()
    #---------------------------------------------------------------------------
    # PARAMETERS
    FILE = joinpath(@__DIR__, "emfield.itp.jld2")
    #---------------------------------------------------------------------------

    # Load the metadata
    metadata = JLD2.load(FILE, "metadata")

    # Validate call signature and variables
    callsig = get(metadata, :call_signature, nothing)
    callsig == "(x, z)" ||
        error("Expected call signature itp(x, z), got $callsig")

    vars = get(metadata, :vars, nothing)
    vars == ["ex", "ey", "ez", "bx", "by", "bz"] ||
        error("Expected electromagnetic field variables, got $vars")

    # Load the interpolator
    stats = @timed JLD2.load(FILE, "itp")

    @info """EM-field loaded:
    Elapsed time: $(stats.time) seconds
    Bytes allocated: $(@sprintf "%.1E" stats.bytes)"""

    # Validate type
    itp = stats.value
    itp isa Interpolations.AbstractInterpolation ||
        error("""
itp is not an Interpolations.AbstractInterpolation. It may be a reconstruction
which throws method-errors when called. The loaded type is: $(typeof(itp))"""
        )

    # Return wrapped functor
    return EMFieldWrapper(itp)
end

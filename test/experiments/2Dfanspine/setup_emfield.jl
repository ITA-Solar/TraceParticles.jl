
using Interpolations
using JLD2
using Printf
using StaticArrays

const STATIC = false
const XZ = false
const EMFIELD_FILE = "test_t1_BE_si_cubicspline.itp.jld2"

include("toplayer.jl")

# Construct a wrapper that wraps the interpolator in a function with a
# call signature (x, y, z, t)
itp_wrapper = construct_itpwrapper(STATIC, XZ)

#=
Now load the electromagnetic field interpolation object.
This step assumes that the file `emfield_file` points to a JLD2-file containing
a vector of callable objects that each return a component of the electric
or magnetic vector field. Specifically it assumes the vector to contain the
components in the following order
    vec[1] -> bx
    vec[2] -> by
    vec[3] -> bz
    vec[4] -> ex
    vec[5] -> ey
    vec[6] -> ez
=#
time_taken = @timed fields_itp = ElectromagneticFieldInterpolator(
    [itp_wrapper(component) for component in load_object(EMFIELD_FILE)]
)
@info """EM-field loaded:
Elapsed time: $(time_taken.time) seconds
Bytes allocated: $(@sprintf "%.1E" time_taken.bytes)"""

fields_itp

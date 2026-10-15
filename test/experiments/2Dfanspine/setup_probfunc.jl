using Interpolations
using JLD2
using Printf
using TraceParticles: MHDSampler
#-------------------------------------------------------------------------------

const STATIC = false
const XZ = false
const TG_FILE = "test_t1_tg_si_cubicspline.itp.jld2"
const TARGET_DISTR_FILE = "test_t1_r_si_cubicspline.itp.jld2"
const PROPOSAL_DISTR_FILE = nothing
const XMIN = 1.5822270393371582e7
const XMAX = 1.6318359375e7
const YMIN = 1_000_000
const YMAX = 1e6
const ZMIN = -6.598462104797363e6
const ZMAX = -6.102306842803955e6
const TMIN = 0.0
const TMAX = 0

include("toplayer.jl")

# Construct a wrapper that wraps the interpolator in a function with a
# call signature (x, y, z, t)
itp_wrapper = construct_itpwrapper(STATIC, XZ)


time_taken = @timed tg = itp_wrapper(load_object(TG_FILE))
@info """Temperature-field loaded:
Elapsed time: $(time_taken.time) seconds
Bytes allocated: $(@sprintf "%.1E" time_taken.bytes)"""

time_taken = @timed target_distr = itp_wrapper(load_object(TARGET_DISTR_FILE))
@info """Target distribution loaded:
Elapsed time: $(time_taken.time) seconds
Bytes allocated: $(@sprintf "%.1E" time_taken.bytes)"""

if isnothing(PROPOSAL_DISTR_FILE) || PROPOSAL_DISTR_FILE == ""
    proposal_distr = target_distr
else
    time_taken = @timed proposal_distr = itp_wrapper(load_object(PROPOSAL_DISTR_FILE))
    @info """Proposal distribution loaded:
    Elapsed time: $(time_taken.time) seconds
    Bytes allocated: $(@sprintf "%.1E" time_taken.bytes)"""
end

MHDSampler(;
    temperature=tg,
    target_distr=target_distr,
    xmin=XMIN,
    xmax=XMAX,
    ymin=YMIN,
    ymax=YMAX,
    zmin=ZMIN,
    zmax=ZMAX,
    tmin=TMIN,
    tmax=TMAX,
    proposal_distr=proposal_distr,
    maxval=maximum(proposal_distr)
)

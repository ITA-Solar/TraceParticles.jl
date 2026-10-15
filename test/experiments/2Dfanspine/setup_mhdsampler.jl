using Interpolations
using JLD2
using Printf
using TraceParticles: MHDSampler

"""
    StaticXZInterpolation{T}
        itp::T
An interpolator that calls `(x, z)` if called with `(x, y, z, t)`.
"""
struct StaticXZInterpolation{T}
    itp::T
end
function (self::StaticXZInterpolation)(x::Real, y::Real, z::Real, t::Real)
    self.itp(x, z)
end

function setup_mhdsampler()
    #---------------------------------------------------------------------------
    # PARAMETERS
    TG_FILE = joinpath(@__DIR__, "temperature.itp.jld2")
    TARGET_DISTR_FILE = joinpath(@__DIR__, "density.itp.jld2")
    PROPOSAL_DISTR_FILE = nothing
    XMIN = 1.5822270393371582e7
    XMAX = 1.6318359375e7
    YMIN = 1e6
    YMAX = 1e6
    ZMIN = -6.598462104797363e6
    ZMAX = -6.102306842803955e6
    TMIN = 0.0
    TMAX = 0.0
    #---------------------------------------------------------------------------

    tg = StaticXZInterpolation(JLD2.load(TG_FILE, "itp"))
    target_distr = StaticXZInterpolation(JLD2.load(TARGET_DISTR_FILE, "itp"))
    if isnothing(PROPOSAL_DISTR_FILE) || PROPOSAL_DISTR_FILE == ""
        proposal_distr = target_distr
    else
        proposal_distr = StaticXZInterpolation(
            JLD2.load(PROPOSAL_DISTR_FILE,"itp")
        )
    end
    return MHDSampler(
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
        maxval=maximum(proposal_distr.itp.itp.itp.coefs)
    )
end


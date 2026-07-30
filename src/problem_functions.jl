# Created 08.02.24
# Author: e.s.oyre@astro.uio.no
#
#           problem_functions.jl
#-------------------------------------------------------------------------------
# Contains problem functions to be used as the keyword-argument `prob_func` in 
# an `EnsembleProblem` from SciMl.
#-------------------------------------------------------------------------------

"""
    struct PredefinedICs
        u0
        mu0
        tspan
        params
# Methods
    (::PredefinedICs)(prob, ctx)
Remakes the problem with predefined initial conditions `u0[ctx.sim_id]`, time
span `tspan[ctx.sim_id]`, and parameters `params[ctx.sim_id]`.
"""
struct PredefinedICs{
    T1<:AbstractVector,
    T2<:AbstractVector,
    T3<:AbstractVector,
}
    u0::T1
    tspan::T2
    params::T3
    function PredefinedICs(
        u0::T1,
        tspan::T2,
        params::T3
    ) where {
        T1<:AbstractVector,
        T2<:AbstractVector,
        T3<:AbstractVector,
    }
        if length(u0) == length(tspan) == length(params)
            new{T1,T2,T3}(u0, tspan, params)
        else
            error("u0, tspan and params must have the same length.")
        end
    end
end
function (self::PredefinedICs)(prob, ctx)
    remake(
        prob,
        u0=self.u0[ctx.sim_id],
        tspan=self.tspan[ctx.sim_id],
        p=self.params[ctx.sim_id],
    )
end

"""
    mhdsampling(
        mass::Real,
        rng::AbstractRNG,
        proposal_distr,
        target_distr,
        tg_itp,
        xbounds::Tuple{Real,Real},
        ybounds::Tuple{Real,Real},
        zbounds::Tuple{Real,Real},
        temporalbounds::Tuple{Real,Real},
        max_value::Real,
    )
Draw `x`, `y`, `z` positions and time `t` from a `proposal_distr`ibution using
rejection sampling. Draws the 3D velocity vector from a Maxwellian distribution
with a temperature given by `tg_itp(x,y,z)`. Evaluates the statistical
importance weight of the particle using the `target_distr`ibution.
"""
function mhdsampling(
    mass::Real,
    rng::AbstractRNG,
    proposal_distr,
    target_distr,
    tg_itp,
    xbounds::T,
    ybounds::T,
    zbounds::T,
    temporalbounds::T,
    max_value::Real,
) where {T<:Tuple{Real,Real}}
    # Extract spatial and temporal limits
    x0, xf = xbounds
    y0, yf = ybounds
    z0, zf = zbounds
    tmin, tmax = temporalbounds
    # Sample osition and time from the proposal distrubution using rejection
    # sampling.
    x, y, z, t, nrejections = rejectionsample(
        rng,
        proposal_distr,
        max_value,
        x0, xf,
        y0, yf,
        z0, zf,
        tmin, tmax,
    )
    weight = eltype(mass)(
        target_distr(x, y, z, t) / proposal_distr(x, y, z, t)
    )
    # Sample velocity components from a Maxwell-Boltzmann distribution.
    temperature = tg_itp(x, y, z, t)
    vx = maxwellianvelocitysample(rng, temperature, mass)
    vy = maxwellianvelocitysample(rng, temperature, mass)
    vz = maxwellianvelocitysample(rng, temperature, mass)
    return x, y, z, vx, vy, vz, t, weight, nrejections
end

"""
    getu0_guidingcentre!(
        u0::AbstractVector,
        x::Real,
        y::Real,
        z::Real,
        vx::Real,
        vy::Real,
        vz::Real,
        t::Real,
        mass::Real,
        charge::Real,
        electromagneticfield,
    )
Calculate the guiding centre position and parallel velocity from the full
orbit coordinates and store them in `u0`. Returns the magnetic moment.
"""
function getu0_guidingcentre!(
    u0::AbstractVector,
    x::Real,
    y::Real,
    z::Real,
    vx::Real,
    vy::Real,
    vz::Real,
    t::Real,
    mass::Real,
    charge::Real,
    electromagneticfield,
)
    efield_at_pos, bfield_at_pos = electromagneticfield(x, y, z, t)
    R, vparal, magneticmoment = get_guidingcentre(
        SVector(x, y, z),
        SVector(vx, vy, vz),
        bfield_at_pos,
        efield_at_pos,
        charge,
        mass
    )
    u0[1:4] .= [R[1], R[2], R[3], vparal]
    return magneticmoment
end

"""
    getu0_fullorbit!(
        u0::AbstractVector,
        x::Real,
        y::Real,
        z::Real,
        vx::Real,
        vy::Real,
        vz::Real,
        t::Real,
        mass::Real,
        charge::Real,
        electromagneticfield,
    )
Store the full orbit position and velocity in `u0`. Returns 0.0 as a dummy-
magnetic moment.
"""
function getu0_fullorbit!(
    u0::AbstractVector,
    x::Real,
    y::Real,
    z::Real,
    vx::Real,
    vy::Real,
    vz::Real,
    _::Real,
    _::Real,
    _::Real,
    _,
)
    u0 .= [x, y, z, vx, vy, vz]
    return 0.0
end

"""
    struct MHDSamplingGCA
        proposal_distr
        target_distr
        tg_itp
        xbounds
        ybounds
        zbounds
        t0bounds
        tf
        max_value

# Methods
    (::MHDSamplingGCA)(prob, _)
Draw `x`, `y`, `z` positions and time `t` from a `proposal_distr`ibution using
rejection sampling. Draws the 3D velocity vector from a Maxwellian distribution
with a temperature given by `tg_itp(x,y,z)`. Evaluates the importance weight of
the particle using the `target_distr`ibution. Calculates the corresponding
guiding centre, magnetic moment, and parallel velocity and remakes the problem
with `GCAParams`.
"""
struct MHDSamplingGCA{
    T1,T2,T3,T4<:Tuple{Real,Real},T5<:Real,T6<:Real,
}
    proposal_distr::T1
    target_distr::T2
    tg_itp::T3
    xbounds::T4
    ybounds::T4
    zbounds::T4
    t0bounds::T4
    tf::T5
    max_value::T6
end
function (self::MHDSamplingGCA)(prob, ctx)
    x, y, z, vx, vy, vz, t, weight, nrejections = mhdsampling(
        prob.p.mass,
        ctx.rng,
        self.proposal_distr,
        self.target_distr,
        self.tg_itp,
        self.xbounds,
        self.ybounds,
        self.zbounds,
        self.t0bounds,
        self.max_value,
    )
    u0 = Vector{Float64}(undef, 4)
    magneticmoment = getu0_guidingcentre!(
        u0, x, y, z, vx, vy, vz, t,
        prob.p.mass, prob.p.charge, prob.p.electromagneticfield
    )
    #=
    # This is not thread-safe if `safetycopy=false`
    # because I'm modifying the arguments of the problem function. This leads
    # to race conditions with multiple threads. However, this mutating
    # could work on distributed systems without threads
    @. prob.u0 = [R[1], R[2], R[3], vparal]
    prob.p.particledata.terminationcode = TerminationCode.NotTerminated
    prob.p.particledata.weight = weight
    prob.p.particledata.nrejections = nrejections
    return prob
    =#

    # This remaking allocates memory when creating the new params, which
    # still points to the same actual values, but it is thread-safe because
    # as long as the parameters that points to `prob.p` is not mutated.
    return remake(
        prob;
        u0=u0,
        tspan=(t, self.tf),
        p=GCAParams(
            charge=prob.p.charge,
            mass=prob.p.mass,
            electromagneticfield=prob.p.electromagneticfield,
            magneticmoment=magneticmoment,
            weight=weight,
            nrejections=nrejections,
        )
    )
end

"""
    struct MHDSamplingHybrid
        proposal_distr
        target_distr
        tg_itp
        xbounds
        ybounds
        zbounds
        t0bounds
        tf
        max_value
        initialeomid

# Methods
    (::MHDSamplingHybrid)(prob, _)
Draw `x`, `y`, `z` positions and time `t` from a `proposal_distr`ibution using
rejection sampling. Draws the 3D velocity vector from a Maxwellian distribution
with a temperature given by `tg_itp(x,y,z)`. Evaluates the importance weight of
the particle using the `target_distr`ibution. Remakes the problem with
`HybridParams`.
"""
struct MHDSamplingHybrid{
    T1,T2,T3,T4<:Tuple{Real,Real},T5<:Real,T6<:Real,T7<:Int
}
    proposal_distr::T1
    target_distr::T2
    tg_itp::T3
    xbounds::T4
    ybounds::T4
    zbounds::T4
    t0bounds::T4
    tf::T5
    max_value::T6
    initialeomid::T7
end
function (self::MHDSamplingHybrid)(prob, ctx)
    x, y, z, vx, vy, vz, t, weight, nrejections = mhdsampling(
        prob.p.mass,
        ctx.rng,
        self.proposal_distr,
        self.target_distr,
        self.tg_itp,
        self.xbounds,
        self.ybounds,
        self.zbounds,
        self.t0bounds,
        self.max_value,
    )
    u0 = zeros(Float64, 6)
    if self.initialeomid == 1
        magneticmoment = getu0_guidingcentre!(
            u0, x, y, z, vx, vy, vz, t,
            prob.p.mass, prob.p.charge, prob.p.electromagneticfield
        )
    elseif self.initialeomid == 2
        magneticmoment = getu0_fullorbit!(
            u0, x, y, z, vx, vy, vz, t,
            prob.p.mass, prob.p.charge, prob.p.electromagneticfield
        )
    else
        error("Invalid initialeomid: $(self.initialeomid). Must be 1 or 2.")
    end
    # This remaking allocates memory when creating the new params, which
    # still points to the same actual values, but it is thread-safe because
    # as long as the parameters that points to `prob.p` is not mutated.
    return remake(
        prob;
        u0=u0,
        tspan=(t, self.tf),
        p=HybridParams(
            charge=prob.p.charge,
            mass=prob.p.mass,
            electromagneticfield=prob.p.electromagneticfield,
            magneticmoment=magneticmoment,
            weight=weight,
            nrejections=nrejections,
            eomid=self.initialeomid,
            rng=ctx.rng
        )
    )
end

"""
    initialconditions_mhdsampling(
        npart::Int,
        electromagneticfield,
        EoM::Function,
        mass::Real,
        charge::Real,
        args...;
        initialeomid::Int=1
    )
Generate initial conditions for `npart` particles using `mhdsampling`.
Determines whether to calculate the initial conditions as a guiding centre or
full orbit based on the provided `EoM` function and `initialeomid`.
"""
function initialconditions_mhdsampling(
    npart::Int,
    electromagneticfield,
    EoM::Function,
    mass::Real,
    charge::Real,
    args...;
    initialeomid::Int=1
)
    u0s = [zeros(6) for _ in 1:npart]
    magneticmoments = Vector{Float64}(undef, npart)
    t0s = Vector{Float64}(undef, npart)
    weights = Vector{Float64}(undef, npart)
    nrejections = Vector{Int}(undef, npart)

    if (
        EoM == lorentzforce! ||
        (EoM == hybridgcafo! && initialeomid == 2)
    )
        u0func! = getu0_fullorbit!
    elseif (
        EoM == guidingcentreapproximation! ||
        (EoM == hybridgcafo! && initialeomid == 1)
    )
        u0func! = getu0_guidingcentre!
    end
    Threads.@threads for i in 1:npart
        x, y, z, vx, vy, vz, t, weights[i], nrejections[i] = mhdsampling(
            mass, args...
        )
        magneticmoments[i] = u0func!(
            u0s[i], x, y, z, vx, vy, vz, t,
            mass, charge, electromagneticfield
        )
        t0s[i] = t
    end
    return u0s, magneticmoments, t0s, weights, nrejections
end

# Created 08.02.24
# Author: e.s.oyre@astro.uio.no
#
#           problem_functions.jl
#-------------------------------------------------------------------------------
# Contains problem functions to be used as the keyword-argument `prob_func` in 
# an DifferentialEquations.EnsembleProblem.
#-------------------------------------------------------------------------------

"""
    struct PredefinedICs
        u0
        mu0
        tspan
        paramstruct
        kwargs...
# Methods
    (::PredefinedICs)(prob, i, _)
Remakes the problem with predefined initial conditions `u0[i]`, magnetic
moment `mu0[i]`, and time span `tspan[i]`. The `paramstruct` is called with
the keyword arguments in `kwargs...` to create the parameter struct for the
new problem.
"""
struct PredefinedICs{T1<:AbstractVector,T2<:AbstractVector,T3<:AbstractVector,T4,T5}
    u0::T1
    mu0::T2
    tspan::T3
    paramstruct::T4
    kwargs::T5

    function PredefinedICs(
        u0::T1,
        mu0::T2,
        tspan::T3,
        paramstruct::T4;
        kwargs...
    ) where {T1<:AbstractVector,T2<:AbstractVector,T3<:AbstractVector,T4}
        if length(u0) == length(mu0) == length(tspan)
            new{T1,T2,T3,T4,typeof(kwargs)}(u0, mu0, tspan, paramstruct, kwargs)
        else
            error("u0, mu0, and tspan must have the same length.")
        end
    end
end
function (self::PredefinedICs)(prob, i, _)
    remake(
        prob,
        f=prob.f,
        u0=self.u0[i],
        tspan=self.tspan[i],
        p=self.paramstruct(
            charge=prob.p.charge,
            mass=prob.p.mass,
            magneticmoment=self.mu0[i],
            electromagneticfield=prob.p.electromagneticfield,
            terminationcode=TerminationCode.NotTerminated,
            self.kwargs...
        )
    )
end

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
        rng
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
    (::MHDSamplingGCA)(prob, _, _)
Draw `x`, `y`, and `z` positions from a given `proposal_distr`ibution using
rejection sampling. Draws the 3D velocity vector from a Maxwellian distribution
with a temperature given by `tg_itp(x,y,z)`. Evaluates the statistical weight of
the particle using the `target_distr`ibution. Calculates the corresponding
guiding centre, magnetic moment, and parallel velocity and remakes the problem
with `GCAParams`.
"""
struct MHDSamplingGCA{
    T1<:AbstractRNG,T2,T3,T4,T5<:Tuple{Real,Real},T6<:Real,T7<:Real,
}
    rng::T1
    proposal_distr::T2
    target_distr::T3
    tg_itp::T4
    xbounds::T5
    ybounds::T5
    zbounds::T5
    t0bounds::T5
    tf::T6
    max_value::T7
end
function (self::MHDSamplingGCA)(prob, _, _)
    x, y, z, vx, vy, vz, t, weight, nrejections = mhdsampling(
        prob.p.mass,
        self.rng,
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
            terminationcode=TerminationCode.NotTerminated
        )
    )
end

"""
    struct MHDSamplingHybrid
        rng
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
    (::MHDSamplingHybrid)(prob, _, _)
Draw `x`, `y`, and `z` positions from a given `proposal_distr`ibution using
rejection sampling. Draws the 3D velocity vector from a Maxwellian distribution
with a temperature given by `tg_itp(x,y,z)`. Evaluates the statistical weight of
the particle using the `target_distr`ibution. Remakes the problem with
`HybridParams`.
"""
struct MHDSamplingHybrid{
    T1<:AbstractRNG,T2,T3,T4,T5<:Tuple{Real,Real},T6<:Real,T7<:Real,T8<:Int
}
    rng::T1
    proposal_distr::T2
    target_distr::T3
    tg_itp::T4
    xbounds::T5
    ybounds::T5
    zbounds::T5
    t0bounds::T5
    tf::T6
    max_value::T7
    initialeomid::T8
end
function (self::MHDSamplingHybrid)(prob, _, _)
    x, y, z, vx, vy, vz, t, weight, nrejections = mhdsampling(
        prob.p.mass,
        self.rng,
        self.proposal_distr,
        self.target_distr,
        self.tg_itp,
        self.xbounds,
        self.ybounds,
        self.zbounds,
        self.t0bounds,
        self.max_value,
    )
    u0 = Vector{Float64}(undef, 6)
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
            terminationcode=TerminationCode.NotTerminated,
            initialeomid=self.initialeomid # initialised with the EoM `lorentzforce!`
        )
    )
end

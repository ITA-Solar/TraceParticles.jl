# Created 08.02.24
# Author: e.s.oyre@astro.uio.no
#
#           problem_functions.jl
#-------------------------------------------------------------------------------
# Contains problem functions to be used as the keyword-argument `prob_func` in 
# an DifferentialEquations.EnsembleProblem.
#-------------------------------------------------------------------------------

"""
Draw initial conditions from an array.
"""
struct PposPvel
    u0::Vector{<:Vector{<:Real}}
    mu0::Vector{<:Real}
    tspan::Vector{<:Tuple{Real,Real}}
    paramstruct::Any
end
function (self::PposPvel)(prob, i, repeat)
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
            terminationcode=TerminationCode.NotTerminated
        )
    )
end

"""
    struct DposMBvelGCA
        rng
        proposal_distr
        target_distr
        tg_itp
        domain
        max_value

# Methods
    (::DposMBvelGCA)(prob, _, _)
Draw `x`, `y`, and `z` positions from a given `proposal_distr`ibution using
rejection sampling. Draws the 3D velocity vector from a Maxwellian distribution
with a temperature given by `tg_itp(x,y,z)`. Evaluates the statistical weight of
the particle using the `target_distr`ibution. Calculates the corresponding
guiding centre, magnetic moment, and parallel velocity and remakes the problem
with `GCAParams`.
"""
struct DposMBvelGCA{T1<:AbstractRNG,T2,T3,T4<:AbstractVector,T5<:Number}
    rng::T1
    proposal_distr::T2
    target_distr::T3
    tg_itp::T3
    domain::T4
    max_value::T5
    time::T5
end
function (self::DposMBvelGCA)(prob, _, _)
    if length(self.domain) == 2
        x, z, nrejections = rejectionsample(
            self.rng,
            self.proposal_distr,
            self.max_value,
            self.domain[1][1],
            self.domain[1][2],
            self.domain[2][1],
            self.domain[2][2],
        )
        y = 1e6 # !!! Hardcoded y-value
    else
        x, y, z, nrejections = rejectionsample(
            self.rng,
            self.proposal_distr,
            self.max_value,
            self.domain[1][1],
            self.domain[1][2],
            self.domain[2][1],
            self.domain[2][2],
            self.domain[3][1],
            self.domain[3][2],
        )
    end
    weight = self.target_distr(x, y, z) / self.proposal_distr(x, y, z)

    # Velocity
    charge = prob.p.charge
    mass = prob.p.mass
    temperature = self.tg_itp(x, y, z)
    vel = @SVector [
        maxwellianvelocitysample(self.rng, temperature, mass) for _ in 1:3
    ]
    efield_at_pos, bfield_at_pos = prob.p.electromagneticfield(
        x, y, z, self.time
    )
    R, vparal, magneticmoment = get_guidingcentre(
        SVector(x, y, z),
        vel,
        bfield_at_pos,
        efield_at_pos,
        charge,
        mass
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
        u0=[R[1], R[2], R[3], vparal],
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
    struct DposMBvelHybrid
        rng
        proposal_distr
        target_distr
        tg_itp
        domain
        max_value

# Methods
    (::DposMBvelHybrid)(prob, _, _)
Draw `x`, `y`, and `z` positions from a given `proposal_distr`ibution using
rejection sampling. Draws the 3D velocity vector from a Maxwellian distribution
with a temperature given by `tg_itp(x,y,z)`. Evaluates the statistical weight of
the particle using the `target_distr`ibution. Remakes the problem with
`HybridParams`.
"""
struct DposMBvelHybrid{
    T1<:AbstractRNG,T2,T3,T4<:Tuple{Real,Real},T5<:Number,T6<:Int
}
    rng::T1
    proposal_distr::T2
    target_distr::T3
    tg_itp::T3
    xbounds::T4
    ybounds::T4
    zbounds::T4
    max_value::T5
    time::T5
    initialeomid::T6
end
function (self::DposMBvelHybrid)(prob, _, _)
    x, y, z, nrejections = rejectionsample(
        self.rng,
        self.proposal_distr,
        self.max_value,
        self.xbounds[1],
        self.xbounds[2],
        self.ybounds[1],
        self.ybounds[2],
        self.zbounds[1],
        self.zbounds[2],
    )
    weight = eltype(prob.p.mass)(
        self.target_distr(x, y, z) / self.proposal_distr(x, y, z)
    )

    # Velocity
    temperature = self.tg_itp(x, y, z)
    vel = @SVector [
        maxwellianvelocitysample(self.rng, temperature, prob.p.mass)
        for _ in 1:3
    ]
    if self.initialeomid == 1
        efield_at_pos, bfield_at_pos = prob.p.electromagneticfield(
            x, y, z, self.time
        )
        R, vparal, magneticmoment = get_guidingcentre(
            SVector(x, y, z),
            vel,
            bfield_at_pos,
            efield_at_pos,
            prob.p.charge,
            prob.p.mass
        )
        u0 = [R[1], R[2], R[3], vparal, 0.0, 0.0]
    elseif self.initialeomid == 2
        magneticmoment = 0.0
        u0 = [x, y, z, vel[1], vel[2], vel[3]]
    else
        error("Invalid initialeomid: $(self.initialeomid). Must be 1 or 2.")
    end
    # This remaking allocates memory when creating the new params, which
    # still points to the same actual values, but it is thread-safe because
    # as long as the parameters that points to `prob.p` is not mutated.
    return remake(
        prob;
        u0=u0,
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

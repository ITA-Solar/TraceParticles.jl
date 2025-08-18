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
    paramstruct::DataType
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
Draws x, y, and z positions from a given distribution using rejection sampling.
"""
struct DposMBvel{T1<:AbstractRNG,T2,T3,T4<:AbstractVector,T5<:Number}
    rng::T1
    proposal_distr::T2
    target_distr::T3
    tg_itp::T3
    domain::T4
    max_value::T5
end
function (self::DposMBvel)(prob, i, repeat)
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
    efield_at_pos, bfield_at_pos = prob.p.electromagneticfield(x, y, z)
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
    hybrid = prob.f == guidingcentreapproximation!
    paramstruct, u0 = hybrid ?
        (GCAParams, [R[1], R[2], R[3], vparal]) :
        (HybridParams, [R[1], R[2], R[3], vparal, 0.0, 0.0])
    return remake(
        prob;
        u0=u0,
        p=paramstruct(
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

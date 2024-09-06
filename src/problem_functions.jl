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
    tspan::Vector{<:Tuple{Real, Real}}
end
function (self::PposPvel)(prob, i, repeat)
    remake(
        prob,
        f = prob.f,
        u0 = self.u0[i],
        tspan = self.tspan[i],
        p = (
            charge = prob.p.charge,
            mass = prob.p.mass,
            magneticmoment = self.mu0[i],
            fields = prob.p.fields,
            weight = 1.0,
            userdata = UserData("none")
        )
    )
end



"""
Draw uniform position, Maxwell-Boltzmann distributed velocity.
"""
struct UposMBvel{RealT}
    rng    ::AbstractRNG
    xbounds::Tuple{RealT, RealT}
    ybounds::Tuple{RealT, RealT}
    zbounds::Tuple{RealT, RealT}
    tg_itp::AbstractInterpolation
end 
function (self::UposMBvel{<:Real})(prob, i, repeat)
    wp_part = eltype(self.xbounds)
    charge = prob.p.charge
    mass = prob.p.mass
    x0 = rand(self.rng, wp_part, self.xbounds...)
    y0 = rand(self.rng, wp_part, self.ybounds...)
    z0 = rand(self.rng, wp_part, self.zbounds...)
    vel = maxwellianvelocitysample(self.rng, self.tg_itp, mass, x0, z0)
    # Has type Vector{Vector}
    bfield_at_pos = [prob.p.fields[i](x0, z0) for i in 1:3]
    efield_at_pos = [prob.p.fields[i](x0, z0) for i in 4:6]
    R, vparal, magneticmoment = get_guidingcentre(
        [x0, y0, z0],
        vel,
        bfield_at_pos,
        efield_at_pos,
        charge,
        mass
        )
    remake(prob, u0=[R; vparal],
    p=(
        charge = charge,
        mass = mass, 
        magneticmoment = magneticmoment,
        fields = prob.p.fields,
        weight = 1.0,
        userdata = UserData("none")
        )
    )
end


"""
Draws x, y, and z positions from a given distribution using the rejection
sampling method.
"""
struct DposMBvel
    rng::AbstractRNG
    proposal_distr::Any
    target_distr::Any
    domain::Vector{Tuple{<:Real, <:Real}}
    max_value::Real
    tg_itp::AbstractInterpolation
end
function (self::DposMBvel)(prob, i, repeat)
    (Rx, Rz), nrejections = rejectionsample(
        self.proposal_distr,
        self.max_value,
        self.domain,
        self.rng
        )
    acceptanceratio = 1/(nrejections + 1)
    weight = self.target_distr(Rx, Rz) / self.proposal_distr(Rx, Rz)
    R = [Rx, 1e6, Rz] # !!! Hardcoded y-value

    # Velocity
    charge = prob.p.charge
    mass = prob.p.mass
    vel = maxwellianvelocitysample(self.rng, self.tg_itp, mass, Rx, Rz)
    bfield_at_pos = [prob.p.fields[i](Rx, Rz) for i in 1:3]
    efield_at_pos = [prob.p.fields[i](Rx, Rz) for i in 4:6]
    vparal, magneticmoment, _ = get_gca_velocities(
        vel,
        bfield_at_pos,
        efield_at_pos,
        mass
        )
    remake(prob, u0=[R; vparal],
    p=(
        charge = charge,
        mass = mass, 
        magneticmoment = magneticmoment,
        fields = prob.p.fields,
        weight = weight,
        acceptanceratio = acceptanceratio,
        userdata = UserData("none")
        )
    )
end 


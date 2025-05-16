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
end
function (self::PposPvel)(prob, i, repeat)
    remake(
        prob,
        f=prob.f,
        u0=self.u0[i],
        tspan=self.tspan[i],
        p=(
            charge=prob.p.charge,
            mass=prob.p.mass,
            magneticmoment=self.mu0[i],
            electromagneticfield=prob.p.electromagneticfield,
            weight=1.0,
            userdata=UserData("None")
        )
    )
end



"""
Draw uniform position, Maxwell-Boltzmann distributed velocity.
"""
struct UposMBvel
    rng::AbstractRNG
    xbounds::Tuple{<:Real,<:Real}
    ybounds::Tuple{<:Real,<:Real}
    zbounds::Tuple{<:Real,<:Real}
    tg_itp::AbstractInterpolation
end
function (self::UposMBvel)(prob, i, repeat)
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
            charge=charge,
            mass=mass,
            magneticmoment=magneticmoment,
            fields=prob.p.fields,
            weight=1.0,
            userdata=UserData("None")
        )
    )
end


"""
Draws x, y, and z positions from a given distribution using rejection sampling.
"""
struct DposMBvel
    rng::AbstractRNG
    proposal_distr::Any
    target_distr::Any
    domain::Vector{Tuple{<:Real,<:Real}}
    max_value::Real
    tg_itp::Any
end
function (self::DposMBvel)(prob, i, repeat)
    pos, nrejections = rejectionsample(
        self.proposal_distr,
        self.max_value,
        self.domain,
        self.rng
    )
    if length(pos) == 2
        pos = [pos[1], 1e6, pos[2]] # !!! Hardcoded y-value
    end
    acceptanceratio = 1 / (nrejections + 1)
    weight = self.target_distr(pos...) / self.proposal_distr(pos...)

    # Velocity
    charge = prob.p.charge
    mass = prob.p.mass
    vel = maxwellianvelocitysample(self.rng, self.tg_itp, mass, pos...)
    efield_at_pos, bfield_at_pos = prob.p.electromagneticfield(pos...)
    u = zeros(Float64, 4)
    magneticmoment = get_guidingcentre!(
        u,
        pos,
        vel,
        bfield_at_pos,
        efield_at_pos,
        charge,
        mass
    )
    vparal = u[4]
    remake(prob, u0=[pos; vparal],
        p=(
            charge=charge,
            mass=mass,
            magneticmoment=magneticmoment,
            electromagneticfield=prob.p.electromagneticfield,
            weight=weight,
            acceptanceratio=acceptanceratio,
            userdata=UserData("None")
        )
    )
end


"""
Draws x and z positions from a given distribution using rejection sampling.
"""
struct DposMBvel_2Dxz
    rng::AbstractRNG
    proposal_distr::Any
    target_distr::Any
    domain::Vector{Tuple{<:Real,<:Real}}
    max_value::Real
    tg_itp::Any
end
function (self::DposMBvel_2Dxz)(prob, i, repeat)
    (Rx, Rz), nrejections = rejectionsample(
        self.proposal_distr,
        self.max_value,
        self.domain,
        self.rng
    )
    acceptanceratio = 1 / (nrejections + 1)
    weight = self.target_distr(Rx, Rz) / self.proposal_distr(Rx, Rz)
    R = [Rx, 1e6, Rz] # !!! Hardcoded y-value

    # Velocity
    charge = prob.p.charge
    mass = prob.p.mass
    vel = maxwellianvelocitysample(self.rng, self.tg_itp, mass, Rx, Rz)
    efield_at_pos, bfield_at_pos = prob.p.electromagneticfield(Rx, Rz)
    u = zeros(Float64, 4)
    magneticmoment = get_guidingcentre!(
        u,
        R,
        vel,
        bfield_at_pos,
        efield_at_pos,
        charge,
        mass
    )
    vparal = u[4]

    remake(prob, u0=[R; vparal],
        p=(
            charge=charge,
            mass=mass,
            magneticmoment=magneticmoment,
            electromagneticfield=prob.p.electromagneticfield,
            weight=weight,
            acceptanceratio=acceptanceratio,
            userdata=UserData("None")
        )
    )
end

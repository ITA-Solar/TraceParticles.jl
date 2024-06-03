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

"""
Draw uniform position, Maxwell-Boltzmann distributed velocity.
"""
struct UposMBvel{RealT}
    rng    ::AbstractRNG
    xbounds::Tuple{RealT, RealT}
    ybounds::Tuple{RealT, RealT}
    zbounds::Tuple{RealT, RealT}
    temp_itp::AbstractInterpolation
end 
function (self::UposMBvel{<:Real})(prob, i, repeat)
    print("Particle $i \r")
    wp_part = eltype(self.xbounds)
    charge = prob.p.charge
    mass = prob.p.mass
    x0 = rand(self.rng, wp_part, self.xbounds...)
    y0 = rand(self.rng, wp_part, self.ybounds...)
    z0 = rand(self.rng, wp_part, self.zbounds...)
    vel = maxwellianvelocitysample(self.rng, self.temp_itp, mass, x0, z0)
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
        userdata = UserData("none")
        )
    )
end


"""
Draws x, y, and z positions from a given distribution using the rejection
sampling method.
"""
struct DposMBvel{RealT}
    rng::AbstractRNG
    proposal_distr::Function
    target_distr::Function
    domain::Vector{Tuple{RealT, RealT}}
    max_value::RealT
    temp_itp::AbstractInterpolation
end
function (self::DposMBvel)(prob, i, repeat)
    (Rx, Rz), acceptanceratio = rejectionsample(
        self.proposal_distr,
        self.max_value,
        self.domain,
        self.rng
        )
    weight = self.target_distr(Rx, Rz) / self.proposal_distr(Rx, Rz)
    R = [Rx, 1e6, Rz]
    
    # Velocity
    charge = prob.p.charge
    mass = prob.p.mass
    vel = maxwellianvelocitysample(self.rng, self.temp_itp, mass, Rx, Rz)
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
        temp = prob.p.temp,
        userdata = UserData("none")
        )
    )
end 


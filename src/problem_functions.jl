# Created 08.02.24
# Author: e.s.oyre@astro.uio.no
#
#           problem_functions.jl
#-------------------------------------------------------------------------------
# Contains problem functions to be used as the keyword-argument `prob_func` in
# an `EnsembleProblem` from SciML.
#
# Also contains functions that intitialises problem functions from
# `TraceParticlesParameters`. These functions are called by the
# `TraceParticlesProblem` constructor.
#
# Adding a new way of generating initial conditions to `TraceParticleProblem`
# means adding a specification struct and one `init_probfunc` method.
#
# Available problem functors:
# - PredefinedICs
#-------------------------------------------------------------------------------

"""
    struct PredefinedICs
        u0
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

#-------------------------------------------------------------------------------
# Problem functions built on a sampler
# One per equation of motion.
#-------------------------------------------------------------------------------

"""
    SampleFullOrbit(sampler, tf)

An `EnsembleProblem` problem function for full-orbit particles whose initial
state vector is drawn from `sampler`, integrated until `tf`.
"""
struct SampleFullOrbit{T,TF<:Real}
    sampler::T
    tf::TF
end
function (self::SampleFullOrbit)(prob, ctx)
    x, y, z, vx, vy, vz, t, weight, nrejections =
        self.sampler(prob.p.mass, ctx.rng)
    u0 = [x, y, z, vx, vy, vz]
    params = FullOrbitParams(
        charge=prob.p.charge,
        mass=prob.p.mass,
        electromagneticfield=prob.p.electromagneticfield,
        weight=weight,
        nrejections=nrejections,
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
        p=params
    )
end

"""
    SampleGCA(sampler, tf)

An `EnsembleProblem` problem function for guiding-centre particles whose
initial state vector is drawn from `sampler`, integrated until `tf`.
"""
struct SampleGCA{T,TF<:Real}
    sampler::T
    tf::TF
end
function (self::SampleGCA)(prob, ctx)
    x, y, z, vx, vy, vz, t, weight, nrejections =
        self.sampler(prob.p.mass, ctx.rng)
    u0 = Vector{typeof(x)}(undef, 4)
    magneticmoment = getu0_guidingcentre!(
        u0, x, y, z, vx, vy, vz, t,
        prob.p.mass, prob.p.charge, prob.p.electromagneticfield
    )
    params = GCAParams(
        charge=prob.p.charge,
        mass=prob.p.mass,
        electromagneticfield=prob.p.electromagneticfield,
        magneticmoment=magneticmoment,
        weight=weight,
        nrejections=nrejections,
    )
    return remake(
        prob;
        u0=u0,
        tspan=(t, self.tf),
        p=params
    )
end


"""
    SampleHybrid(sampler, tf, initialeomid, hybridscheme)

An `EnsembleProblem` problem function for hybrid GCA/full-orbit particles whose
initial state vector is drawn from `sampler`, integrated until `tf`.

`initialeomid` selects which equations of motion a particle starts on, and
`hybridscheme` whether the switches are recorded (see
[`HybridParamsWithDetection`](@ref)) or not.
"""
struct SampleHybrid{T,TF<:Real,T2<:EoMID.T,T3<:HybridScheme.T}
    sampler::T
    tf::TF
    initialeomid::T2
    hybridscheme::T3
end
function (self::SampleHybrid)(prob, ctx)
    x, y, z, vx, vy, vz, t, weight, nrejections =
        self.sampler(prob.p.mass, ctx.rng)
    # Hybrid scheme with GCA equations initially will only set the
    # first 4 variables of the state vector, possibly leaving the
    # fifth and sixth variable NaN, which the DiffEq-solver does not
    # like in its state vector, hence we initialise the state-vector to
    # zero in this case.
    u0 = zeros(typeof(x), 6)

    if self.initialeomid == EoMID.FullOrbit
        u0 .= [x, y, z, vx, vy, vz]
        #= `HybridParams` requires the magnetic moment to have the same type as
        the charge and the mass, so this placeholder must not be a bare
        `Float64` `NaN`. =#
        magneticmoment = oftype(prob.p.mass, NaN)
    elseif self.initialeomid == EoMID.GuidingCentreApproximation
        magneticmoment = getu0_guidingcentre!(
            u0, x, y, z, vx, vy, vz, t,
            prob.p.mass, prob.p.charge, prob.p.electromagneticfield
        )
    else
        error("Invalid initialeomid: $(self.initialeomid).")
    end

    if self.hybridscheme == HybridScheme.Default
        params = HybridParams(
            charge=prob.p.charge,
            mass=prob.p.mass,
            electromagneticfield=prob.p.electromagneticfield,
            magneticmoment=magneticmoment,
            weight=weight,
            nrejections=nrejections,
            initialeomid=self.initialeomid,
            rng=ctx.rng
        )
    elseif self.hybridscheme == HybridScheme.SaveSwitches
        params = HybridParamsWithDetection(
            charge=prob.p.charge,
            mass=prob.p.mass,
            electromagneticfield=prob.p.electromagneticfield,
            magneticmoment=magneticmoment,
            weight=weight,
            nrejections=nrejections,
            initialeomid=self.initialeomid,
            rng=ctx.rng
        )
    else
        error("Invalid hybridscheme: $(self.hybridscheme).")
    end
    return remake(
        prob;
        u0=u0,
        tspan=(t, self.tf),
        p=params
    )
end


"""
    init_probfunc(
        probfunc_specifier,
        p::TraceParticlesParameters,
        electromagneticfield
    )

Resolve `probfunc_specifier` into a `prob_func` that `EnsembleProblem` can call as
`(prob, ctx)`.

Anything that is not a recognised specification is returned unchanged, so a
user-supplied problem function needs no method here.
See [`Rerun`](@ref) and [`ICsFromFile`](@ref) examples of  specifications
that are recognised.
"""
function init_probfunc end

init_probfunc(pfs::Any,  _, _) = pfs


"""
    ICsFromFile(; ic_file::String)
Take the initial conditions from `ic_file`, an HDF5 file. Resolves into a
[`PredefinedICs`](@ref).

## Arguments
-`ic_file` is assumed to be a JLD2 file containing
the initial positions, initial velocities, initial times, charge, mass, magnetic
moment, statistical weight, and number of rejections (from sampling the condition)
for the particles to simulate.
"""
Base.@kwdef struct ICsFromFile
    ic_file::String
end

function init_probfunc(probfuncspec::ICsFromFile, params, electromagneticfield)
    u0, tspans, odeparams = create_diffeq_ic(
        probfuncspec.ic_file,
        params.tf,
        electromagneticfield;
        f=params.eom
    )
    return PredefinedICs(u0, tspans, odeparams)
end

Base.@kwdef struct Rerun{T<:AbstractVector}
    datafile::String
    idxs::T
end

function init_probfunc(pfs::Rerun, params, electromagneticfield)
    x0 = h5_getdataset(pfs.datafile, "x0")
    y0 = h5_getdataset(pfs.datafile, "y0")
    z0 = h5_getdataset(pfs.datafile, "z0")
    t0 = h5_getdataset(pfs.datafile, "t0")
    tf = h5_getdataset(pfs.datafile, "tf")
    tspans = [(t0[i], tf[i]) for i in pfs.idxs]
    if params.eom == guidingcentreapproximation!
        vparal0 = h5_getdataset(pfs.datafile, "vparal0")
        u0 = [[x0[i], y0[i], z0[i], vparal0[i]] for i in pfs.idxs]
        mu0 = h5_getdataset(pfs.datafile, "magneticmoment")[pfs.idxs]
    elseif params.eom == lorentzforce!
        vx0 = h5_getdataset(pfs.datafile, "vx0")
        vy0 = h5_getdataset(pfs.datafile, "vy0")
        vz0 = h5_getdataset(pfs.datafile, "vz0")
        mu0 = Vector{Bool}(undef, length(pfs.idxs))
        u0 = [
           [x0[i], y0[i], z0[i], vx0[i], vy0[i], vz0[i]]
           for i in pfs.idxs
        ]
    elseif params.eom == hybridgcafo!
        x0 = h5_getdataset(pfs.datafile, "x0")
        y0 = h5_getdataset(pfs.datafile, "y0")
        z0 = h5_getdataset(pfs.datafile, "z0")
        vx0 = h5_getdataset(pfs.datafile, "vx0")
        vy0 = h5_getdataset(pfs.datafile, "vy0")
        vz0 = h5_getdataset(pfs.datafile, "vz0")
        mu0 = h5_getdataset(pfs.datafile, "magneticmoment")[pfs.idxs]
        u0 = [
           [x0[i], y0[i], z0[i], vx0[i], vy0[i], vz0[i]]
           for i in pfs.idxs
        ]
    else
        @error "Cannot initialise `Rerun` problem function with uknown EoM:
            $(params.eom)"
    end
    odeparamstype = typeof(define_parameters(params, electromagneticfield))
    odeparams = Vector{odeparamstype}(undef, length(pfs.idxs))
    for i in eachindex(pfs.idxs)
        odeparams[i] = define_parameters(
            params, electromagneticfield; magneticmoment=mu0[i]
        )
    end
    return PredefinedICs(u0, tspans, odeparams)
end


"""
    init_probfun(sampler::AbstractParticleSampler, p, _)

Initialise a problem function which uses a `AbstractParticleSampler`.

Depending on the equations of motion, constructs either a 
[`SampleFullOrbit`](@ref), [`SampleGCA`](@ref) or [`SampleHybrid`](@ref),
instance.

See also [`MHDSampler`](@ref)
"""
function init_probfunc(sampler::AbstractParticleSampler, p, _)
    if p.eom == lorentzforce!
        return SampleFullOrbit(sampler, p.tf)
    elseif p.eom == guidingcentreapproximation!
        return SampleGCA(sampler, p.tf)
    elseif p.eom == hybridgcafo!
        scheme = p.save_switchinfo ?
            HybridScheme.SaveSwitches : HybridScheme.Default
        return SampleHybrid(sampler, p.tf, p.initialeomid, scheme)
    else
        throw(ArgumentError(
            "`SampleICsFromMHD` has no problem function for eom=$(p.eom)."
        ))
    end
end
#-------------------------------------------------------------------------------
# Functions for saving initial conditions to file
#-------------------------------------------------------------------------------


"""
    initialconditions_mhdsampling(
        npart::Int,
        electromagneticfield,
        EoM::Function,
        mass::Real,
        charge::Real,
        args...;
        initialeomid::Int=EoMID.GuidingCentreApproximation
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
    initialeomid::Int=EoMID.GuidingCentreApproximation
)
    u0s = [zeros(6) for _ in 1:npart]
    magneticmoments = Vector{Float64}(undef, npart)
    t0s = Vector{Float64}(undef, npart)
    weights = Vector{Float64}(undef, npart)
    nrejections = Vector{Int}(undef, npart)

    if (
        EoM == lorentzforce! ||
        (EoM == hybridgcafo! && initialeomid == EoMID.FullOrbit)
    )
        u0func! = getu0_fullorbit!
    elseif (
        EoM == guidingcentreapproximation! ||
        (EoM == hybridgcafo! && initialeomid == EoMID.GuidingCentreApproximation)
    )
        u0func! = getu0_guidingcentre!
    end
    Threads.@threads for i in 1:npart
        x, y, z, vx, vy, vz, t, weights[i], nrejections[i] = mhdsample(
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
    return float(typeof(x))(NaN)
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

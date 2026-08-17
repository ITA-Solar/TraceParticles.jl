"""
    TraceParticlesParameters(; kwargs...)

Parameters needed to assemble a [`TraceParticlesProblem`](@ref).

All fields are keyword arguments.

Optional callbacks are specified with a tuple of [`CallbackSpec`](@ref)s.
E.g. see [`OutOfBounds`](@ref), [`Relativistic`](@ref) or
[`HighMagneticGradient`](@ref).

Providing single-precision input consistently gives a single-precision run.

# Mandatory fields
- `eom`: equations of motion,
- `npart`: number of particles,
- `tf`: end time of the integration,
- `mass`, `charge`: particle mass and charge,
- `alg`: the `OrdinaryDiffEq` algorithm to integrate with,
- `maxiters`, `reltol`, `abstol`: solver controls,
- `expname`: experiment name; also the output subdirectory and file stem,
- `emfield_file`: JLD2 file holding a 6-element vector of callables in the
  order `bx, by, bz, ex, ey, ez`,
- `prob_func`: how each particle gets its initial state. Either a specification
  such as [`MHDSample`](@ref) or [`ICsFromFile`](@ref), or any function that
  `EnsembleProblem` can call as `(prob, ctx)`.

# Optional fields (not a complete list)
- `output_func`: The output the integrator is to save from each particle
   trajectory. If given it is directly passed as a keyword argument to the
   `SciMLBase.EnsembleProblem`. The default is to use
   `output_func_lightweight` or `output_func_lightweight_hybrid` depending on
   the equations of motion.
- `reduction`: The data reduction of each batch of particles. The default is
   `SaveBatchAsHDF5`, which requires the batch to be transformable into a
   dataframe via `DataFrame(batch)`.
- `callback_specs`: tuple of callback specifications, empty by default,
- `batchsize`: number of particles per reduction batch, defaults to `npart`,
- `datadir`: where the experiment directory is created,
- `seed`: random seed for the ensemble.

# Example
```julia
params = TraceParticlesParameters(
    eom = guidingcentreapproximation!,
    npart = 1000,
    tf = 1.0,
    mass = TraceParticles.m_e,
    charge = -TraceParticles.e,
    alg = Tsit5(),
    maxiters = 10_000_000,
    reltol = 3e-8,
    abstol = 1e0,
    expname = "myrun",
    emfield_file = "BE.jld2",
    prob_func = SampleICsFromMHD(
        tg_file = "tg.jld2",
        target_distr_file = "r.jld2",
        xbounds = (1.58e7, 1.63e7),
        ybounds = (1e6, 1e6),
        zbounds = (-6.6e6, -6.1e6),
    ),
    callback_specs = (
        OutOfBounds(x = (1.58e7, 1.63e7), z = (-6.6e6, -6.1e6)),
        Relativistic(fraction = 0.02),
    ),
)
```

See also [`TraceParticlesProblem`](@ref), [`validate`](@ref).
"""
Base.@kwdef struct TraceParticlesParameters
    eom::Function
    npart::Int
    tf::Real
    mass::Real
    charge::Real
    alg::SciMLBase.AbstractODEAlgorithm
    maxiters::Int
    reltol::Real
    abstol::Real
    expname::String
    datadir::String = "data"
    emfield_file::String
    prob_func::Any
    # Optional parameters
    output_func::Any = nothing
    reduction::Any = nothing
    callback_specs::Tuple = ()
    batchsize::Int = npart
    safetycopy::Bool = false
    save_max_observables::Bool = false
    verbose::Bool = true
    loglevel::Logging.LogLevel = Logging.Info
    static_itp_wrapper::Bool = false
    xz_itp_wrapper::Bool = false
    seed::Int = 1
    compute_initial_and_final_energy::Bool = true
    #= Parameters for the hybrid GCA/full-orbit switch. Unlike the callback
    specifications these are not optional add-ons: they are the switch itself,
    and apply whenever `eom == hybridgcafo!`. =#
    initialeomid::EoMID.T = EoMID.FullOrbit
    gradient_tolerance::Real = 1e-4
    curvature_tolerance::Real = 1e-4
    eparal_ratio::Real = Inf
    switchback_tolerance::Real = 0.9
    switch_save_positions::Tuple{Bool,Bool} = (true, true)
    save_switchinfo::Bool = false

    #= This inner constructor intercepts the positional call that
    `Base.@kwdef` generates, so that every instance is validated without having
    to spell out all the fields. =#
    function TraceParticlesParameters(args...)
        p = new(args...)
        validate(p)
        return p
    end
end

"""
    validate(p::TraceParticlesParameters)
Check the parameter combinations that the type system does not already rule
out. Throws an `ArgumentError` describing the problem, and returns `nothing`
when `p` is consistent.

Most requirements are encoded in the types instead: e.g. an
[`OutOfBounds`](@ref) cannot be constructed without bounds. What is left are
the checks that span several fields.

Called automatically on construction.
"""
function validate(p::TraceParticlesParameters)
    if p.npart < 1
        throw(ArgumentError("`npart` must be positive, got $(p.npart)."))
    end
    if p.batchsize < 1
        throw(ArgumentError("`batchsize` must be positive, got $(p.batchsize)."))
    end
    if !isfinite(p.tf)
        throw(ArgumentError("`tf` must be finite, got $(p.tf)."))
    end

    notaspec = filter(s -> !(s isa CallbackSpec), collect(p.callback_specs))
    if !isempty(notaspec)
        throw(ArgumentError(
            "`callback_specs` must hold `CallbackSpec`s, got "*
            "$(join(typeof.(notaspec), ", "))."
        ))
    end

    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", p::TraceParticlesParameters)
    names = fieldnames(TraceParticlesParameters)
    pad = maximum(length ∘ string, names)
    println(io, "TraceParticlesParameters:")
    for name in names
        println(io, "  ", rpad(name, pad), " = ", repr(getfield(p, name)))
    end
end

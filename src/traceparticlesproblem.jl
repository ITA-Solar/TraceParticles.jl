"""
    TraceParticlesProblem(p::TraceParticlesParameters; logger)

An assembled test-particle simulation, ready to be handed to `solve`.

The constructor is the assembly step: it loads the electromagnetic field from
`p.emfield_file`, wraps every field component so that it takes `(x, y, z, t)`,
initialises the `prob_func`, picks the `output_func` and `reduction`, and
creates callbacks implied by `p`. The result is a SciML `EnsembleProblem` plus
the pieces `solve` needs alongside it.

# Fields
- `ensemble_prob`: the `EnsembleProblem` to integrate.
- `callbackset`: the `CallbackSet` built from the termination and equation-of-
  motion switching parameters in `p`.
- `trajectories`: number of particles, i.e. `p.npart`.
- `batch_size`: trajectories per batch, capped at `trajectories`.
- `logger`: the logger the simulation writes to. Defaults to one that tees to
  `stderr` and to a log file in the experiment directory.

See also [`TraceParticlesParameters`](@ref).
"""
@kwdef struct TraceParticlesProblem{
    ProbType<:EnsembleProblem,
    CallbacksetType<:CallbackSet,
    I<:Int,
    T,
    L<:Logging.AbstractLogger
}
    ensemble_prob::ProbType
    callbackset::CallbacksetType
    trajectories::I
    batch_size::T
    logger::L
end

function TraceParticlesProblem(
    p::TraceParticlesParameters;
    logger=current_logger()
)
    create_expdir(p.datadir, p.expname)
    logger = determine_logger(logger, p)

    with_logger(logger) do
        electromagneticfield = setup_emfield(p.electromagneticfield)
        prob_func = invokelatest(setup_probfunc, p, electromagneticfield)
        return invokelatest(build_problem, p, logger, electromagneticfield, prob_func)
    end
end

function build_problem(
    p::TraceParticlesParameters,
    logger::AbstractLogger,
    electromagneticfield::Any,
    prob_func::Any
)
    cbs = create_callbacks(p)
    # Determine output functions;
    # the data footprint and format of output data for each particle
    output_func = determine_output_func(p)
    # Determine reduction
    reduction, effective_batchsize = determine_reduction(p)
    # Define the type of problem to solve
    odeprob = define_odeproblem(p, electromagneticfield)

    ## Finally define the ensemble problem
    ensemble_prob = EnsembleProblem(
        odeprob;
        prob_func = prob_func,
        safetycopy = p.safetycopy,
        output_func = output_func,
        reduction = reduction
    )

    # Return the initialised problem
    return TraceParticlesProblem(
        ensemble_prob=ensemble_prob,
        callbackset=CallbackSet(cbs...),
        trajectories=p.npart,
        batch_size=effective_batchsize,
        logger=logger
    )
end

"""
    setup_emfield(emfield::Any)
Expect `emfield` to be the electromagneticfield-callable itself.
"""
setup_emfield(emfield::Any) = emfield
"""
    setup_emfield(prob_func_setup::AbstractString)
Include the setup file and expect a returned electromagnetic field.
"""
function setup_emfield(emfield_setup::AbstractString)
    @info "Setting up EM-field..."
    time_taken = @timed emfield = include_setup(emfield_setup)
    logtimed(time_taken, "EM-field set up")
    return emfield
end

"""
    setup_probfunc(params, electromagneticfield)
Include the setup file and expect a returned problem function to be initialised
with `init_probfunc`.
"""
function setup_probfunc(params, electromagneticfield)
    if params.prob_func isa String
        @info "Setting up problem function..."
        time_taken = @timed prob_func = init_probfunc(
            include_setup(params.prob_func),
            params,
            electromagneticfield
        )
        logtimed(time_taken, "prob_func set up")
        return prob_func
    else
        return init_probfunc(params.prob_func, params, electromagneticfield)
    end
    return prob_func
end

"""
    create_callbacks(p::TraceParticlesParameters)
Build the callbacks implied by `p`: one per entry of `p.callback_specs`, plus
the GCA/full-orbit switch callback when the equations of motion are hybrid.
"""
function create_callbacks(p::TraceParticlesParameters)
    # Each callback spec knows how to build its own callback.
    cbs = Any[create_callback(spec, p) for spec in p.callback_specs]
    # The hybrid switch is not an optional add-on but part of the equations of
    # motion, so it is not a callback spec.
    if p.eom == hybridgcafo!
        affect = p.save_switchinfo ?
            hybridswitchaffect_withdetection! :
            hybridswitchaffect!
        switch_cb = DiscreteCallback(
            HybridSwitchCondition(
                p.gradient_tolerance,
                p.curvature_tolerance,
                p.eparal_ratio,
                p.switchback_tolerance
            ),
            affect,
            save_positions=p.switch_save_positions,
        )
        push!(cbs, switch_cb)
    end
    return cbs
end

function determine_output_func(p::TraceParticlesParameters)
    isnothing(p.output_func) || return p.output_func
    if p.eom == hybridgcafo!
        return p.save_switchinfo ?
            (sol, _) -> (sol, false) :
            output_func_lightweight_hybrid
    elseif p.eom == guidingcentreapproximation!
        return p.save_max_observables ?
            output_func_max_lightweight : output_func_lightweight_gca
    elseif p.eom == lorentzforce!
        return output_func_lightweight_fo
    else
        throw(ArgumentError(
        """Could not infer `output_func` from the equations of motion.
           Please provide one as a parameter."""
        ))
    end
end

function determine_reduction(p::TraceParticlesParameters)
    # Set batchsize to npart if the npart < batchsize
    effective_batchsize = p.npart > p.batchsize ? p.batchsize : p.npart
    isnothing(p.reduction) || return p.reduction, effective_batchsize
    nbatches = ceil(Int, p.npart / effective_batchsize)
    reduction = (p.eom == hybridgcafo! && p.save_switchinfo) ?
        (u, data, _) -> (append!(u, data), false) :
        SaveBatchAsHDF5(
            expdir=p.datadir,
            expname=p.expname,
            batchsize=effective_batchsize,
            nbatches=nbatches,
            verbose=p.verbose
        )
    return reduction, effective_batchsize
end

function define_parameters(
    p::TraceParticlesParameters, electromagneticfield;
    magneticmoment=zero(p.mass)
)
    if p.eom == guidingcentreapproximation!
        return GCAParams(
            charge=p.charge,
            mass=p.mass,
            electromagneticfield=electromagneticfield,
            magneticmoment = magneticmoment
        )
    elseif p.eom == lorentzforce!
        return FullOrbitParams(
           charge=p.charge,
           mass=p.mass,
           electromagneticfield=electromagneticfield,
        )
    elseif p.eom == hybridgcafo! && p.save_switchinfo
        return HybridParamsWithDetection(
            charge=p.charge,
            mass=p.mass,
            electromagneticfield=electromagneticfield,
            magneticmoment = magneticmoment,
            initialeomid=p.initialeomid
        )
    elseif p.eom == hybridgcafo!
        return HybridParams(
            charge=p.charge,
            mass=p.mass,
            electromagneticfield=electromagneticfield,
            magneticmoment = magneticmoment,
            initialeomid=p.initialeomid
        )
    else
        throw(ArgumentError("Unknown prob_func for f=$(p.eom)."))
    end
end

function define_odeproblem(p::TraceParticlesParameters, electromagneticfield)
    odeparams = define_parameters(p, electromagneticfield)
    ## Define the type of ODEProblem
    ## I think the important thing here is the type of the problem, not
    ## the values because they will be overwritten by the problem function
    ## passed to the EnsembleProblem.
    return ODEProblem(
        p.eom,
        # Only the type of the time span matters here; `prob_func` overwrites
        # the values with each particle's own initial time.
        zeros(typeof(p.mass), 1),
        (zero(p.tf), p.tf),
        odeparams
    )
end

function determine_logger(logger, p::TraceParticlesParameters)
    if p.log2file
        return logger2fileandstderr(
            p.loglevel,
            joinpath(p.datadir, p.expname, p.expname*".log.txt")
        )
    else
        return logger
    end
end

function get_eom(prob::TraceParticlesProblem)
    return prob.ensemble_prob.prob.f.f
end

function get_electromagneticfield(prob::TraceParticlesProblem)
    return prob.ensemble_prob.prob.p.electromagneticfield
end

function get_prob_func(prob::TraceParticlesProblem)
    return prob.ensemble_prob.prob_func
end

function get_output_func(prob::TraceParticlesProblem)
    return prob.ensemble_prob.output_func
end

function get_reduction(prob::TraceParticlesProblem)
    return prob.ensemble_prob.reduction
end

function get_filename(prob::TraceParticlesProblem)
    return get_filename(prob.ensemble_prob.reduction)
end

"""
    experimentname(prob::TraceParticlesProblem)
The name of the experiment, taken from the file its reduction writes. A
reduction that writes no file of its own, has no name to give.
"""
function experimentname(prob::TraceParticlesProblem)
    reduction = get_reduction(prob)
    applicable(get_filename, reduction) || return "(unnamed)"
    return first(splitext(basename(get_filename(reduction))))
end

function get_reductiontype(prob::TraceParticlesProblem)
    return typeof(prob.ensemble_prob.reduction)
end

function get_nofbatches(prob::TraceParticlesProblem)
    return ceil(Int, prob.trajectories/prob.batch_size)
end

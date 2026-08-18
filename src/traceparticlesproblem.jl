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

        # Construct a wrapper that wraps the interpolator in a function with a
        # call signature (x, y, z, t)
        itp_wrapper = construct_itpwrapper(p)

        #=
        Now load the electromagnetic field interpolation object.
        This step assumes that the file `emfield_file` points to a JLD2-file containing
        a vector of callable objects that each return a component of the electric
        or magnetic vector field. Specifically it assumes the vector to contain the
        components in the following order
            vec[1] -> bx
            vec[2] -> by
            vec[3] -> bz
            vec[4] -> ex
            vec[5] -> ey
            vec[6] -> ez
        =#
        time_taken = @timed fields_itp = ElectromagneticFieldInterpolator(
            [itp_wrapper(component) for component in load_object(p.emfield_file)]
        )
        logtimed(time_taken, "EM-field loaded")

        # Resolve `p.prob_func` into the function that gives each particle its
        # initial state.
        prob_func = init_probfunc(p.prob_func, p, itp_wrapper, fields_itp)

        cbs = create_callbacks(p)
        # Set batchsize to npart if the npart < batchsize
        effective_batchsize = p.npart > p.batchsize ? p.batchsize : p.npart
        # Determine output functions;
        # the data footprint and format of output data for each particle
        output_func = determine_output_func(p)
        # Determine reduction
        reduction = determine_reduction(p, effective_batchsize)
        # Define the type of problem to solve
        ## I think the important thing here is the type of the problem, not
        ## the values because they will be overwritten by the problem function
        ## passed to the EnsembleProblem.
        odeprob = define_odeproblem(p, fields_itp)

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
end

"""
    get_itpwrapper(p::TraceParticlesParams)
Construct a function to wrap the interpolation objects into a new
object with appropriate call signatures.

E.g:
An interpolation object with original call signature `(x, z)` should be
wrapped in the `StaticXZInterpolation` struct such that it can be called with
the signature `itp(x, y, z, t) -> itp(x, z)`.
"""
function construct_itpwrapper(p::TraceParticlesParameters)
    if p.static_itp_wrapper
        if p.xz_itp_wrapper
            return (obj) -> StaticXZInterpolation(obj)
        else
            return (obj) -> StaticInterpolation(obj)
        end
    elseif p.xz_itp_wrapper
        return (obj) -> XZInterpolation(obj)
    else
       return (obj) -> obj
    end
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
            output_func = (sol, ctx) -> (sol, false) :
            output_func_lightweight_hybrid
    elseif p.eom == guidingcentreapproximation!
        return p.save_max_observables ?
            output_func_max_lightweight : output_func_lightweight
    elseif p.eom == lorentzforce!
        return output_func_lightweight_fo
    else
        throw(ArgumentError(
        """Could not infer `output_func` from the equations of motion.
           Please provide one as a parameter."""
        ))
    end
end

function determine_reduction(p::TraceParticlesParameters, effective_batchsize)
    isnothing(p.reduction) || return p.reduction
    nbatches = ceil(Int, p.npart / effective_batchsize)
    return (p.eom == hybridgcafo! && p.save_switchinfo) ?
        (u, data, I) -> (append!(u, data), false) :
        SaveBatchAsHDF5(
            expdir=p.datadir,
            expname=p.expname,
            batchsize=effective_batchsize,
            nbatches=nbatches,
            verbose=p.verbose
        )
end

function define_parameters(
    p::TraceParticlesParameters, fields_itp;
    magneticmoment=zero(p.mass)
)
    if p.eom == guidingcentreapproximation!
        return GCAParams(
            charge=p.charge,
            mass=p.mass,
            electromagneticfield=fields_itp,
            magneticmoment = magneticmoment
        )
    elseif p.eom == lorentzforce!
        return FullOrbitParams(
           charge=p.charge,
           mass=p.mass,
           electromagneticfield=fields_itp,
        )
    elseif p.eom == hybridgcafo! && p.save_switchinfo
        return HybridParamsWithDetection(
            charge=p.charge,
            mass=p.mass,
            electromagneticfield=fields_itp,
            magneticmoment = magneticmoment,
            initialeomid=p.initialeomid
        )
    elseif p.eom == hybridgcafo!
        return HybridParams(
            charge=p.charge,
            mass=p.mass,
            electromagneticfield=fields_itp,
            magneticmoment = magneticmoment,
            initialeomid=p.initialeomid
        )
    else
        throw(ArgumentError("Unknown prob_func for f=$(p.eom)."))
    end
end

function define_odeproblem(p::TraceParticlesParameters, fields_itp)
    odeparams = define_parameters(p, fields_itp)
    ## Define the type of ODEProblem
    odeprob = ODEProblem(
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

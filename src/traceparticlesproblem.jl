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
    batch_size::T = nothing
    logger::L = global_logger()
end

function TraceParticlesProblem(
    p::TraceParticlesParameters;
    logger=nothing
)
    create_expdir(p.datadir, p.expname)
    logger = if isnothing(logger)
        logger2fileandstderr(
            p.loglevel,
            joinpath(p.datadir, p.expname, p.expname*".log.txt")
        )
    else
        logger
    end
    with_logger(logger) do

        # Construct a wrapper that wraps the interpolator in a function with a
        # call signature (x, y, z, t)
        itp_wrapper = construct_itpwrapper(p)

        # Create problem function for computing initial conditions on-the-fly
        if p.ic_onthefly
            prob_func = construct_ic_onthefly_prob_func(p, itp_wrapper)
        end

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
        time_taken = @timed fields_itp = itp_wrapper(
            ElectromagneticFieldInterpolator(
                load_object(p.emfield_file)
            )
        )
        @info """EM-field loaded:
        Elapsed time: $(time_taken.time) seconds
        Allocations: $(@sprintf "%i" time_taken.bytes/1e6) MB"""

        # Create problem function for using pre-defined initial conditions
        if !p.ic_onthefly
            prob_func = construct_ic_predefined_prob_func()
        end

        cbs = construct_callbacks(p)
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
    return p.static_itp_wrapper ?
        begin
            p.xz_itp_wrapper ?
                (obj) -> StaticXZInterpolation(obj) :
                (obj) -> StaticInterpolation(obj)
        end :
       (obj) -> obj
end

function construct_ic_onthefly_prob_func(
    p::TraceParticlesParameters,
    itp_wrapper
)
    (; eom, tf) = p
    # Load temperature, density and assign proposal and target distributions
    tg_itp = itp_wrapper(load_object(p.tg_file))
    # Target distribution for positional sampling
    target_distr = itp_wrapper(load_object(p.target_distr_file))
    if isempty(p.proposal_distr_file)
        prop_distr = target_distr
    else
        prop_distr = itp_wrapper(load_object(p.proposal_distr_file))
    end
    maxval = maximum(prop_distr) # Needed by the rejection algorithm.
    # Create the problem function
    # determines the initial conditions of the particles,
    # including statistical weight and initial and final simulation times.
    xbounds = (p.xstart_ic, p.xend_ic)
    ybounds = (p.ystart_ic, p.yend_ic)
    zbounds = (p.zstart_ic, p.zend_ic)
    t0bounds = (p.t0start, p.t0end)
    (eomid, hybridscheme) = 
        if eom == lorentzforce!
            (EoMID.FullOrbit, HybridScheme.No)
        elseif eom == guidingcentreapproximation!
            (EoMID.GuidingCentreApproximation, HybridScheme.No)
        elseif eom == hybridgcafo! && p.save_switchinfo
            (p.initialeomid, HybridScheme.YesSaveSwitches)
        elseif eom == hybridgcafo!
            (p.initialeomid, HybridScheme.Yes)
        end
    prob_func = MHDSampling(
        prop_distr,
        target_distr,
        tg_itp,
        xbounds,
        ybounds,
        zbounds,
        t0bounds,
        tf,
        maxval,
        eomid,
        hybridscheme
    )
end

function construct_callbacks(p::TraceParticlesParameters)
    cbs = []
    # If terminate particles out of bounds
    if p.kill_oob
        (; xkill, ykill, zkill) = p
        (; xlowerbound, xupperbound) = p
        (; ylowerbound, yupperbound) = p
        (; zlowerbound, zupperbound) = p

        kill_at_axis = (xkill, ykill, zkill)
        axisbounds = (
            (xlowerbound, xupperbound),
            (ylowerbound, yupperbound),
            (zlowerbound, zupperbound),
        )
        axes_idxs = [i for i in 1:3 if kill_at_axis[i]]
        axes_bounds = Tuple(axisbounds[i] for i in axes_idxs)

        out_cb = DiscreteCallback(
            OutOfBoundsCondition(axes_bounds, axes_idxs),
            outofboundsaffect!,
            save_positions = p.oob_save_positions
        )
        push!(cbs, out_cb)
    end
    # If kill relativistic particles
    (; mass, eom) = p
    if p.kill_relativistic
        relativistic_condition = 
            if eom == hybridgcafo!
                RelativisticConditionHybrid
            elseif eom == guidingcentreapproximation!
                RelativisticConditionGCA
            elseif eom == lorentzforce!
                RelativisticCondition
            else
                throw(ArgumentError("Unknown prob_func for f=$eom."))
            end
        rel_cb = DiscreteCallback(
            relativistic_condition(; mass=mass, fraction=p.relativistic_fraction),
            relativisticaffect!,
            save_positions = p.rel_save_positions
        )
        push!(cbs, rel_cb)
    end
    # If kill when high gradient in magnetic field
    if p.kill_high_gradb
        gca_cb = DiscreteCallback(
            MagneticGradientCondition(p.gradient_tolerance),
            magneticgradientaffect!,
            save_positions = p.gradb_save_positions
        )
        push!(cbs, gca_cb)
    end
    # If hybrid scheme is used, we need to create the switch callback
    if eom == hybridgcafo!
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
            save_positions=(true,true),
        )
        push!(cbs, switch_cb)
    end
    return cbs
end

function construct_ic_predefined_prob_func(
    p::TraceParticlesParameters,
    fields_itp
)
    u0, tspan, params = create_diffeq_ic(
        p.ic_file,
        p.tf,
        fields_itp;
        f=p.eom
    )
    return PredefinedICs(u0, tspan, params)
end

function determine_output_func(p::TraceParticlesParameters)
    if p.eom == hybridgcafo!
        return p.save_switchinfo ?
            output_func = (sol, ctx) -> (sol, false) :
            output_func_lightweight_hybrid
    elseif p.eom == guidingcentreapproximation!
        return p.save_max_observables ?
            output_func_max_lightweight : output_func_lightweight
    elseif p.eom == lorentzforce!
        return output_func_lightweight_fo
    end
end

function determine_reduction(p::TraceParticlesParameters, effective_batchsize)
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

function define_parameters(p::TraceParticlesParameters, fields_itp)
    if p.eom == guidingcentreapproximation!
        return GCAParams(
            charge=p.charge,
            mass=p.mass,
            electromagneticfield=fields_itp,
            magneticmoment = p.precision(0.0)
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
            magneticmoment = p.precision(0.0),
            initialeomid=p.initialeomid
        )
    elseif p.eom == hybridgcafo!
        return HybridParams(
            charge=p.charge,
            mass=p.mass,
            electromagneticfield=fields_itp,
            magneticmoment = p.precision(0.0),
            initialeomid=p.initialeomid
        )
    else
        throw(ArgumentError("Unknown prob_func for f=$eom."))
    end
end

function define_odeproblem(p::TraceParticlesParameters, fields_itp)
    odeparams = define_parameters(p, fields_itp)
    ## Define the type of ODEProblem
    odeprob = ODEProblem(
        p.eom,
        # I think the important thing is the type of u0, not the length
        zeros(p.precision, 1),
        (p.t0start, p.tf),
        odeparams
    )
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

function get_reductiontype(prob::TraceParticlesProblem)
    return typeof(prob.ensemble_prob.reduction)
end

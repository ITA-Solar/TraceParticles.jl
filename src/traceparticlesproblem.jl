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
        #=
        First, construct a function to wrap the interpolation objects into a new
        object with appropriate call signatures.

        E.g:
        An interpolation object with original call signature `(x, z)` should be
        wrapped in the `StaticXZInterpolation` struct such that it can be called with
        the signature `itp(x, y, z, t) -> itp(x, z)`.
        =#
        itp_wrapper = p.static_itp_wrapper ?
            begin
                p.xz_itp_wrapper ?
                    (obj) -> StaticXZInterpolation(obj) :
                    (obj) -> StaticInterpolation(obj)
            end :
           (obj) -> obj

        #=
        Next, if the particle's initial conditions are sampled on-the-fly we need to
        create the appropriate problem function.
        =#
        (; ic_onthefly, eom, tf) = p
        if ic_onthefly
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
            rng = Xoshiro(p.seed) # Random number generator
            # Create the problem function
            # determines the initial conditions of the particles,
            # including statistical weight and initial and final simulation times.
            xbounds = (p.xstart_ic, p.xend_ic)
            ybounds = (p.ystart_ic, p.yend_ic)
            zbounds = (p.zstart_ic, p.zend_ic)
            t0bounds = (p.t0start, p.t0end)
            if eom == hybridgcafo!
                prob_func = MHDSamplingHybrid(
                    rng,
                    prop_distr,
                    target_distr,
                    tg_itp,
                    xbounds,
                    ybounds,
                    zbounds,
                    t0bounds,
                    tf,
                    maxval,
                    p.initialeomid
                )
            elseif eom == guidingcentreapproximation!
                prob_func = MHDSamplingGCA(
                    rng,
                    prop_distr,
                    target_distr,
                    tg_itp,
                    xbounds,
                    ybounds,
                    zbounds,
                    t0bounds,
                    tf,
                    maxval,
                )
            else
                throw(ArgumentError("Unknown prob_func for f=$eom."))
            end
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
        Allocations: $(time_taken.bytes) bytes"""

        #=
        Create the appropriate problem function for the case where the initial conditions
        are predefined.
        =#
        if !ic_onthefly
            u0, tspan, params = create_diffeq_ic(
                p.ic_file,
                tf,
                p.fields_itp;
                f=eom
            )
            prob_func = PredefinedICs(u0, tspan, params)
        end

        # Define Callbacks
        cbs = []
        # If terminate particles out of bounds
        if p.kill_oob
            (; xkill, ykill, zkill) = p
            (; xlowerbound, xupperbound) = p
            (; ylowerbound, yupperbound) = p
            (; zlowerbound, zupperbound) = p
            if xkill && !ykill && !zkill
                axes_bounds = ((xlowerbound, xupperbound))
                axes_idxs = [1]
            elseif !xkill && ykill && !zkill
                axes_bounds = ((ylowerbound, yupperbound))
                axes_idxs = [2]
            elseif !xkill && !ykill && zkill
                axes_bounds = ((zlowerbound, zupperbound))
                axes_idxs = [3]
            elseif xkill && ykill && !zkill
                axes_bounds = (
                    (xlowerbound, xupperbound),
                    (ylowerbound, yupperbound)
                )
                axes_idxs = 1:2
            elseif xkill && !ykill && zkill
                axes_bounds = (
                    (xlowerbound, xupperbound),
                    (zlowerbound, zupperbound)
                )
                axes_idxs = 1:2:3
            elseif !xkill && ykill && zkill
                axes_bounds = (
                    (ylowerbound, yupperbound),
                    (zlowerbound, zupperbound)
                )
                axes_idxs = 2:3
            elseif xkill && ykill && zkill
                axes_bounds = (
                    (xlowerbound, xupperbound),
                    (ylowerbound, yupperbound),
                    (zlowerbound, zupperbound)
                )
                axes_idxs = 1:3
            end
            out_cb = DiscreteCallback(
                OutOfBoundsCondition(axes_bounds, axes_idxs),
                outofboundsaffect!,
                save_positions = p.oob_save_positions
            )
            push!(cbs, out_cb)
        end
        # If kill relativistic particles
        (; mass) = p
        if p.kill_relativistic
            rel_cb = DiscreteCallback(
            RelativisticConditionGCA(; mass=mass, fraction=p.relativistic_fraction),
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
            switch_cb = DiscreteCallback(
                HybridSwitchCondition(
                    p.gradient_tolerance,
                    p.curvature_tolerance,
                    p.eparal_ratio,
                    p.switchback_tolerance
                ),
                hybridswitchaffect!,
                save_positions=(true,true),
            )
            push!(cbs, switch_cb)
        end
        # Set batchsize to npart if the npart < batchsize
        (; npart, batchsize) = p
        BATCHSIZE = npart > batchsize ? batchsize : npart
        nbatches = ceil(Int, npart / BATCHSIZE)

        # Determine output functions
        # The data footprint and format of output data for each particle
        if eom == hybridgcafo!
            output_func = output_func_lightweight_hybrid
        elseif eom == guidingcentreapproximation!
            output_func =
                p.save_max_observables ?
                    output_func_max_lightweight : output_func_lightweight
        end

        # Define the type of problem to solve
        ## I think the important thing here is the type of the problem, not
        ## the values because they will be overwritten by the problem function
        ## passed to the EnsembleProblem.

        ## Define the type of problem parameters
        (; precision, charge) = p

        if eom == hybridgcafo!
            ODEParamStruct = HybridParams
        elseif eom == guidingcentreapproximation!
            ODEParamStruct = GCAParams
        else
           throw(ArgumentError("Unknown prob_func for f=$eom."))
        end
        odeparams = ODEParamStruct(
            charge=charge,
            mass=mass,
            electromagneticfield=fields_itp,
            magneticmoment = precision(0.0)
        )
        ## Define the type of ODEProblem
        odeprob = ODEProblem(
            eom,
            # I think the important thing is the type of u0, not the length
            zeros(p.precision, 1),
            (p.t0start, tf),
            odeparams
        )

        ## Finally define the ensemble problem
        ensemble_prob = EnsembleProblem(
            odeprob;
            prob_func = prob_func,
            output_func = output_func,
            reduction = SaveBatchAsHDF5(
                expdir=p.datadir,
                expname=p.expname,
                batchsize=BATCHSIZE,
                nbatches=nbatches,
                verbose=p.verbose
            ),
            safetycopy = p.safetycopy,
        )

        # Return the initialised problem
        return TraceParticlesProblem(
            ensemble_prob=ensemble_prob,
            callbackset=CallbackSet(cbs...),
            trajectories=npart,
            batch_size=BATCHSIZE,
            logger=logger
        )
    end
end

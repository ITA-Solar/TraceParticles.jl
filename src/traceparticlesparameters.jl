@with_kw struct TraceParticlesParameters{
        D<:AbstractFloat,
        I<:Int,
        F<:Function,
        S<:String,
        B<:Bool,
        TB<:Tuple{B, B},
        DT<:DataType,
        AlgType<:SciMLBase.AbstractODEAlgorithm,
        L<:Logging.LogLevel
    } @deftype D
    eom::F
    npart::I
    t0start
    t0end
    tf
    mass
    charge
    alg::AlgType
    maxiters::I
    reltol
    abstol
    expname::S
    datadir::S = "data"
    batchsize::I
    emfield_file::S
    ic_onthefly::B
    precision::DT
    # Optional parameters
    safetycopy::B = false
    save_max_observables::B = false
    verbose::B = true
    loglevel::L = Logging.Info
    static_itp_wrapper::B = false
    xz_itp_wrapper::B = false
    initialeomid::I = 2
    tg_file::S = ""
    target_distr_file::S = ""
    proposal_distr_file::S = ""
    xstart_ic = precision(NaN)
    xend_ic = precision(NaN)
    ystart_ic = precision(NaN)
    yend_ic = precision(NaN)
    zstart_ic = precision(NaN)
    zend_ic = precision(NaN)
    seed::I = 1
    ic_file::S = ""
    ## Parameters for callbacks
    ### Hybrid switch
    gradient_tolerance = precision(1e-4)
    curvature_tolerance = precision(1e-4)
    eparal_ratio = precision(Inf)
    switchback_tolerance = precision(0.9)
    switch_save_positions::TB = (true, true)
    ### For terminating the particle when it is out-of-bounds
    kill_oob::B = false
    xkill::B = false
    ykill::B = false
    zkill::B = false
    xlowerbound = precision(NaN)
    xupperbound = precision(NaN)
    ylowerbound = precision(NaN)
    yupperbound = precision(NaN)
    zlowerbound = precision(NaN)
    zupperbound = precision(NaN)
    oob_save_positions::TB = (false, true)
    ### For terminating the particle when it becomes relativistic
    kill_relativistic::B = false
    relativistic_fraction = precision(0.12)
    rel_save_positions::TB = (true, false)
    ### For terminating particles in areas of high gradients in the
    ### magnetic field
    kill_high_gradb::B = false
    gradb_save_positions::TB = (true, false)

    # Requirements
    ## Certain features requires a set of otherwise optional parameters.
    ## Check for missing parameters in these cases.
    @assert begin
        if kill_oob && (
            (xkill && (isnan(xlowerbound) || isnan(xupperbound))) ||
            (ykill && (isnan(ylowerbound) || isnan(yupperbound))) ||
            (zkill && (isnan(zlowerbound) || isnan(zupperbound))) ||
            (!xkill && !ykill && !zkill)
        )
            false
        else
            true
        end
    end """Terminating out of bounds particles requires defining:
            xkill::Bool, xlowerbound, xupperbound, and/or
            ykill::Bool, ylowerbound, yupperbound, and/or
            zkill::Bool, zlowerbound, zupperbound
    """

    @assert begin
        if ic_onthefly && (
            isempty(tg_file) ||
            isempty(target_distr_file) ||
            isnan(xstart_ic) ||
            isnan(xend_ic) ||
            isnan(ystart_ic) ||
            isnan(yend_ic) ||
            isnan(zstart_ic) ||
            isnan(zend_ic)
        )
            false
        else
            true
        end
    end """Calculating initial conditions on-the-fly requires defining:
            tg_file
            target_distr_file
            xstart_ic
            xend_ic
            ystart_ic
            yend_ic
            zstart_ic
            zend_ic
        For importance sampling, also define:
            proposal_distr_file
    """

    @assert begin
        if !ic_onthefly && isempty(ic_file)
            false
        else
            true
        end
    end """Using pre-calculated initial conditions requires defining:
        ic_file"""
end

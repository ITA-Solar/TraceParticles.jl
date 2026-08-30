#-------------------------------------------------------------------------------
#           test_traceparticlesproblem.jl
#-------------------------------------------------------------------------------
# Unit tests for the assembly step: the helpers that pick the interpolator
# wrapper, callbacks, problem function, output function, reduction and ODE
# problem from a `TraceParticlesParameters`, plus an end-to-end assembly of a
# `TraceParticlesProblem` from synthetic on-disk inputs.
#
# The synthetic fields are homogeneous and written to JLD2 files, so these
# tests need no Bifrost snapshot.
#-------------------------------------------------------------------------------

using HDF5
using Interpolations
using JLD2
using Logging
using OrdinaryDiffEq
using Random
using SciMLBase
using Test

import TraceParticles as TP
using TraceParticles:
    ElectromagneticFieldInterpolator,
    EoMID,
    FullOrbitParams,
    GCAParams,
    HighMagneticGradient,
    HybridParams,
    HybridParamsWithDetection,
    HybridScheme,
    ICsFromFile,
    SampleICsFromMHD,
    MHDSampler,
    OutOfBounds,
    PredefinedICs,
    Relativistic,
    SampleFullOrbit,
    SampleGCA,
    SampleHybrid,
    SaveBatchAsHDF5,
    StaticInterpolation,
    StaticXZInterpolation,
    TraceParticlesParameters,
    TraceParticlesProblem,
    XZInterpolation,
    construct_itpwrapper,
    create_callbacks,
    define_odeproblem,
    define_parameters,
    determine_output_func,
    determine_reduction,
    guidingcentreapproximation!,
    hybridgcafo!,
    init_probfunc,
    lorentzforce!,
    output_func_lightweight,
    output_func_lightweight_fo,
    output_func_lightweight_hybrid,
    output_func_max_lightweight

if !isdefined(Main, :verbose)
    verbose = 4
end

"""
A serialisable homogeneous field. An anonymous function cannot be
round-tripped through JLD2, so the emfield and distribution files hold these
instead. `maximum` is what the rejection algorithm needs from a proposal
distribution.
"""
struct ConstantField{T<:Real}
    value::T
end
(self::ConstantField)(x, y, z, t) = self.value
(self::ConstantField)(x, y, z) = self.value
Base.maximum(self::ConstantField) = self.value

"""
Write a JLD2 file holding the 6-element `[bx, by, bz, ex, ey, ez]` vector of
callables that `TraceParticlesProblem` expects, and return its path.
"""
function write_test_emfield(dir; b=(0.0, 0.0, 1e-4), e=(0.0, 0.0, 0.0))
    path = joinpath(dir, "emfield.jld2")
    save_object(path, [ConstantField.(b)..., ConstantField.(e)...])
    return path
end

"""
Write the temperature and density files an `MHDSampler` needs, and return an
`SampleICsFromMHD` spec pointing at them.
"""
function write_test_sample(dir; bounds=(0.0, 1.0), t0bounds=(0.0, 0.0))
    tg_file = joinpath(dir, "tg.jld2")
    target_file = joinpath(dir, "r.jld2")
    save_object(tg_file, ConstantField(1e6))      # temperature [K]
    save_object(target_file, ConstantField(1e-10)) # mass density
    return SampleICsFromMHD(
        tg_file = tg_file,
        target_distr_file = target_file,
        xbounds = bounds,
        ybounds = bounds,
        zbounds = bounds,
        t0bounds = t0bounds,
    )
end

"""
Write an HDF5 file of pre-calculated initial conditions for `npart` particles,
in the layout `create_diffeq_ic` reads, and return its path.
"""
function write_test_icfile(dir, npart; eomid=nothing)
    path = joinpath(dir, "ic.h5")
    h5open(path, "w") do fid
        fid["x0"] = collect(range(0.0, 1.0, npart))
        fid["y0"] = zeros(npart)
        fid["z0"] = zeros(npart)
        fid["vx0"] = ones(npart)
        fid["vy0"] = zeros(npart)
        fid["vz0"] = ones(npart)
        fid["t0"] = zeros(npart)
        fid["charge"] = fill(-TP.e, npart)
        fid["mass"] = fill(TP.m_e, npart)
        fid["magneticmoment"] = zeros(npart)
        fid["weight"] = ones(npart)
        fid["nrejections"] = zeros(Int, npart)
        if !isnothing(eomid)
            fid["eomid"] = fill(eomid, npart)
        end
    end
    return path
end

@testset verbose = verbose ≥ 4 "traceparticlesproblem.jl" begin

    # Keywords shared by the helper tests below. These never touch the disk;
    # only the end-to-end testset writes files.
    baseparams = (;
        npart = 10,
        tf = 1.0,
        mass = TP.m_e,
        charge = -TP.e,
        alg = Tsit5(),
        maxiters = 100,
        reltol = 1e-8,
        abstol = 1e-8,
        expname = "unittest",
        batchsize = 5,
        emfield_file = "emfield.jld2",
        prob_func = ICsFromFile(ic_file="ic.h5"),
        verbose = false,
        log2file=false,
    )
    makeparams(; kwargs...) = TraceParticlesParameters(;
        baseparams..., eom=lorentzforce!, kwargs...
    )

    # A stand-in for a loaded interpolator; only its identity matters here.
    dummy_itp = ConstantField(1.0)
    dummy_field = ElectromagneticFieldInterpolator(
        [ConstantField(0.0) for _ in 1:6]
    )

    @testset verbose = verbose ≥ 5 "construct_itpwrapper" begin
        # The four flag combinations map to four different wrappers.
        plain = construct_itpwrapper(makeparams())
        @test plain(dummy_itp) === dummy_itp

        xz = construct_itpwrapper(makeparams(xz_itp_wrapper=true))
        @test xz(dummy_itp) isa XZInterpolation

        static = construct_itpwrapper(makeparams(static_itp_wrapper=true))
        @test static(dummy_itp) isa StaticInterpolation

        staticxz = construct_itpwrapper(
            makeparams(static_itp_wrapper=true, xz_itp_wrapper=true)
        )
        @test staticxz(dummy_itp) isa StaticXZInterpolation
    end

    @testset verbose = verbose ≥ 5 "create_callbacks" begin
        # No callback specs and a non-hybrid equation of motion: no callbacks.
        # Adding a criterion never has to touch this function.
        @test isempty(create_callbacks(makeparams()))

        # One callback per spec, in the order given, carrying that spec's
        # `save_positions`.
        cbs = create_callbacks(makeparams(callback_specs=(
            OutOfBounds(x=(0.0, 1.0), save_positions=(true, false)),
        )))
        @test length(cbs) == 1
        @test cbs[1] isa DiscreteCallback
        @test cbs[1].save_positions == [true, false]

        @test length(create_callbacks(makeparams(
            callback_specs=(Relativistic(),)
        ))) == 1
        @test length(create_callbacks(makeparams(
            callback_specs=(HighMagneticGradient(),)
        ))) == 1
        @test length(create_callbacks(makeparams(callback_specs=(
            OutOfBounds(x=(0.0, 1.0)),
            Relativistic(),
            HighMagneticGradient(),
        )))) == 3

        # A hybrid run always gets the switch callback on top of the specs.
        @test length(create_callbacks(makeparams(eom=hybridgcafo!))) == 1
        @test length(create_callbacks(
            makeparams(eom=hybridgcafo!, callback_specs=(HighMagneticGradient(),))
        )) == 2
        # The switch callback honours `switch_save_positions`.
        switchcb = only(create_callbacks(
            makeparams(eom=hybridgcafo!, switch_save_positions=(false, true))
        ))
        @test switchcb.save_positions == [false, true]

        # The relativistic condition needs a known equation of motion.
        notaneom!(du, u, p, t) = nothing
        @test_throws ArgumentError create_callbacks(
            makeparams(eom=notaneom!, callback_specs=(Relativistic(),))
        )
    end

    @testset verbose = verbose ≥ 5 "SampleICsFromMHD" begin
        mktempdir() do dir

            spec = write_test_sample(dir)
            @test spec isa SampleICsFromMHD
            # No importance sampling by default.
            @test isempty(spec.proposal_distr_file)
            @test spec.t0bounds == (0.0, 0.0)

            # Bounds are mandatory keywords; a half-specified sampling setup
            # cannot be constructed.
            @test_throws UndefKeywordError SampleICsFromMHD(
                tg_file="tg.jld2", target_distr_file="r.jld2",
                xbounds=(0.0, 1.0), ybounds=(0.0, 1.0),
            )
            @test_throws UndefKeywordError SampleICsFromMHD(
                target_distr_file="r.jld2",
                xbounds=(0.0, 1.0), ybounds=(0.0, 1.0), zbounds=(0.0, 1.0),
            )
            # Inverted bounds are rejected.
            @test_throws ArgumentError SampleICsFromMHD(
                tg_file="tg.jld2", target_distr_file="r.jld2",
                xbounds=(1.0, 0.0), ybounds=(0.0, 1.0), zbounds=(0.0, 1.0),
            )

            #= The bounds are promoted among themselves, and that type is what
            the sampler draws positions in =#
            mixed = SampleICsFromMHD(
                tg_file="tg.jld2", target_distr_file="r.jld2",
                xbounds=(0, 1), ybounds=(0.0, 1.0), zbounds=(0.0, 1.0),
            )
            @test mixed isa SampleICsFromMHD
            @test SampleICsFromMHD(
                tg_file="tg.jld2", target_distr_file="r.jld2",
                xbounds=(0.0f0, 1.0f0), ybounds=(0.0f0, 1.0f0),
                zbounds=(0.0f0, 1.0f0), t0bounds=(0.0f0, 0.0f0),
            ) isa SampleICsFromMHD

            # Resolving loads the interpolators and the proposal maximum. With
            # no proposal file the target doubles as the proposal.

            target_distr = load_object(spec.target_distr_file)

            sampler = MHDSampler(
                target_distr,
                target_distr,
                load_object(spec.tg_file),
                spec.xbounds[1], spec.xbounds[2],
                spec.ybounds[1], spec.ybounds[2],
                spec.zbounds[1], spec.zbounds[2],
                spec.t0bounds[1], spec.t0bounds[2],
                maximum(target_distr),
                spec.precision
            )
            @test sampler.proposal_distr === sampler.target_distr
            @test sampler.max_value == 1e-10
            @test sampler.xmin === 0.0
            @test sampler.xmax === 1.0
            @test sampler.precision === Float64

            # Calling it draws a position, a velocity and a weight.
            x, y, z, vx, vy, vz, t, weight, nrejections =
                sampler(TP.m_e, Xoshiro(1))
            @test 0.0 ≤ x ≤ 1.0
            @test 0.0 ≤ y ≤ 1.0
            @test 0.0 ≤ z ≤ 1.0
            @test t == 0.0
            # Target used as its own proposal means unit importance weights.
            @test weight ≈ 1.0
            @test nrejections isa Int
            @test all(isfinite, (vx, vy, vz))
        end
    end

    @testset verbose = verbose ≥ 5 "init_probfunc" begin
        # Anything unrecognised is passed straight through, so a user-supplied
        # problem function needs no method.
        custom = (prob, ctx) -> prob
        @test init_probfunc(custom, makeparams(), identity, dummy_field) === custom

        mktempdir() do dir
            spec = write_test_sample(dir)

            # The equations of motion pick the problem function.
            fo = init_probfunc(
                spec, makeparams(eom=lorentzforce!), identity, dummy_field
            )
            @test fo isa SampleFullOrbit
            @test fo.tf == 1.0
            @test fo.sampler isa MHDSampler

            @test init_probfunc(
                spec, makeparams(eom=guidingcentreapproximation!),
                identity, dummy_field,
            ) isa SampleGCA

            hyb = init_probfunc(
                spec,
                makeparams(eom=hybridgcafo!, initialeomid=EoMID.FullOrbit),
                identity, dummy_field,
            )
            @test hyb isa SampleHybrid
            @test hyb.initialeomid == EoMID.FullOrbit
            # `save_switchinfo` is the single source of truth for the scheme.
            @test hyb.hybridscheme == HybridScheme.Default
            @test init_probfunc(
                spec, makeparams(eom=hybridgcafo!, save_switchinfo=true),
                identity, dummy_field,
            ).hybridscheme == HybridScheme.SaveSwitches

            notaneom!(du, u, p, t) = nothing
            @test_throws ArgumentError init_probfunc(
                spec, makeparams(eom=notaneom!), identity, dummy_field
            )
        end
    end

    @testset verbose = verbose ≥ 5 "determine_output_func" begin
        @test determine_output_func(makeparams(eom=lorentzforce!)) ===
            output_func_lightweight_fo
        @test determine_output_func(makeparams(eom=guidingcentreapproximation!)) ===
            output_func_lightweight
        @test determine_output_func(
            makeparams(eom=guidingcentreapproximation!, save_max_observables=true)
        ) === output_func_max_lightweight
        @test determine_output_func(makeparams(eom=hybridgcafo!)) ===
            output_func_lightweight_hybrid
        # Saving switch info keeps the whole solution, so it uses a closure
        # rather than one of the named lightweight output functions.
        withswitches = determine_output_func(
            makeparams(eom=hybridgcafo!, save_switchinfo=true)
        )
        @test withswitches ∉ (
            output_func_lightweight_hybrid,
            output_func_lightweight,
            output_func_lightweight_fo,
        )

        # Test custom output_function
        custom = (sol, ctx) -> (sol, false)
        @test determine_output_func(makeparams(output_func = custom)) === custom
    end

    @testset verbose = verbose ≥ 5 "determine_reduction" begin
        red = determine_reduction(makeparams(npart=10, batchsize=5), 5)
        @test red isa SaveBatchAsHDF5
        @test red.batchsize == 5
        @test red.nbatches == 2
        @test red.expname == "unittest"
        @test !red.verbose

        # A final partial batch is rounded up, not truncated.
        @test determine_reduction(makeparams(npart=11, batchsize=5), 5).nbatches == 3

        # Saving switch info bypasses the HDF5 writer.
        @test !(determine_reduction(
            makeparams(eom=hybridgcafo!, save_switchinfo=true), 5
        ) isa SaveBatchAsHDF5)

        # Test custom output_function
        custom = (u, data, I) -> (append!(u, data), false)
        @test determine_reduction(makeparams(reduction = custom), 1) === custom
    end

    @testset verbose = verbose ≥ 5 "define_parameters" begin
        @test define_parameters(makeparams(eom=lorentzforce!), dummy_field) isa
            FullOrbitParams
        gca = define_parameters(
            makeparams(eom=guidingcentreapproximation!), dummy_field
        )
        @test gca isa GCAParams
        @test gca.charge == -TP.e
        @test gca.mass == TP.m_e
        @test gca.electromagneticfield === dummy_field
        @test gca.magneticmoment == 0.0

        hybrid = define_parameters(
            makeparams(eom=hybridgcafo!, initialeomid=EoMID.GuidingCentreApproximation),
            dummy_field,
        )
        @test hybrid isa HybridParams
        @test hybrid.eomid == EoMID.GuidingCentreApproximation
        @test define_parameters(
            makeparams(eom=hybridgcafo!, save_switchinfo=true), dummy_field
        ) isa HybridParamsWithDetection

        #= The template magnetic moment follows `mass`, so a single-precision
        mass is not promoted back to `Float64` by the container. =#
        gca32 = define_parameters(
            makeparams(
                eom=guidingcentreapproximation!,
                mass=Float32(TP.m_e),
                charge=Float32(-TP.e),
            ),
            dummy_field,
        )
        @test gca32.magneticmoment === 0.0f0
        @test gca32.mass isa Float32

        notaneom!(du, u, p, t) = nothing
        @test_throws ArgumentError define_parameters(
            makeparams(eom=notaneom!), dummy_field
        )
    end

    @testset verbose = verbose ≥ 5 "define_odeproblem" begin
        prob = define_odeproblem(makeparams(eom=lorentzforce!), dummy_field)
        @test prob isa ODEProblem
        @test prob.f.f === lorentzforce!
        @test eltype(prob.u0) === Float64
        @test prob.tspan == (0.0, 1.0)
        @test prob.p isa FullOrbitParams

        # `u0` follows `mass` and the time span follows `tf`, so consistent
        # single-precision input gives a single-precision template.
        prob32 = define_odeproblem(
            makeparams(
                eom=lorentzforce!,
                tf=1.0f0,
                mass=Float32(TP.m_e),
                charge=Float32(-TP.e),
            ),
            dummy_field,
        )
        @test eltype(prob32.u0) === Float32
        @test prob32.tspan === (0.0f0, 1.0f0)
    end

    @testset verbose = verbose ≥ 5 "TraceParticlesProblem assembly" begin
        mktempdir() do dir
            npart = 6
            batchsize = 4
            params = TraceParticlesParameters(;
                baseparams...,
                eom = guidingcentreapproximation!,
                npart = npart,
                batchsize = batchsize,
                datadir = dir,
                emfield_file = write_test_emfield(dir),
                prob_func = ICsFromFile(ic_file=write_test_icfile(dir, npart)),
                callback_specs = (OutOfBounds(x=(-1.0, 2.0)),),
            )
            tpprob = TraceParticlesProblem(params; logger=current_logger())

            @test tpprob isa TraceParticlesProblem
            @test tpprob.trajectories == npart
            # `batchsize` is capped at `npart`; here it is already smaller.
            @test tpprob.batch_size == batchsize
            @test tpprob.ensemble_prob isa EnsembleProblem
            @test tpprob.callbackset isa CallbackSet
            @test length(tpprob.callbackset.discrete_callbacks) == 1

            # Pre-calculated initial conditions resolve into a `PredefinedICs`
            # problem function holding one entry per particle.
            @test tpprob.ensemble_prob.prob_func isa PredefinedICs
            @test length(tpprob.ensemble_prob.prob_func.u0) == npart
            @test length(tpprob.ensemble_prob.prob_func.params) == npart
            @test tpprob.ensemble_prob.reduction isa SaveBatchAsHDF5

            # The field survives the JLD2 round-trip and the wrapping, and
            # returns (E, B) at a point.
            emfield = tpprob.ensemble_prob.prob.p.electromagneticfield
            @test emfield isa ElectromagneticFieldInterpolator
            efield, bfield = emfield(0.5, 0.0, 0.0, 0.0)
            @test bfield ≈ [0.0, 0.0, 1e-4]
            @test efield ≈ [0.0, 0.0, 0.0]

            # The experiment directory and its log file are created.
            @test isdir(joinpath(dir, "unittest"))
        end

        # The sampling path resolves its filenames at assembly, so the
        # parameters could be constructed before the files existed.
        mktempdir() do dir
            npart = 4
            params = TraceParticlesParameters(;
                baseparams...,
                eom = guidingcentreapproximation!,
                npart = npart,
                batchsize = npart,
                datadir = dir,
                emfield_file = write_test_emfield(dir),
                prob_func = write_test_sample(dir),
            )
            tpprob = TraceParticlesProblem(params; logger=current_logger())
            @test tpprob.ensemble_prob.prob_func isa SampleGCA
            @test tpprob.ensemble_prob.prob_func.sampler isa MHDSampler
        end

        # `batchsize` larger than `npart` is capped at `npart`.
        mktempdir() do dir
            npart = 3
            params = TraceParticlesParameters(;
                baseparams...,
                eom = guidingcentreapproximation!,
                npart = npart,
                batchsize = 100,
                datadir = dir,
                emfield_file = write_test_emfield(dir),
                prob_func = ICsFromFile(ic_file=write_test_icfile(dir, npart)),
            )
            @test TraceParticlesProblem(
                    params; logger=current_logger()
                ).batch_size == npart

        end
    end

end # testset traceparticlesproblem.jl

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

using TraceParticles:
    electron_mass,
    electron_charge,
    EoMID,
    FullOrbitParams,
    GCAParams,
    HighMagneticGradient,
    HybridParams,
    HybridParamsWithDetection,
    HybridScheme,
    ICsFromFile,
    MHDSampler,
    OutOfBounds,
    PredefinedICs,
    Relativistic,
    RelativisticConditionHybrid,
    SampleFullOrbit,
    SampleGCA,
    SampleHybrid,
    SaveBatchAsHDF5,
    TraceParticlesParameters,
    TraceParticlesProblem,
    create_callbacks,
    define_odeproblem,
    define_parameters,
    determine_output_func,
    determine_reduction,
    guidingcentreapproximation!,
    hybridgcafo!,
    init_probfunc,
    lorentzforce!,
    output_func_lightweight_gca,
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

"""
Write a JLD2 file holding the 6-element `[bx, by, bz, ex, ey, ez]` vector of
callables that `TraceParticlesProblem` expects, and return its path.
"""
function write_test_emfield(dir; b=(0.0, 0.0, 1e-4), e=(0.0, 0.0, 0.0))
    path = joinpath(dir, "emfield.jld2")
    save_object(path, [ConstantField.(b)..., ConstantField.(e)...])
    return path
end

function fixture_mhdsampler(; lb=0.0, ub=1.0, tmin=0.0, tmax=0.0)
    temperature = ConstantField(1e6) # K
    target = ConstantField(1e-10) # kg/m^3
    return MHDSampler(;
        temperature=temperature,
        target_distr=target,
        xmin=lb, ymin=lb, zmin=lb, tmin=tmin,
        xmax=ub, ymax=ub, zmax=ub, tmax=tmax,
        maxval=1e-10,
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
        fid["charge"] = fill(electron_charge, npart)
        fid["mass"] = fill(electron_mass, npart)
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
    dummy_field(x, y, z, t) = [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]

    # Keywords shared by the helper tests below. These never touch the disk;
    # only the end-to-end testset writes files.
    baseparams = (;
        npart = 10,
        tf = 1.0,
        mass = electron_mass,
        charge = electron_charge,
        alg = Tsit5(),
        maxiters = 100,
        reltol = 1e-8,
        abstol = 1e-8,
        expname = "test",
        batchsize = 5,
        electromagneticfield = dummy_field,
        prob_func = ICsFromFile(ic_file="ic.h5"),
    )
    makeparams(; kwargs...) = TraceParticlesParameters(;
        baseparams..., eom=lorentzforce!, kwargs...
    )

    # A stand-in for a loaded interpolator; only its identity matters here.

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
        @test length(create_callbacks(makeparams(callback_specs=(
            OutOfBounds(x=(0.0, 1.0)),
            Relativistic(),
            HighMagneticGradient(),
        )))) == 3

        # A hybrid run always gets the switch callback on top of the specs.
        @test length(create_callbacks(makeparams(eom=hybridgcafo!))) == 1
        cbs = create_callbacks(
            makeparams(
                eom=hybridgcafo!,
                switch_save_positions=(false, true),
                callback_specs=(Relativistic(),)
            )
        )
        @test cbs[1].condition isa RelativisticConditionHybrid
        # The switch callback honours `switch_save_positions`.
        @test cbs[2].save_positions == [false, true]

        # The relativistic condition needs a known equation of motion.
        notaneom!(args...) = nothing
        @test_throws ArgumentError create_callbacks(
            makeparams(eom=notaneom!, callback_specs=(Relativistic(),))
        )
    end

    @testset verbose = verbose ≥ 5 "determine_output_func" begin
        @test determine_output_func(makeparams(eom=lorentzforce!)) ===
            output_func_lightweight_fo
        @test determine_output_func(makeparams(eom=guidingcentreapproximation!)) ===
            output_func_lightweight_gca
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
            output_func_lightweight_gca,
            output_func_lightweight_fo,
        )

        # Test custom output_function
        custom = (sol, ctx) -> (sol, false)
        @test determine_output_func(makeparams(output_func = custom)) === custom
    end

    @testset verbose = verbose ≥ 5 "determine_reduction" begin
        red, effective_batchsize = determine_reduction(
            makeparams(npart=10, batchsize=5, verbose=false)
        )
        @test effective_batchsize == 5
        @test red isa SaveBatchAsHDF5
        @test red.batchsize == 5
        @test red.nbatches == 2
        @test red.expname == "test"
        @test !red.verbose

        # A final partial batch is rounded up, not truncated.
        @test determine_reduction(makeparams(npart=11, batchsize=5))[1].nbatches == 3

        # Saving switch info bypasses the HDF5 writer.
        @test !(determine_reduction(
            makeparams(eom=hybridgcafo!, save_switchinfo=true)
        ) isa SaveBatchAsHDF5)

        # Test custom output_function
        custom = (u, data, I) -> (append!(u, data), false)
        @test determine_reduction(makeparams(reduction = custom))[1] === custom
    end

    @testset verbose = verbose ≥ 5 "define_parameters" begin
        @test define_parameters(makeparams(eom=lorentzforce!), dummy_field) isa
            FullOrbitParams
        gca = define_parameters(
            makeparams(eom=guidingcentreapproximation!), dummy_field
        )
        @test gca isa GCAParams
        @test gca.charge == makeparams().charge
        @test gca.mass == makeparams().mass
        @test gca.electromagneticfield === dummy_field

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
                mass=Float32(electron_mass),
                charge=Float32(electron_charge),
            ),
            dummy_field,
        )
        @test gca32.magneticmoment isa Float32
        @test gca32.mass isa Float32

        notaneom!(du, u, p, t) = nothing
        @test_throws ArgumentError define_parameters(
            makeparams(eom=notaneom!), dummy_field
        )
    end

    @testset verbose = verbose ≥ 5 "define_odeproblem" begin
        # `u0` follows `mass` and the time span follows `tf`, so consistent
        # single-precision input gives a single-precision template.
        prob32 = define_odeproblem(
            makeparams(
                eom=lorentzforce!,
                tf=1.0f0,
                mass=Float32(electron_mass),
                charge=Float32(electron_charge),
            ),
            dummy_field,
        )
        @test prob32.p isa FullOrbitParams
        @test prob32 isa ODEProblem
        @test prob32.f.f === lorentzforce!
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
                electromagneticfield = dummy_field,
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
            @test tpprob.ensemble_prob.prob.p.electromagneticfield ==
                dummy_field
            # The experiment directory and its log file are created.
            @test isdir(joinpath(dir, "test"))
        end

        # The sampling path resolves its filenames at assembly, so the
        # parameters could be constructed before the files existed.
        mktempdir() do dir
            npart = 4
            params = TraceParticlesParameters(;
                baseparams...,
                eom = guidingcentreapproximation!,
                npart = npart,
                batchsize = 100,
                datadir = dir,
                electromagneticfield = joinpath(dir, "emfield.jl"),
                prob_func = joinpath(dir, "prob_func.jl"),
                log2file = true
            )

            emfield_script = raw"""
            struct EMField{T}
                E::T
                B::T
            end
            (self::EMField)(x, y, z, t) = self.E, self.B

            myfield = EMField([0,1,0], [0,0,1])
            """
            write(joinpath(dir, "emfield.jl"), emfield_script)
            prob_func_script = raw"""
            using TraceParticles: MHDSampler

            tg(x, y, z, t) = 1e6
            density(x, y, z, t) = 1e-10

            MHDSampler(
                temperature = tg,
                target_distr = density,
                xmin = 0, xmax = 1.0,
                ymin = 0, ymax = 1,
                zmin = 0, zmax = 1,
                tmin = 0, tmax = 1,
                maxval = 1e-10
            )
            """
            write(joinpath(dir, "prob_func.jl"), prob_func_script)

            tpprob = TraceParticlesProblem(params)
            invokelatest() do
                @test tpprob.batch_size == npart
                @test tpprob.ensemble_prob.prob_func isa SampleGCA
                @test tpprob.ensemble_prob.prob_func.sampler isa MHDSampler
                @test tpprob.ensemble_prob.prob_func.sampler.target_distr ===
                    tpprob.ensemble_prob.prob_func.sampler.proposal_distr 
                @test tpprob.ensemble_prob.prob_func.sampler.zmin === 0.0
                @test tpprob.ensemble_prob.prob.p.electromagneticfield isa
                    parentmodule(
                        typeof(tpprob.ensemble_prob.prob.p.electromagneticfield)
                    ).EMField
                @test tpprob.ensemble_prob.prob.p.
                    electromagneticfield(0,0,0,0) == ([0,1,0], [0,0,1])
                @test isfile(joinpath(dir, "test", "test.log.txt"))
            end
        end
    end
end # testset traceparticlesproblem.jl

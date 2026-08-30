using Logging
using OrdinaryDiffEq
using Test

using TraceParticles:
    EoMID,
    HighMagneticGradient,
    ICsFromFile,
    SampleICsFromMHD,
    OutOfBounds,
    Relativistic,
    TraceParticlesParameters,
    guidingcentreapproximation!,
    hybridgcafo!,
    lorentzforce!,
    setparameters,
    validate

if !isdefined(Main, :verbose)
    verbose = 4
end

@testset verbose = verbose ≥ 4 "traceparticlesparameters.jl" begin

    # A minimal, valid set of keywords. Uses pre-calculated initial conditions
    # so that no sampling inputs are needed; nothing here touches the disk,
    # because the filenames are only resolved when a problem is assembled.
    baseparams = (;
        eom = lorentzforce!,
        npart = 10,
        tf = 1.0,
        mass = 1.0,
        charge = -1.0,
        alg = Tsit5(),
        maxiters = 100,
        reltol = 1e-8,
        abstol = 1e-8,
        expname = "unittest",
        emfield_file = "emfield.jld2",
        prob_func = ICsFromFile(ic_file="ic.h5"),
    )

    @testset verbose = verbose ≥ 5 "Field types and defaults" begin
        p = TraceParticlesParameters(; baseparams...)
        @test p isa TraceParticlesParameters
        @test p.npart isa Int
        @test p.maxiters isa Int
        @test p.seed isa Int
        @test p.expname isa String
        @test p.datadir isa String
        @test p.loglevel isa Logging.LogLevel
        @test p.initialeomid isa EoMID.T
        @test p.switch_save_positions isa Tuple{Bool,Bool}

        @test p.datadir == "data"
        @test p.seed == 1
        @test p.verbose
        @test p.initialeomid == EoMID.FullOrbit
        @test !p.save_switchinfo
        # No callbacks unless asked for.
        @test p.callback_specs == ()
        # `batchsize` defaults to the whole ensemble.
        @test p.batchsize == p.npart == 10
        @test TraceParticlesParameters(; baseparams..., batchsize=3).batchsize == 3

        # Mandatory fields have no default.
        for missingfield in (:eom, :npart, :prob_func, :emfield_file, :tf, :alg)
            @test_throws UndefKeywordError TraceParticlesParameters(;
                Base.structdiff(baseparams, NamedTuple{(missingfield,)}((nothing,)))...
            )
        end
    end

    @testset verbose = verbose ≥ 5 "Precision is carried, not converted" begin
        p = TraceParticlesParameters(; baseparams...)
        for name in (:tf, :mass, :charge, :reltol, :abstol)
            @test typeof(getfield(p, name)) === Float64
        end

        p32 = TraceParticlesParameters(;
            baseparams...,
            tf = 1.0f0, mass = 1.0f0, charge = -1.0f0,
            reltol = 1.0f-8, abstol = 1.0f-8,
        )
        for name in (:tf, :mass, :charge, :reltol, :abstol)
            @test typeof(getfield(p32, name)) === Float32
        end

        mixed = TraceParticlesParameters(; baseparams..., mass=1.0f0)
        @test mixed.mass === 1.0f0
        @test mixed.gradient_tolerance === 1e-4
    end

    @testset verbose = verbose ≥ 5 "validate" begin
        @test validate(TraceParticlesParameters(; baseparams...)) === nothing

        @test_throws ArgumentError TraceParticlesParameters(; baseparams..., npart=0)
        @test_throws ArgumentError TraceParticlesParameters(; baseparams..., npart=-1)
        @test_throws ArgumentError TraceParticlesParameters(; baseparams..., batchsize=0)
        @test_throws ArgumentError TraceParticlesParameters(; baseparams..., tf=Inf)
        @test_throws ArgumentError TraceParticlesParameters(; baseparams..., tf=NaN)

        # `callback_specs` must hold callback specifications.
        @test_throws ArgumentError TraceParticlesParameters(;
            baseparams..., callback_specs=(OutOfBounds(x=(0.0, 1.0)), "not a spec")
        )
        @test TraceParticlesParameters(;
            baseparams...,
            callback_specs=(
                OutOfBounds(x=(0.0, 1.0)),
                Relativistic(fraction=0.02),
                HighMagneticGradient(tolerance=1e-3),
            ),
        ) isa TraceParticlesParameters
    end

    @testset verbose = verbose ≥ 5 "prob_func" begin
        # A prob_func specification is stored as given; it is resolved only when
        # a problem is assembled.
        p = TraceParticlesParameters(; baseparams...)
        @test p.prob_func isa ICsFromFile
        @test p.prob_func.ic_file == "ic.h5"

        sampled = TraceParticlesParameters(;
            baseparams...,
            prob_func = SampleICsFromMHD(
                tg_file = "tg.jld2",
                target_distr_file = "r.jld2",
                xbounds = (0.0, 1.0),
                ybounds = (2.0, 3.0),
                zbounds = (4.0, 5.0),
            ),
        )
        @test sampled.prob_func isa SampleICsFromMHD

        # Any callable of `(prob, ctx)` is accepted and passed through
        # untouched, so a user-supplied problem function needs no spec type.
        custom = (prob, ctx) -> prob
        @test TraceParticlesParameters(;
            baseparams..., prob_func=custom
        ).prob_func === custom
    end

    @testset verbose = verbose ≥ 5 "output_func" begin
        p = TraceParticlesParameters(; baseparams...)
        @test isnothing(p.output_func)

        custom = (sol, ctx) -> (sol, false)
        @test TraceParticlesParameters(;
            baseparams..., output_func=custom
        ).output_func === custom
    end

    @testset verbose = verbose ≥ 5 "reduction" begin
        p = TraceParticlesParameters(; baseparams...)
        @test isnothing(p.reduction)

        custom = (u, data, I) -> (append!(u, data), false)
        @test TraceParticlesParameters(;
            baseparams..., reduction=custom
        ).reduction === custom
    end

    @testset verbose = verbose ≥ 5 "setparameters" begin
        p = TraceParticlesParameters(; baseparams..., npart=10, batchsize=5)

        # Only the named fields change; everything else is carried over.
        varied = setparameters(p; npart=20, eom=guidingcentreapproximation!)
        @test varied.npart == 20
        @test varied.eom === guidingcentreapproximation!
        @test varied.batchsize == 5
        @test varied.expname == p.expname
        @test varied.prob_func === p.prob_func
        # The original is left alone, so it can be varied again.
        @test p.npart == 10
        @test p.eom === lorentzforce!

        # No changes gives back the same values.
        unchanged = setparameters(p)
        for name in fieldnames(TraceParticlesParameters)
            @test getfield(unchanged, name) == getfield(p, name)
        end

        # A misspelled field is caught here rather than silently ignored.
        @test_throws ArgumentError setparameters(p; nparticles=20)
        # The copy is validated like any other instance.
        @test_throws ArgumentError setparameters(p; npart=0)
    end

    @testset verbose = verbose ≥ 5 "show" begin
        p = TraceParticlesParameters(; baseparams...)
        str = sprint(show, MIME"text/plain"(), p)
        @test occursin("TraceParticlesParameters", str)
        for name in fieldnames(TraceParticlesParameters)
            @test occursin(string(name), str)
        end
    end

end # testset traceparticlesparameters.jl

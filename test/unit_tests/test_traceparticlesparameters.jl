using Logging
using OrdinaryDiffEq
using Test

using TraceParticles:
    EoMID,
    HighMagneticGradient,
    ICsFromFile,
    OutOfBounds,
    Relativistic,
    TraceParticlesParameters,
    guidingcentreapproximation!,
    lorentzforce!,
    setparameters,
    validate

if !isdefined(Main, :verbose)
    verbose = 4
end

@testset verbose = verbose ≥ 4 "traceparticlesparameters.jl" begin

    # A valid set of keywords. Uses pre-calculated initial conditions
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
        expname = "test",
        electromagneticfield = "emfield.jl",
        prob_func = ICsFromFile(ic_file="ic.h5"),
        output_func = (sol, ctx) -> (sol, false),
        reduction = (u, data, I) -> (append!(u, data), false),
        solver_options = (save_everystep=false, saveat=[0.1,0.2])
    )

    @testset verbose = verbose ≥ 5 "Defaults" begin
        p = TraceParticlesParameters(; baseparams...)
        # No callbacks unless asked for.
        @test p.callback_specs == ()
        # `batchsize` defaults to the whole ensemble.
        @test p.batchsize == p.npart == 10
        # Mandatory fields have no default.
        for missingfield in (:eom, :npart, :prob_func, :electromagneticfield, :tf, :alg)
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

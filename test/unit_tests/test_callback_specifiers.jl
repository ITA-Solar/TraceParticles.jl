using OrdinaryDiffEq
using SciMLBase
using Test

using TraceParticles:
    CallbackSpec,
    HighMagneticGradient,
    OutOfBounds,
    Relativistic,
    create_callback,
    guidingcentreapproximation!,
    hybridgcafo!,
    lorentzforce!

if !isdefined(Main, :verbose)
    verbose = 4
end

@testset verbose = verbose ≥ 4 "callback_specifiers.jl" begin

    @testset verbose = verbose ≥ 5 "OutOfBounds" begin
        s = OutOfBounds(x=(0.0, 1.0), z=(-1.0, 2.0))
        @test s isa CallbackSpec
        @test s.u_idxs === (1, 3)
        @test s.bounds == ((0.0, 1.0), (-1.0, 2.0))
        @test s.save_positions == (false, true)

        # Axes are recorded in x, y, z order regardless of keyword order.
        @test OutOfBounds(z=(0.0, 1.0), x=(0.0, 1.0)).u_idxs === (1, 3)
        @test OutOfBounds(y=(0.0, 1.0)).u_idxs === (2,)
        @test OutOfBounds(x=(0.0, 1.0), y=(0.0, 1.0), z=(0.0, 1.0)).u_idxs === (1, 2, 3)

        # Bounds are mandatory: an unbounded `OutOfBounds` cannot be built at
        # all, so this needs no run-time check elsewhere.
        @test_throws ArgumentError OutOfBounds()
        @test_throws ArgumentError OutOfBounds(save_positions=(true, true))
        # Inverted or empty bounds are rejected.
        @test_throws ArgumentError OutOfBounds(x=(1.0, 0.0))
        @test_throws ArgumentError OutOfBounds(x=(1.0, 1.0))

        #= The bounds are promoted among themselves so that mixed input is
        still representable. They are only ever compared against the state
        vector, so this cannot promote the integration. =#
        @test OutOfBounds(x=(0, 1)).bounds == ((0.0, 1.0),)
        @test eltype(OutOfBounds(x=(0, 1)).bounds[1]) === Float64
        @test eltype(OutOfBounds(x=(0.0f0, 1.0f0)).bounds[1]) === Float32
        @test eltype(OutOfBounds(x=(0.0f0, 1.0f0), z=(0, 1)).bounds[1]) === Float32
    end

    @testset verbose = verbose ≥ 5 "Relativistic and HighMagneticGradient" begin
        @test Relativistic() isa CallbackSpec
        @test Relativistic().fraction == 0.12
        @test Relativistic(fraction=0.02).fraction == 0.02
        @test Relativistic().save_positions == (true, false)
        # The given precision is kept rather than converted.
        @test Relativistic(fraction=0.02f0).fraction === 0.02f0

        @test HighMagneticGradient() isa CallbackSpec
        @test HighMagneticGradient().tolerance == 1e-4
        @test HighMagneticGradient(tolerance=1e-3).tolerance == 1e-3
        @test HighMagneticGradient(tolerance=1.0f-3).tolerance === 1.0f-3
    end

    @testset verbose = verbose ≥ 5 "create_callback" begin
        # `create_callback` only reads `mass` and `eom` off the parameters, so
        # a NamedTuple stands in for them here.
        for eom in (lorentzforce!, guidingcentreapproximation!, hybridgcafo!)
            p = (; mass=1.0, eom=eom)
            cb = create_callback(Relativistic(fraction=0.02), p)
            @test cb isa DiscreteCallback
            @test cb.save_positions == [true, false]
        end

        # The relativistic condition is chosen per equation of motion.
        notaneom!(du, u, p, t) = nothing
        @test_throws ArgumentError create_callback(
            Relativistic(), (; mass=1.0, eom=notaneom!)
        )

        p = (; mass=1.0, eom=lorentzforce!)
        oob = create_callback(
            OutOfBounds(x=(0.0, 1.0), save_positions=(true, true)), p
        )
        @test oob isa DiscreteCallback
        @test oob.save_positions == [true, true]
        @test oob.condition.bounds == ((0.0, 1.0),)
        @test oob.condition.u_idxs == [1]

        gradb = create_callback(HighMagneticGradient(tolerance=1e-3), p)
        @test gradb isa DiscreteCallback
        @test gradb.condition.tolerance == 1e-3
    end

end # testset callback_specifiers.jl

#-------------------------------------------------------------------------------
# Created 10.01.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
# Script for running fast tests.
#-------------------------------------------------------------------------------

@testset verbose = verbose ≥ 2 "Fast tests" begin

    @testset verbose = verbose ≥ 3 "Unit tests" begin
        include("test_physics.jl")
        include("test_callbacks.jl")
    end

    @testset verbose = verbose ≥ 3 "Regression tests" begin

        @testset verbose = verbose ≥ 4 "ExB-drift" begin
            include("experiments/ExBdrift.jl")
        end # testset ExB-drift

        @testset verbose = verbose ≥ 4 "∇B-drift" begin
            include("experiments/gradBdrift.jl")
        end # testset ∇B-drift

        @testset verbose = verbose ≥ 4 "Mirroring" begin
            include("experiments/mirroring.jl")
        end # testset mirroring

        @testset verbose = verbose ≥ 4 "Speiser (1965)" begin
            include("experiments/speiser1965.jl")
        end # testset Speiser 1965

        @testset verbose = verbose ≥ 4 "Hybrid switch" begin
            include("experiments/hybridswitch.jl")
        end # testset Hybrid switch

    end # testset Experiments

end # testset All test

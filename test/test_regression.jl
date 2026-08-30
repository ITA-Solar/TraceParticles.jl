@testset verbose = verbose ≥ 2 "Long tests" begin

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

        @testset verbose = verbose ≥ 4 "Speiser (1965) w/ hybrid" begin
            include("experiments/speiser1965_hybrid.jl")
        end # testset Hybrid switch

        @testset verbose = verbose ≥ 4 "2D fan-spine" begin
            include("experiments/2Dfanspine/test_2Dfanspine.jl")
        end # testset 2D coronal jet

        @testset verbose = verbose ≥ 4 "Dipole" begin
            include("experiments/dipoleloop.jl")
        end # testset dipole

    end # testset Experiments

end # testset long tests

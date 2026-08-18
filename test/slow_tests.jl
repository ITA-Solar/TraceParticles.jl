@testset verbose = verbose ≥ 2 "Long tests" begin

    @testset verbose = verbose ≥ 3 "Experiments" begin

        @testset verbose = verbose ≥ 4 "Dipole" begin
            include("experiments/dipoleloop.jl")
        end # testset dipole

        @testset verbose = verbose ≥ 4 "2D fan-spine" begin
            include("experiments/2Dfanspine/test_2Dfanspine.jl")
        end # testset 2D coronal jet

    end # testset Experiments

end # testset long tests

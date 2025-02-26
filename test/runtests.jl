if "fast_tests" in ARGS
    include("Jenkins_short.jl")
else
    include("Jenkins.jl")
end
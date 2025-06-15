using BenchmarkTools

if !isempty(ARGS)
    localARGS = copy(ARGS)
end

resultsname = try
    localARGS[1]
catch
    error("Please provide a benchmarkname as the first argument.")
end

paramsfile = try
    localARGS[2]
catch
    defaultparameters = "params-34cd008-macbook20250605.json"
    @info "Using default parameters: $defaultparameters"
    defaultparameters
end

baselinefile = try
    localARGS[3]
catch
    defaultbaseline = "baseline-34cd008-macbook20250605.json"
    @info "Using default baseline: $defaultbaseline"
    defaultbaseline
end

sourcedir = Base.source_dir()
resultsname = joinpath(sourcedir, resultsname)
baselinefile = joinpath(sourcedir, baselinefile)
paramsfile = joinpath(sourcedir, paramsfile)

include("benchmarksuite.jl")

loadparams!(suite, BenchmarkTools.load(paramsfile)[1], :evals, :samples)

results = run(suite, verbose=true)
BenchmarkTools.save(resultsname * "_trials.json", results)

baseline = BenchmarkTools.load(baselinefile)[1]
judgement = judge(minimum(results), minimum(baseline); time_tolerance=0.05)

open(resultsname * "_judgement.txt", "w") do fid
    write(
        fid,
        """
--------------------------------------------------------------------------------
Machine name: $(gethostname())

Parameters used:
$(paramsfile)

Results compared with:
$(baselinefile)
--------------------------------------------------------------------------------

$(string(judgement))
"""
    )
end

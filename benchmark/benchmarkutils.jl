using BenchmarkTools
using Statistics

function tune_suite!(
    suite::BenchmarkGroup;
    fname=joinpath(Base.source_dir(), "params.json")
)
    tune!(suite)
    BenchmarkTools.save(fname, params(suite))
end

function run_suite(
    suite::BenchmarkGroup;
    fname=joinpath(Base.source_dir(), "results.json"),
    verbose=true
)
    results = run(suite, verbose=verbose)
    BenchmarkTools.save(fname, results)
    return results
end

function judge_suite(
    results::BenchmarkTools.Trial,
    baseline::BenchmarkTools.Trial;
    metric::Function=minimum,
    fname=joinpath(Base.source_dir(), "judgement.json"),
    time_tolerance=0.05
)
    judgement = judge(minimum(results), minimum(baseline); time_tolerance=time_tolerance)
    BenchmarkTools.save(fname, judgement)
    return judgement
end

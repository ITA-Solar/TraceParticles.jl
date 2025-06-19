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
    results,
    baselinefile::String;
    metric::Function=minimum,
    time_tolerance=0.15,
    group=nothing,
    subgroup=nothing
)
    baseline = BenchmarkTools.load(baselinefile)[1]
    if isnothing(group)
        judgement = judge(
            metric(results), metric(baseline);
            time_tolerance=time_tolerance
        )
    else
        if isnothing(subgroup)
            judgement = judge(
                metric(results), metric(baseline[group]);
                time_tolerance=time_tolerance
            )
        else
            judgement = judge(
                metric(results), metric(baseline[group][subgroup]);
                time_tolerance=time_tolerance
            )
        end
    end
end

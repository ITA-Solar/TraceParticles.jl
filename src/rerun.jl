"""
    rerun(idxs, datafile::String, params::TraceParticlesParameters;
          paramsbackup=true, otherparamchanges...)

Re-integrate a chosen set of particles in the electromagnetic field of an
existing experiment, keeping their **full trajectories** instead of the
lightweight endpoint summary the ensemble run stores.

The particles are the ones at index `idxs` in `datafile`, the HDF5-file the
experiment wrote, and `params` are the parameters it was run with. Any of those
parameters can be changed for the rerun by passing it as a keyword argument.

The particle data in `datafile` is assumed to callable with
[`h5_getinitialstate`](@ref) and [`h5_getinitialstate_gca`](@ref).

This is the tool for asking *how* a particle gained its energy: the returned
solutions can be fed to [`get_observable`](@ref) to time-integrate the
guiding-centre work channels (`:parallelpower`, `:fermi`, `:betatron`,
`:polarisationacc`), which endpoint data cannot answer.

#### Example
```julia
sol = rerun(1:100, "run.h5", runparams; alg=Vern9())
```
"""
function rerun(
    idxs::AbstractVector,
    datafile::String,
    params::TraceParticlesParameters;
    paramsbackup=true,
    otherparamchanges...
)
    new_params = setparameters(params;
        prob_func = Rerun(datafile, idxs),
        output_func = SciMLBase.DEFAULT_OUTPUT_FUNC,
        reduction = SciMLBase.DEFAULT_REDUCTION,
        expname = params.expname*"-rerun-"*string(now()),
        paramsbackup = paramsbackup,
        # Only the selected particles are integrated, so the ensemble is as
        # large as the selection, not as the experiment.
        npart = length(idxs),
        batchsize = length(idxs),
        otherparamchanges...
    )

    prob = TraceParticlesProblem(new_params)
    return solve(prob, new_params)
end

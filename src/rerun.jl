"""
    rerun(u0, mu0, nparticles, fname; kwargs)
Reruns test particles in the experiment defined by `fname`. The initial
positions and velocities of the particles are arbitrary, defined in `u0`,
and their magnetic moments are defined in `mu0`. Some parameters of the
experiment may be adjusted through keyword arguments. These are

#### Keyword arguments
- `reltol`
- `abstol`
- `maxiters`
- `alg`: solver algorithm
- `gcatol`: Tolerance for the guiding centre approximation.

#### Example use
```julia
fname = "testparticles.h5"
npart = 10
sol = rerun(find_nonthermals(fname, npart)[1:2]..., npart, fname)
"""
function rerun(
    u0::Vector{<:Vector{<:Real}},
    mu0::Vector{<:Real},
    fname::String
    ;
    reltol=nothing,
    abstol=nothing,
    maxiters=nothing,
    alg=Rosenbrock23(),
    gcatol=nothing,
    relativistic_tol=nothing,
    onlymagneticfield=false,
    emfields=nothing,
    EoM=nothing,
    timespan=nothing,
    parallelisation=EnsembleSerial(),
    kwargs...
)
    nparticles = length(u0)
    expname, _ = splitext(fname)
    include(expname * ".jl")

    if !isnothing(timespan)
        tspans = [timespan for _ in 1:nparticles]
    else
        tspans = [tspan for _ in 1:nparticles]
    end

    for cb in solve_kwargs[:callback].discrete_callbacks
        if cb.condition isa MagneticGradientCondition && !isnothing(gcatol)
            set_tol!(cb.condition, gcatol)
        elseif (
            (cb.condition isa RelativisticConditionGCA) &&
            !isnothing(relativistic_tol)
        )
            set_limit!(cb.condition, mass, relativistic_tol)
        end
    end

    s_kwargs = Dict(
        :reltol => isnothing(reltol) ? solve_kwargs[:reltol] : reltol,
        :abstol => isnothing(abstol) ? solve_kwargs[:abstol] : abstol,
        :maxiters => isnothing(maxiters) ? solve_kwargs[:maxiters] : maxiters,
        :trajectories => nparticles,
        :callback => solve_kwargs[:callback],
    )

    if isnothing(emfields)
        fields = fields_itp
    else
        fields = emfields
    end
    if onlymagneticfield
        fi = Vector{Any}(undef, 6)
        for i in 1:3
            fi[i] = fields[i]
        end
        for i in 4:6
            fi[i] = (x, z) -> 0.0
        end
        fields = fi
    end
    if isnothing(EoM)
        eqs = eom
    else
        eqs = EoM
    end

    return rerun(
        u0,
        mu0,
        eqs,
        charge,
        mass,
        fields,
        tspans,
        merge(
            Dict(:safetycopy => false),
            Dict{Symbol,Any}(k => v for (k, v) in kwargs)
        ),
        s_kwargs;
        alg=alg,
        parallelisation=parallelisation
    )
end

"""
    rerun(
        u0::Vector{<:Vector{<:Real}},
        mu0::Vector{<:Real},
        charge::Real,
        mass::Real,
        fields::Vector{Any},
        tspans::Vector{<:Tuple{Real, Real}},
        ensemble_prob_kwargs::Dict{Symbol, Any},
        solve_kwargs::Dict{Symbol, Any};
        alg=Rosenbrock23(),
        parallelisation=EnsembleSerial(),
    )
Rerun GCA particles with initial condition `u0` and magnetic moments `mu0`,
using the parameters `charge`, `mass`, `fields`, and `tspans`. The arguments
`ensemble_prob_kwargs` and `solve_kwargs` are dictionaries containing keyword
arguments for the `EnsembleProblem` and `DifferentialEquations.solve` function.
"""
function rerun(
    u0::Vector{<:Vector{<:Real}},
    mu0::Vector{<:Real},
    eom::Function,
    charge::Real,
    mass::Real,
    electromagneticfield::Any,
    tspans::Vector{<:Tuple{Real,Real}},
    ensemble_prob_kwargs::Dict{<:Symbol,<:Any},
    solve_kwargs::Dict{<:Symbol,<:Any};
    alg=Rosenbrock23(),
    parallelisation=EnsembleSerial(),
)
    nparticles = length(u0)

    prob = ODEProblem(
        eom,
        u0[1],
        tspans[1],
        (
            charge=charge,
            mass=mass,
            electromagneticfield=electromagneticfield,
            magneticmoment=mu0[1],
        )
    )
    if eom == guidingcentreapproximation!
        paramstruct = GCAParams
    elseif eom == hybridgcafo!
        paramstruct = HybridParams
    else
        @error "Don't know which parameter-struct to use."
    end
    ensemble_prob = EnsembleProblem(
        prob
        ;
        prob_func=PposPvel(u0, mu0, tspans, paramstruct),
        ensemble_prob_kwargs...
    )

    @info "Re-running with $nparticles particles..."
    @time sol = DifferentialEquations.solve(
        ensemble_prob,
        alg,
        parallelisation
        ;
        solve_kwargs...
    )
    return sol
end

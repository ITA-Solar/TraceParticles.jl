"""
    rerun(u0, mu0, nparticles, fname; kwargs)
Reruns test particles in the experiment defined by `fname`. The initial
positions and velocities of the particles are arbitrary, defined in `u0`,
and their magnetic moments are defined in `mu0`. Some parameters of the
experiment may be adjusted through keyword arguments. These are

# Keyword arguments
- `reltol`
- `abstol`
- `maxiters`
- `alg`: solver algorithm
- `gcatol`: Tolerance for the guiding centre approximation.
"""
function rerun(
    u0::Vector{<:Vector{<:Real}},
    mu0::Vector{<:Real},
    nparticles::Int,
    fname::String;
    reltol=nothing,
    abstol=nothing,
    maxiters=nothing,
    alg=Tsit5(),
    gcatol=nothing,
    onlymagneticfield=false,
    emfields=nothing,
    EoM=nothing
)
    include(fname * ".jl")
    tspans = [tspan for _ in 1:nparticles]
    if !isnothing(gcatol)
        solve_kwargs[:callback] = CallbackSet(
            solve_kwargs[:callback].discrete_callbacks[1],
            DiscreteCallback(
                GCABreakDownCondition_2Dxz(gcatol),
                gcabreakdownaffect!
            )
        )
    end
    s_kwargs = (
        reltol=isnothing(reltol) ? solve_kwargs[:reltol] : reltol,
        abstol=isnothing(abstol) ? solve_kwargs[:abstol] : abstol,
        maxiters=isnothing(maxiters) ? solve_kwargs[:maxiters] : maxiters,
        trajectories=nparticles,
        callback=solve_kwargs[:callback],
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
    prob = ODEProblem(
        eqs,
        u0[1],
        tspans[1],
        (
            charge=charge,
            mass=mass,
            fields=fields,
            magneticmoment=mu0[1],
        )
    )
    ensemble_prob = EnsembleProblem(
        prob
        ;
        prob_func=PposPvel(u0, mu0, tspans),
        # Use default output function. I.e. return the whole solution
        safetycopy=false,
    )

    solve_args = (alg, EnsembleSerial())
    @info "Re-running $fname with $nparticles particles..."
    @time sol = DifferentialEquations.solve(
        ensemble_prob,
        solve_args...
        ;
        s_kwargs...
    )
    return sol
end

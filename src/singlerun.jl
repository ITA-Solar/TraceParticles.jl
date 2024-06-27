using Distributed
using DifferentialEquations

ensemble_algorithm = EnsembleSerial()
solve_args = try (alg, ensemble_algorithm)
catch
    (ensemble_algorithm,)
end

try u0
catch
    error("`u0::Vector{Vector{<:Real}}` must be defined.")
end
try mu0
catch
    error("`mu0::Vector{<:Real}` must be defined.")
end

# Adjust keyword arguments to DifferentialEquations.solve
npart = length(u0)
solve_kwargs_internal = (solve_kwargs...,
    trajectories = npart,
    batch_size = npart,
    reltol = try reltol
    catch
        try solve_kwargs.reltol
        catch 
            nothing
        end
    end,
    abstol = try abstol
    catch
        try solve_kwargs.abstol
        catch 
            nothing
        end
    end,
    maxiters = try maxiters
    catch
        try solve_kwargs.maxiters
        catch 
            nothing
        end
    end
    )

# Adjust keyword arguments to EnsembleProb
ensembleprob_kwargs_internal = (ensembleprob_kwargs...,
    prob_func = tp.PposPvel(u0, mu0),
    # Use default output function. I.e. return the whole solution
    output_func = (sol, i) -> (sol, false),
    reduction = (u, data, I) -> (append!(u, data), false)
)


#===============================================================================#

#PROBLEM SETUP
#
"""
"Sketch" the shape of the ODE-problem. This will be remade by the problem_func
for each particle in the EnsembleProblem, but the particle independent `eom`,
`tspan` and `charge`, `mass`, and `fields` may be used.
"""
prob = ODEProblem(
    eom,
    u0,
    tspan,
    (
        charge = charge,
        mass = mass,
        fields = fields_itp,
        mu0,
    )
)
ensemble_prob = EnsembleProblem(
        prob
        ;
        ensembleprob_kwargs_internal...
)

#===============================================================================#
# START OF SIMULATION

println("Running ensemble of $npart particles...")
println(".................................................................",
        "...............")

@time sol = DifferentialEquations.solve(
    ensemble_prob,
    solve_args...
    ;
    solve_kwargs_internal...
);


println("\nPost-processing: Computing GCAStates...")
gcastates = Vector{Vector{GCAState}}(undef, npart)
ntimes = try ntimes
catch
    1000
end
for i in 1:npart
    times = range(sol[i].t[1], sol[i].t[end], length=ntimes)
    gcastates[i] = tp.timeseriesofgcastates(sol[i], times)
end

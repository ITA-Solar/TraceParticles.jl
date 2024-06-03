using DifferentialEquations

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
    zeros(1), # This is a dummy initial condition
    tspan,
    (
        charge = charge,
        mass = mass, 
        fields = fields_itp,
    )
)
ensemble_prob = EnsembleProblem(
        prob
        ;
        ensembleprob_kwargs...
    )

#===============================================================================#
# START OF SIMULATION

println("Running ensemble of $npart particles...")
@time sim = DifferentialEquations.solve(
    ensemble_prob,
    EnsembleSerial()
    ;
    solve_kwargs...
    );

if !(:reduction in keys(ensembleprob_kwargs))
    @warn "No reduction function specified. 
    Saving the full ensemble simulation as 'out.jld2'."
    JLD2.@save "out.jld2" sim
end

using tp
using DifferentialEquations

"""
We will in this example simulate two charged particles embedded in an
electromagnetic field.
"""
eom = lorentzforce! # The equations of motion of charge particles
                    # in a collisionless electromagntic field is the
                    # Lorentz force.
# Here we define the initial conditions of our two particles. The initial
# conditions are 3 position- and 3 velocity-components, which could be
# represented in a 6-component state vector.
ic = [[0.0, 1.0, 0.0, 1.0, 0.0, 0.0], [0.0, 0.5, 0.0, -0.5, 0.0, 0.0]]
# The parameters of the Lorentz equation is the particle's charge, mass and
# functor a giving the electromagnetic fields at the particle position.
mass = 1
charge = [1, -1]
struct StaticUniformElectromagneticField
    bx::Real
    by::Real
    bz::Real
    ex::Real
    ey::Real
    ez::Real
end
function (field::StaticUniformElectromagneticField)(x,y,z)
    return [field.bx, field.by, field.bz, field.ex, field.ey, field.ez]
end
emfield = StaticUniformElectromagneticField(0,0,1,0,1,0)
tspan = (0,10) # Define the simulation time-span
npart = 2     # Define the number of particles in the simulation
prob = ODEProblem(eom, ic[1], tspan) # Create an initial ODE-problem

# Define a problem function, which will define a new ODE-problem for each
# particle.
function prob_func(prob, i, repeat)
    remake(prob, u0=ic[i], p=(mass, charge[i], emfield))
end
# Define an output function - what to be saved from each particle solution
# The ODESolution type contains a lot of metadata, including alogorithm
# choice, number of algorithm switches, etc. In a large ensemble, this
# overhead will quickly consume a lot of memory.
#
# In this output function we save only the initial and final state of the
# particle, in addition to its timesteps. We set the required "rerun" return
# argument to false.
function output_func(sol,i)
    return ([first(sol), last(sol)], sol.t), false
end
# Define an Ensamble-problem, to simulate the ensemble of particles
ensamble_prob = EnsembleProblem(
    prob # Initial problem
    ;
    prob_func=prob_func,
    output_func=output_func
)

# Run the simulation
@time sim = DifferentialEquations.solve(
    ensamble_prob, 
    EnsembleSerial()
    ;
    trajectories=npart
)

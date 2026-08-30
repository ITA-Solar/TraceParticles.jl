using Interpolations
using LinearAlgebra
using OrdinaryDiffEq
using StaticArrays
using Test
using TraceParticles
import TraceParticles as TP


function magneticbottle(x, y, z, _)
    B0 = 10.0
    L = 0.4    # Mirror length
    a = B0 * z / L^2
    # Zero electric field
    return zeros(3), [-x * a, -y * a, B0 + z * a]
end 

dt = 1.e-3
tf = 300.0
tspan = (0, tf)

mass = 1
charge = 1

# INITIAL CONDITIONS
vel0 = [0.0, 0.1, 0.1]
rL = mass * √(vel0[1]^2 + vel0[2]^2) / (charge * B0)
pos0 = [-1rL, 0.0, 0.0]

#-------------------------------------------------------------------------------
# SIMULATION DURATION
#
numsteps = trunc(Int64, tf / dt)   # Number of timesteps in the simulation
#println("Number of time steps = $numsteps.")

#-------------------------------------------------------------------------------
# PARTICLE CREATION
# Set initial position and velocities
E, B = emfields_itp(pos0...)
R0, vparal, μ = TP.get_guidingcentre(
    pos0,
    vel0,
    B,
    E,
    charge,
    mass
)
# To make R0 a vector instead of SVector which crashes the ODE solver with an
# in-place EoM.
R0 = [R0...]
#-------------------------------------------------------------------------------
# CREATE PROBLEM
prob_FO = ODEProblem(
    lorentzforce!,
    [pos0; vel0],
    tspan,
    (charge=charge, mass=mass, electromagneticfield=emfields_itp),
)
prob_GCA = ODEProblem(
    guidingcentreapproximation!,
    [R0; vparal],
    tspan,
    (
        charge=charge,
        mass=mass,
        magneticmoment=μ,
        electromagneticfield=emfields_itp
    )
)

#...............................................................................
# RUN SIMULATION
sol_FO = solve(prob_FO, Tsit5();
    reltol=5e-7, abstol=1e-9)
sol_GCA = solve(prob_GCA, Tsit5();
    reltol=1e-4, abstol=1e-9)

#...............................................................................


#-------------------------------------------------------------------------------
# TESTING
Bmax = B0 * norm(vel0)^2 / (vel0[1]^2 + vel0[2]^2)
zmax = L * √(Bmax / B0 - 1)

function z(t, μ, mass, B0, L, A, ϕ)
    ω = √(2μ * B0 / (mass * L^2))
    return @. A * sin(ω * t + ϕ)
end
ϕ = 0
A = zmax
times = collect(range(0.0, step=dt, length=numsteps + 1))
z_anal_FO = z(sol_FO.t, μ, mass, B0, L, A, ϕ)
z_anal_GCA = z(sol_GCA.t, μ, mass, B0, L, A, ϕ)

times = range(0.0, tf, length=1000)
analytical = z.(times, μ, mass, B0, L, A, ϕ)

z_GCA = [u[3] for u in sol_GCA.u]
z_FO = [u[3] for u in sol_FO.u]
rmse_GCA = √(sum((z_anal_GCA .- z_GCA) .^ 2) / numsteps)
rmse_FO = √(sum((z_anal_FO .- z_FO) .^ 2) / numsteps)
@testset verbose = true "GCA: Default alg." begin
    @test isapprox(rmse_GCA, 0.0, atol=1e-5)
end # testset GCA: Euler
@testset verbose = true "Full orbit: Default alg." begin
    @test isapprox(rmse_FO, 0.0, atol=3e-4)
e

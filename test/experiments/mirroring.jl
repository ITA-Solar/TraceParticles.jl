# Created 02.04.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
"""
    mirroring.jl

This experiments tests the mirroring term in the GCA solver and mirroring in
the full orbit solver by comparing with an analytic expression of the motion
parallel to the magnetic field. The analytic magnetic bottle field is obtained
from Ripperda et al. 2018, and the motion of the particle is assume adiabatic
such that the GCA holds.

The GC is initialised in the centre of the bottle to avoid any other drifts.
"""

#-------------------------------------------------------------------------------
#                            EXPERIMENTAL PARAMETERS
#
#-------------------------------------------------------------------------------
# NUMBER OF PARTICLES, SIMULATION DURATION, TIMESTEP   |
#......................................................|
numparticles = 1  # Number of particles to simulate    |
dt = 1.e-3        # Time step [s]                      |
tf = 300.0 #n=100   # End time of simulation [s]         |
tspan = (0, tf)
#tf = 2.3 #n=10   # End time of simulation [s]         |
#......................................................|

#...............................................
# PARTICLE TYPE
mass = 1
charge = 1

#...............................................
# MAGNETIC FIELD PARAMETERS
B0 = 10.0 # Magnetic field strenth parameter
L = 0.4 # Mirroring length

#...............................................
# INITIAL CONDITIONS
vel0 = [0.0, 0.1, 0.1]
rL = mass*√(vel0[1]^2 + vel0[2]^2)/(charge*B0)
pos0 = [-rL, 0.0, 0.0]

#...............................................
# SPATIAL PARAMETERS (x, y, z)
numdims = 3
# Lower bounds of the three spatial axes
a = 2.1L
xi0 = (-a, -a, -a)
# Upper bound of the three spatial axes
xif = (a, a, a)
# Grid resolution of the axes
#n = (10, 10, 2)
ni = (100, 100, 100)
 
#...............................................
# SOLVER CONDITIONS
# RK4 for GCA? Euler-Cromer for full orbit?

#-------------------------------------------------------------------------------
#                            RUN EXPERIMENT
#
#-------------------------------------------------------------------------------
# COMPUTING THE AXES, MAGNETIC FIELD AND ELECTRIC FIELD
analytical_field = false
if analytical_field
    emfields_itp(x,y,z) = [magneticmirrorfield(x,y,z,B0,L); zeros(3)]
else
    xx, yy, zz, dx, dy, dz = create3Daxes(xi0, xif, ni)
    Bfield = zeros(Float64, numdims, ni[1], ni[2], ni[3])
    Efield = zeros(size(Bfield))
    Bfield = stack(discretise(
        (x, y, z) -> magneticmirrorfield(x,y,z,B0,L),
        xx,
        yy,
        zz
        )
    )
    #emfields = eachslice(vcat(Bfield, Efield), dims=(2,3,4))
    fields = [Bfield; Efield]
    emfields_itp = Vector{AbstractInterpolation}(undef, 6)
    for i = 1:6
        emfields_itp[i] = cubic_spline_interpolation((xx, yy, zz), fields[i,:,:,:],
            extrapolation_bc=Flat()
            )
    end
end
#-------------------------------------------------------------------------------
# SIMULATION DURATION
#
numsteps = trunc(Int64, tf/dt)   # Number of timesteps in the simulation
#println("Number of time steps = $numsteps.")

#-------------------------------------------------------------------------------
# PARTICLE CREATION
# Set initial position and velocities
B⃗ = [itp(pos0...) for itp in emfields_itp][1:3]
E⃗ = zeros(3)
R0, vparal, μ = tp.get_guidingcentre(pos0, vel0, B⃗, E⃗, charge, mass)
#-------------------------------------------------------------------------------
# CREATE PROBLEM
prob_FO = ODEProblem(
    lorentzforce!, 
    [pos0; vel0],
    tspan,
    (charge=charge, mass=mass, fields=emfields_itp),
    )
prob_GCA = ODEProblem(
    guidingcentreapproximation!,
    [R0; vparal],
    tspan,
    (charge=charge, mass=mass, magneticmoment=μ, fields=emfields_itp)
    )

#...............................................................................
# RUN SIMULATION
sol_FO = DifferentialEquations.solve(prob_FO, Tsit5();
    reltol=5e-5)
sol_GCA = DifferentialEquations.solve(prob_GCA, Tsit5();
    reltol=1e-4, abstol=1e-9)

#...............................................................................


#-------------------------------------------------------------------------------
# TESTING
Bmax = B0*norm(vel0)^2/(vel0[1]^2 + vel0[2]^2)
zmax = L*√(Bmax/B0 - 1)

function z(t, μ, mass, B0, L, A, ϕ)
    ω = √(2μ*B0/(mass*L^2))
    return @. A*sin(ω*t + ϕ)
end
ϕ = 0
A = zmax
times = collect(range(0.0, step=dt, length=numsteps+1))
z_anal_FO = z(sol_FO.t, μ, mass, B0, L, A, ϕ)
z_anal_GCA = z(sol_GCA.t, μ, mass, B0, L, A, ϕ)

times = range(0.0, tf, length=1000)
analytical = z.(times, μ, mass, B0, L, A, ϕ)

z_GCA  = [u[3] for u in sol_GCA.u] 
z_FO  = [u[3] for u in sol_FO.u] 
rmse_GCA = √(sum((z_anal_GCA .- z_GCA).^2)/numsteps)
rmse_FO = √(sum((z_anal_FO .- z_FO).^2)/numsteps)
@testset verbose = true "GCA: Default alg." begin
    @test isapprox(rmse_GCA, 0.0, atol=0.001)
end # testset GCA: Euler
@testset verbose = true "Full orbit: Default alg." begin
    @test isapprox(rmse_FO, 0.0, atol=0.001)
end # testset GCA: Euler
    

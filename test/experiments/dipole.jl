# Created 14.04.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 dipole.jl
#
#-------------------------------------------------------------------------------
# Runs a GCA and full orbit particle in a magnetic dipole. The setup is the same
# as in Ripperda et al., 2018. Initial positions is set such that the guiding
# centre of the particle is at x = [1.0, 0.0, 0.0]. See 'dipoleloop.jl' for the
# testsets that compares the result from the solvers with each other and the
# 'T_dipole' approximation. 
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#                            EXPERIMENTAL PARAMETERS
#
#-------------------------------------------------------------------------------
# NUMBER OF PARTICLES, SIMULATION DURATION, TIMESTEP   |
#......................................................|
numparticles = 1  # Number of particles to simulate    |
tf = 80.0 #n=100   # End time of simulation [s]         |
tspan = (0.0, tf)
#tf = 2.3 #n=10   # End time of simulation [s]         |
#......................................................|

#...............................................
# PARTICLE TYPE
mass = 1.0
charge = 1.0

#...............................................
# MAGNETIC FIELD PARAMETERS
#qMm = 20.0
#println("qMm: $qMm")
M = qMm * mass / charge
itp_type = Gridded(Linear())
itp_bc = Flat()

#...............................................
# INITIAL CONDITIONS
# Initial position is such that vperp is 1 and the guiding centre is
# at x⃗ = [1, 0, 0]. It is also assumed that the GCA is valid such that B is equal
# to B(x⃗) = 2M, directed in -ẑ. We may then compare with the Walt (1994)
# approximation of the azimuthal drift-period.
R0 = SA[1.0, 0, 0]
vparal = 0.5
vperp = 1.0

#...............................................
# SPATIAL PARAMETERS (x, y, z)
numdims = 3
# Lower bounds of the three spatial axes
a = 1.5
xi0 = (-a, -a, -a)
# Upper bound of the three spatial axes
xif = (a, a, a)
# Grid resolution of the axes
#n = (10, 10, 2)
ni = (256, 256, 256)

#...............................................
# ELECTRIC FIELD PARAMETERS
# It will be zero

#-------------------------------------------------------------------------------
#                            RUN EXPERIMENT
#
#-------------------------------------------------------------------------------
# COMPUTING THE AXES, MAGNETIC FIELD AND ELECTRIC FIELD
analytical_field = false
if analytical_field
    emfields(x, y, z) = [magneticdipolefield(x, y, z, M); zeros(3)]
    emfields_itp = Vector{Function}(undef, 6)
    for i = 1:6
        emfields_itp[i] = (x, y, z) -> emfields(x, y, z)[i]
    end
else
    xx, yy, zz, dx, dy, dz = create3Daxes(xi0, xif, ni)
    Efield = zeros(Float64, numdims, ni...)
    Bfield = stack(discretise(
        (x, y, z) -> magneticdipolefield(x, y, z, M),
        xx,
        yy,
        zz,
    )
    )
    # Create interpolation objects
    #emfields = eachslice(vcat(Bfield, Efield), dims=(2,3,4))
    fields = [Bfield; Efield]
    emfields_itp = Vector{AbstractInterpolation}(undef, 6)
    for i = 1:6
        emfields_itp[i] = cubic_spline_interpolation((xx, yy, zz), fields[i, :, :, :],
            extrapolation_bc=Flat()
        )
    end
end
emfields_itp = ElectromagneticFieldInterpolator(
    [StaticInterpolation(itp) for itp in emfields_itp]
)

#-------------------------------------------------------------------------------
# PARTICLE CREATION
# Set initial position and velocities
E, B = emfields_itp(R0[1], R0[2], R0[3])
μ = TP.magneticmoment(vperp, mass, norm(B))
u0 = TP.get_fullorbit(B, E, R0, vparal, μ, charge, mass, pi / 2)
u0 = [u0[1], u0[2], u0[3], u0[4], u0[5], u0[6]]  # Make not static

fo_params = (charge=charge, mass=mass, electromagneticfield=emfields_itp)
gca_params = (
    charge=charge,
    mass=mass,
    magneticmoment=μ,
    electromagneticfield=emfields_itp
)

#-------------------------------------------------------------------------------
# CREATE PROBLEM
#fo_prob = ODEProblem(lorentzforce!, [pos0; vel0], tspan, fo_params)
fo_prob = ODEProblem(lorentzforce!, u0, tspan, fo_params)
gca_prob = ODEProblem(
    guidingcentreapproximation!, [R0...; vparal], tspan, gca_params
)

#...............................................................................
# RUN SIMULATION
reltol_fo = 1e-7
abstol_fo = 1e-9
reltol_gca = 1e-4
abstol_gca = 1e-9
fo_sim = solve(
    fo_prob, Tsit5();
    reltol=reltol_fo, abstol=abstol_fo
)
gca_sim = solve(
    gca_prob, Tsit5();
    reltol=reltol_gca, abstol=abstol_gca
)
#...............................................................................


#-------------------------------------------------------------------------------
# TESTING
ϕ_FO = atan(fo_sim.u[end][2] / fo_sim.u[end][1])
ϕ_GCA = atan(gca_sim.u[end][2] / gca_sim.u[end][1])

#v = norm(vel0)
#vparall = abs((B⃗ ⋅ vel0)/B)
#vperp = √(v^2 - vparall^2)

v = norm(u0[4:6])
α = atan(vperp / vparal)
R0x = R0[1]

if !isdefined(Main, :T_dipole)
    function T_dipole(q, m, M, R0, v, α)
        return @. 2π * q * M / (m * v^2 * R0) * (1 - 1 / 3 * sin(α)^0.62)
    end # function T_dipole
end

T = T_dipole(charge, mass, M, R0x, v, α)
angularfreq = 2π / T

times = range(0.0, tf, length=1000)
driftx = R0x * cos.(angularfreq * times)
drifty = R0x * sin.(angularfreq * times)
angle_anal = atan.(drifty, driftx)
angle_fo = atan.(fo_sim(times)[2, :], fo_sim(times)[1, :])
angle_gca = atan.(gca_sim(times)[2, :], gca_sim(times)[1, :])
angle_fo = atan.(fo_sim(times)[2, :], fo_sim(times)[1, :])

global ϕ = 2π * tf / T
global ϕ_FO = angle_fo[end]
global ϕ_gca = angle_gca[end]

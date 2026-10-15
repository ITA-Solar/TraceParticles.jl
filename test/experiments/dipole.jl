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

function dipolesimulation(qMm)
    # EXPERIMENTAL PARAMETERS
    tf = 80.0 # End time of simulation [s]
    tspan = (0.0, tf)
    mass = 1.0
    charge = 1.0

    # MAGNETIC FIELD PARAMETERS
    #qMm = 20.0
    M = qMm * mass / charge
    f = 80.0 # En
    # ELECTRIC FIELD PARAMETERS
    # It will be zero

    # INITIAL CONDITIONS
    # Initial position is such that vperp is 1 and the guiding centre is
    # at x⃗ = [1, 0, 0]. It is also assumed that the GCA is valid such that B is equal
    # to B(x⃗) = 2M, directed in -ẑ. We may then compare with the Walt (1994)
    # approximation of the azimuthal drift-period.
    R0 = SA[1.0, 0, 0]
    vparal = 0.5
    vperp = 1.0

    # SPATIAL PARAMETERS (x, y, z)
    # Lower bounds of the three spatial axes
    a = 1.5
    xi0 = (-a, -a, -a)
    # Upper bound of the three spatial axes
    xif = (a, a, a)
    # Grid resolution of the axes
    ni = (256, 256, 256)

    # CREATING THE ELECTROMAGNETIC FIELD
    xx, yy, zz, dx, dy, dz = create3Daxes(xi0, xif, ni)
    bfield = discretise3D(
        (x, y, z) -> magneticdipolefield(x, y, z, M),
        xx,
        yy,
        zz,
    )
    bfield_itp = cubic_spline_interpolation(
        (xx, yy, zz),
        bfield,
        extrapolation_bc=Flat()
    )
    electromagneticfield(x, y, z, _) = zeros(SVector{3}), bfield_itp(x, y, z)

    # PARTICLE CREATION
    # Set initial position and velocities
    E, B = electromagneticfield(R0[1], R0[2], R0[3], 0)
    μ = magneticmoment(vperp, mass, norm(B))
    u0 = get_fullorbit(B, E, R0, vparal, μ, charge, mass, pi / 2)
    u0 = [u0[1], u0[2], u0[3], u0[4], u0[5], u0[6]]  # Make not static

    fo_params = (charge=charge, mass=mass, electromagneticfield=electromagneticfield)
    gca_params = (
        charge=charge,
        mass=mass,
        magneticmoment=μ,
        electromagneticfield=electromagneticfield
    )


    # CREATE PROBLEM
    #fo_prob = ODEProblem(lorentzforce!, [pos0; vel0], tspan, fo_params)
    fo_prob = ODEProblem(lorentzforce!, u0, tspan, fo_params)
    gca_prob = ODEProblem(
        guidingcentreapproximation!, [R0...; vparal], tspan, gca_params
    )

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

    # AZIMUTHAL ANGLES
    ϕ_FO = atan(fo_sim.u[end][2] / fo_sim.u[end][1])
    ϕ_GCA = atan(gca_sim.u[end][2] / gca_sim.u[end][1])

    # Analytical
    v = norm(u0[4:6])
    α = atan(vperp / vparal)
    R0x = R0[1]
    T = T_dipole(charge, mass, M, R0x, v, α)
    ϕ = 2π * tf / T

    return ϕ, ϕ_FO, ϕ_GCA, fo_sim, gca_sim
end

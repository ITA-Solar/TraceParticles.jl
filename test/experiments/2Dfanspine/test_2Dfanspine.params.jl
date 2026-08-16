using OrdinaryDiffEq
using TraceParticles

# EXPEIRMENT PARAMETERS
#==============================================================================#
params = TraceParticlesParameters(
    emfield_file = Base.source_dir() * "/cs_t1_BE.jld2",
    static_itp_wrapper = false,
    xz_itp_wrapper = false,
    tf = 1e-2 * 100, # dtsnap is hs
    # Particle parameters
    charge = -TraceParticles.e, # Charge of particles
    mass = TraceParticles.m_e, # Mass of particles
    eom = guidingcentreapproximation!, # Equation of motion
    npart = Int(1e2),
    seed = 5,
    # Parameters for saving data
    expname = "testresults",
    batchsize = Int(1e1),
    # Where to save the results
    datadir = Base.source_dir(),
    # Other stuff
    alg = Rosenbrock23(),
    save_max_observables = true,
    reltol = 3e-8,
    abstol = 1e0,
    maxiters = 10_000_000,
    verbose = false,

    # Initial conditions: sampled from the snapshot during the run.
    # The whole snapshot except a 40km edge on both sides.
    prob_func = MHDSample(
        tg_file = Base.source_dir() * "/cs_t1_tg_normfalse.jld2",
        target_distr_file = Base.source_dir() * "/cs_t1_r_normfalse.jld2",
        xbounds = (1.5822270393371582e7, 1.6318359375e7),
        ybounds = (1e6, 1e6),
        zbounds = (-6.598462104797363e6, -6.102306842803955e6),
        t0bounds = (0.0, 0.0),
    ),

    # Termination criteria
    callback_specs = (
        OutOfBounds(
            x = (1.5822270393371582e7, 1.6318359375e7),
            z = (-6.598462104797363e6, -6.102306842803955e6),
            save_positions = (true, true),
        ),
        Relativistic(fraction = 0.02, save_positions = (true, true)),
        HighMagneticGradient(tolerance = 0.001, save_positions = (true, true)),
    ),
)
#==============================================================================#

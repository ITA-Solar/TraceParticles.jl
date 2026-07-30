using OrdinaryDiffEq
using TraceParticles

# EXPEIRMENT PARAMETERS
#==============================================================================#
params = TraceParticlesParameters(
    precision = Float64,
    emfield_file = Base.source_dir() * "/linear_t1_BE.jld2",
    tg_file = Base.source_dir() * "/linear_t1_tg_normfalse.jld2",
    target_distr_file = Base.source_dir() * "/linear_t1_r_normfalse.jld2",
    static_itp_wrapper = false,
    xz_itp_wrapper = false,
    ic_onthefly = true,

    t0start = 0.0,
    t0end = 0.0,
    tf = 1e-2 * 100, # dtsnap is hs
    #xbounds = (14.527344e6, 18.824219e6) # Bounds of of initial position in x
    # x bounds for the whole simulation snapshot
    # The whole snapshot except a 40km edge on both sides.
    xstart_ic = 1.5822270393371582e7,
    xend_ic = 1.6318359375e7,
    ystart_ic = 1e6,
    yend_ic = 1e6,
    zstart_ic = -6.598462104797363e6,
    zend_ic = -6.102306842803955e6,

    # Particle parameters
    charge = -TraceParticles.e, # Charge of particles
    mass = TraceParticles.m_e, # Mass of particles
    eom = guidingcentreapproximation!, # Equation of motion
    npart = Int(1e2),
    seed = 10,
    # Parameters for saving data
    expname = "2Dfanspine_run1",
    batchsize = Int(1e6),
    # Where to save the results
    datadir = Base.source_dir(),
    # Define parameters for callbacks
    kill_oob = true,
    xkill = true,
    zkill = true,
    xlowerbound = 1.5822270393371582e7,
    xupperbound = 1.6318359375e7,
    zlowerbound = -6.598462104797363e6,
    zupperbound = -6.102306842803955e6,
    oob_save_positions = (true, true),
    kill_relativistic = true,
    relativistic_fraction = 0.02,
    rel_save_positions = (true, true),
    kill_high_gradb = true,
    gradient_tolerance = 0.001,
    gradb_save_positions = (true, true),
    # Other stuff
    alg = Rosenbrock23(),
    save_max_observables = false,
    reltol = 3e-8,
    abstol = 1e0,
    maxiters = 10_000_000,
    verbose = false,
)
#==============================================================================#

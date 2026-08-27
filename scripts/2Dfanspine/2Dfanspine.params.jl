using OrdinaryDiffEq
using TraceParticles
import TraceParticles as TP

# EXPERIMENT PARAMETERS
#==============================================================================#
params = TraceParticlesParameters(
    emfield_file = Base.source_dir() * "/linear_t1_BE.jld2",
    static_itp_wrapper = false,
    xz_itp_wrapper = false,
    tf = 1e-2 * 100, # dtsnap is hs
    # Particle parameters
    charge = TP.electron_charge,
    mass = TP.electron_mass,
    eom = hybridgcafo!, # Equation of motion
    npart = Int(1e2),
    seed = 10,
    initialeomid = EoMID.FullOrbit,
    # Parameters for saving data
    expname = "test",
    batchsize = Int(1e6),
    # Where to save the results
    datadir = Base.source_dir(),
     # Solve parameters
    alg = Tsit5(),
    reltol = 3e-8,
    abstol = 1e0,
    maxiters = 10_000,

    # Initial conditions: sampled from the snapshot during the run.
    #xbounds = (14.527344e6, 18.824219e6) # Bounds of of initial position in x
    # x bounds for the whole simulation snapshot
    # The whole snapshot except a 40km edge on both sides.
    prob_func = SampleICsFromMHD(
        tg_file = Base.source_dir() * "/linear_t1_tg_normfalse.jld2",
        target_distr_file = Base.source_dir() * "/linear_t1_r_normfalse.jld2",
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
        ),
        Relativistic(fraction = 0.06)
    ),
)
#==============================================================================#

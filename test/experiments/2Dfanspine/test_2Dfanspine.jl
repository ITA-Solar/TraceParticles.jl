#-------------------------------------------------------------------------------
# Created 06.08.26
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
"""
    test_2Dfanspine.jl

Runs the complete `TraceParticles` pipeline on a small 2D fan-spine Bifrost
snapshot and tests properties that must hold on *any* machine. The exact
trajectories of the particles in the snapshot are chaotic, so we do not compare
the solutions to exact or approxiamate values.

- the Bifrost snapshot is turned into usable field interpolators,
- the sampled initial conditions obey the sampling parameters,
- the ensemble output is laid out and read back consistently,
- every particle terminates for a reason that is verifiable from its final
  state (out of bounds, magnetic-gradient, relativistic, or not at all),
- the energies computed by `save_energy` and `save_gcastates` agree, and the
  guiding-centre states satisfy the defining relations of the GCA, and
- the simulation is reproducible for a given seed.

All output is written to a temporary directory, which is removed when the
process exits.
"""

using BifrostTools
using DataFrames
using Interpolations
using JLD2
using LinearAlgebra
using Logging
using OrdinaryDiffEq
using Random
using SciMLBase
using StaticArrays
using Test
using TraceParticles
import TraceParticles as TP

expdir = Base.source_dir()
outdir = mktempdir()

#-------------------------------------------------------------------------------
# Create interpolation objects from the Bifrost snapshot.
brxp = BifrostExperiment("testsnap", expdir)
create_bifrost_itps(
    brxp,
    1,
    BSpline(Cubic()),
    Flat(),
    ["BE", "tg", "r"];
    name=joinpath(outdir, "test"),
    units="si",
    normalise=[false, false, false],
    destagger=[true, false, false],
    xzinterpolation=true,
)

emfield_file = joinpath(outdir, "test_t1_BE_si_cubicspline.itp.jld2")
tg_file = joinpath(outdir, "test_t1_tg_si_cubicspline.itp.jld2")
target_distr_file = joinpath(outdir, "test_t1_r_si_cubicspline.itp.jld2")

# Bounds of the snapshot, except a 40 km edge on both sides.
xbounds = (1.5822270393371582e7, 1.6318359375e7)
zbounds = (-6.598462104797363e6, -6.102306842803955e6)
ypos = 1e6
ybounds = (ypos, ypos)

# All particles start at `t0` and are integrated until `tf`. 
t0 = 0.0
tf = 1e-2 * 100 # dtsnap is hs
# Termination callback parameters
relativisticfraction = 0.02
gradienttolerance = 0.001

"""
    fanspineparameters(expname; seed)
Parameters of the fan-spine experiment, writing its results to the group
`expname` in the temporary output directory.
"""
fanspineparameters(expname; seed=5) = TraceParticlesParameters(
    emfield_file=emfield_file,
    # The JLD2 files already contain StaticXZInterpolation-objects
    static_itp_wrapper=false,
    xz_itp_wrapper=false,
    tf=tf,

    # Particle parameters
    charge=TP.electron_charge,
    mass=TP.electron_mass,
    eom=guidingcentreapproximation!,
    npart=Int(1e2),
    seed=seed,
    # Parameters for saving data
    expname=expname,
    batchsize=Int(1e1),
    datadir=outdir,
    # Other stuff
    alg=Rosenbrock23(),
    save_max_observables=true,
    reltol=3e-8,
    abstol=1e0,
    maxiters=10_000_000,
    verbose=false,
    log2file=false,
    solver_options = (;
        save_idxs = [1,2,3,4], # Just for testing, this is the default
        save_everystep=true    # Just for testing, this is the default
    ),

    # Initial conditions: sampled from the snapshot during the run.
    prob_func=SampleICsFromMHD(
        tg_file=tg_file,
        target_distr_file=target_distr_file,
        xbounds=xbounds,
        ybounds=ybounds,
        zbounds=zbounds,
        t0bounds=(t0, t0),
    ),

    # Termination criteria
    callback_specs=(
        OutOfBounds(x=xbounds, z=zbounds, save_positions=(true, true)),
        Relativistic(
            fraction=relativisticfraction, save_positions=(true, true)
        ),
        HighMagneticGradient(
            tolerance=gradienttolerance, save_positions=(true, true)
        ),
    ),
)

"""
    runexperiment(params)
Construct and solve the fan-spine problem and return the problem together with
the filename of its results.
"""
function runexperiment(params::TraceParticlesParameters)
    tpprob = TraceParticlesProblem(params; logger=current_logger())
    solve(tpprob, params)
    fname = joinpath(outdir, params.expname, params.expname * ".h5")
    return tpprob, fname
end

params = fanspineparameters("testresults")
tpprob, fname = runexperiment(params)
electromagneticfield = tpprob.ensemble_prob.prob.p.electromagneticfield
TraceParticles.save_gcastates(fname, electromagneticfield)

# The results of the ensemble, and the initial and final energies in eV.
df = h5_getall(fname)
e0, ef = h5_getenergies(fname)

# Energy at which particles are killed by the relativistic callback
relativisticenergy = relativisticfraction * params.mass *
                     TraceParticles.csqrd * TraceParticles.J2eV
# Termination codes as they are stored in the results
notterminated = Int32(TerminationCode.NotTerminated)
outofbounds = Int32(TerminationCode.OutOfBounds)
magneticgradient = Int32(TerminationCode.MagneticGradient)
relativistic = Int32(TerminationCode.Relativistic)

"""
    isinsidedomain(x, z)
Whether the position (`x`, `z`) is inside the domain in which the particles are
traced.
"""
isinsidedomain(x, z) =
    xbounds[1] <= x <= xbounds[2] && zbounds[1] <= z <= zbounds[2]

#-------------------------------------------------------------------------------
#                                   TESTS
#-------------------------------------------------------------------------------
@testset "Field interpolators" begin
    @test isfile(emfield_file)
    @test isfile(tg_file)
    @test isfile(target_distr_file)

    # A point in the middle of the snapshot
    x, z = sum(xbounds) / 2, sum(zbounds) / 2
    E, B = electromagneticfield(x, ypos, z, t0)
    @test length(E) == length(B) == 3
    @test all(isfinite, E)
    @test all(isfinite, B)
    @test norm(B) > 0
    @test norm(E) > 0
    # There should be no electric field parallel to the plane and no magnetic
    # field out of the plane
    @test all(iszero, [E[1:2:3]; B[2]])
    # The snapshot is a 2D (x, z) slice, so the field must be independent of y
    @test electromagneticfield(x, 2ypos, z, t0) == (E, B)

    # Temperature and mass density are positive definite
    @test load_object(tg_file)(x, ypos, z, t0) > 0
    @test load_object(target_distr_file)(x, ypos, z, t0) > 0
end

@testset "Initial condition sampling" begin
    # `SampleICsFromMHD` resolves into a `SampleGCA` problem function, which
    # only uses the random number generator of the ensemble context, so it can
    # be called directly with a given generator.
    prob = tpprob.ensemble_prob.prob
    @test tpprob.ensemble_prob.prob_func isa SampleGCA
    sample(seed) = tpprob.ensemble_prob.prob_func(prob, (sim_id=1, rng=Xoshiro(seed)))

    particle = sample(42)
    @test sample(42).u0 == particle.u0     # The sampling is reproducible
    @test sample(43).u0 != particle.u0     # and does depend on the generator

    # The guiding centre state is (x, y, z, vparal)
    @test length(particle.u0) == 4
    @test particle.tspan == (t0, tf)
    @test particle.p isa GCAParams
    @test particle.p.magneticmoment > 0
    @test particle.p.nrejections >= 0
    # The proposal distribution is the target distribution, so all particles
    # carry unit statistical weight.
    @test particle.p.weight == 1

    # The guiding centre is sampled inside the domain, up to the Larmor radius
    # separating it from the particle itself.
    x, y, z, vparal = particle.u0
    r_L = TraceParticles.larmorradius(
        SVector(x, y, z), t0, particle.p.magneticmoment,
        params.charge, params.mass, electromagneticfield
    )
    # The particles are sampled in the plane y = `ypos`, so the offset of the
    # guiding centre from that plane is the y-component of the Larmor radius
    # vector. It is evaluated at the guiding centre rather than at the particle
    # itself, hence the 10 % margin.
    @test xbounds[1] - r_L <= x <= xbounds[2] + r_L
    @test zbounds[1] - r_L <= z <= zbounds[2] + r_L
    @test abs(y - ypos) <= 1.1r_L
    @test isfinite(vparal)
end

@testset "Ensemble results layout" begin
    batches = h5_getbatchnames(fname)
    @test length(batches) == params.npart ÷ params.batchsize
    @test nrow(df) == params.npart
    @test issubset(
        ["x0", "y0", "z0", "vparal0", "xf", "yf", "zf", "vparalf", "t0", "tf",
            "charge", "mass", "magneticmoment", "weight", "nrejections", "nt",
            "retcode", "terminationcode", "maxrl", "minlb", "maxscalesratio"],
        names(df)
    )
    # `h5_getall` reads the batches in the order returned by
    # `h5_getbatchnames`, so the first batch comes first.
    @test h5_getall(fname, 1) == df[1:params.batchsize, :]
    @test h5_getdataset(fname, "x0") == df.x0

    x0, y0, z0, vparal0, mu0 = h5_getinitialstate_gca(fname)
    @test (x0, y0, z0, vparal0, mu0) ==
          (df.x0, df.y0, df.z0, df.vparal0, df.magneticmoment)
    xf, yf, zf, vparalf, muf = h5_getfinalstate_gca(fname)
    @test (xf, yf, zf, vparalf, muf) ==
          (df.xf, df.yf, df.zf, df.vparalf, df.magneticmoment)

    # Particle properties are those requested in the parameters
    @test all(==(params.charge), df.charge)
    @test all(==(params.mass), df.mass)
    @test all(==(1), df.weight)
    @test all(>=(0), df.nrejections)
    @test all(>=(2), df.nt)
end

@testset "Sampled initial conditions" begin
    @test all(==(t0), df.t0)
    @test all(>(0), df.magneticmoment)
    # The guiding centres are sampled inside the domain, up to the Larmor
    # radius separating them from the particles themselves. `maxrl` is the
    # largest Larmor radius along the trajectory, and hence an upper bound of
    # the initial one.
    @test all(@. xbounds[1] - df.maxrl <= df.x0 <= xbounds[2] + df.maxrl)
    @test all(@. zbounds[1] - df.maxrl <= df.z0 <= zbounds[2] + df.maxrl)
    @test all(@. abs(df.y0 - ypos) <= 1.1df.maxrl)
    @test all(isfinite, df.vparal0)
end

@testset "Termination of particles" begin
    # Particles either run to completion or are stopped by a callback
    @test all(in((notterminated, outofbounds, magneticgradient, relativistic)),
        df.terminationcode)
    @test (df.retcode .== Int32(SciMLBase.ReturnCode.Success)) ==
          (df.terminationcode .== notterminated)
    @test all(in((Int32(SciMLBase.ReturnCode.Success),
            Int32(SciMLBase.ReturnCode.Terminated))), df.retcode)
    @test all(@. t0 < df.tf <= tf)

    @testset "Particles that were not terminated" begin
        m = df.terminationcode .== notterminated
        @test any(m)
        @test all(df.tf[m] .== tf)
        @test all(isinsidedomain.(df.xf[m], df.zf[m]))
    end

    @testset "Particles terminated out of bounds" begin
        m = df.terminationcode .== outofbounds
        @test any(m)
        @test all(.!isinsidedomain.(df.xf[m], df.zf[m]))
    end

    @testset "Particles terminated by the magnetic gradient" begin
        m = df.terminationcode .== magneticgradient
        # The ratio between the Larmor radius and the length scale of the
        # magnetic field exceeds the tolerance at the final position, which is
        # exactly the condition the callback triggers on.
        ratio = [
            TraceParticles.scalesratio(
                SVector(df.xf[i], df.yf[i], df.zf[i]), df.tf[i], df.mass[i],
                df.charge[i], df.magneticmoment[i], electromagneticfield
            ) for i in findall(m)
        ]
        @test all(>=(gradienttolerance), ratio)
        # The largest ratio along the trajectory is the one that triggered the
        # callback, up to the error of the interpolated solution.
        @test all(df.maxscalesratio[m] .>= ratio .* (1 - 1e-6))
    end

    @testset "Particles terminated by the relativistic condition" begin
        m = df.terminationcode .== relativistic
        # Particles are killed above the energy limit, so those that survived
        # the callback must be below it.
        @test all(<(relativisticenergy), ef[.!m])
        @test all(>=(relativisticenergy), ef[m])
    end
end

@testset "Energies" begin
    # `save_energy` and `save_gcastates` must agree on the energies
    @test (e0, ef) == h5_getenergies_gca(fname)
    @test all(isfinite, e0)
    @test all(isfinite, ef)
    @test all(>(0), e0)
    @test all(>(0), ef)
    # The particles are sampled from a Maxwellian in a corona of a few MK,
    # which is nowhere near relativistic.
    @test all(<(relativisticenergy), e0)
end

@testset "Guiding centre states" begin
    gcastates = h5_getall(joinpath(
        outdir, params.expname, params.expname * "_gcastates0.h5"
    ))
    B = @. sqrt(gcastates.bx^2 + gcastates.by^2 + gcastates.bz^2)
    @test all(>(0), B)
    # The defining relations of the guiding centre approximation
    @test gcastates.vperp ≈
          TraceParticles.perpendicular_velocity.(df.magneticmoment, df.mass, B)
    @test gcastates.r_L ≈
          TraceParticles.larmorradius.(df.mass, gcastates.vperp, df.charge, B)
    @test gcastates.ω_c ≈ TraceParticles.gyrofrequency.(df.mass, df.charge, B)
    @test gcastates.L_B ≈ TraceParticles.characteristicfieldlength.(
        [SVector(df.x0[i], df.y0[i], df.z0[i]) for i in 1:nrow(df)],
        df.t0,
        Ref(electromagneticfield)
    )
    # The cosine of the pitch angle
    @test all(@. -1 <= gcastates.β <= 1)
    # The energy of the guiding centre states is what `h5_getenergies` reports
    @test gcastates.energy .* TraceParticles.J2eV ≈ e0
    # The drifts are well defined everywhere
    for drift in ("exbdrift", "∇Bdrift", "curvaturedrift", "polarisationdrift")
        for component in ("x", "y", "z")
            @test all(isfinite, gcastates[!, drift*"_"*component])
        end
    end
    # The E cross B-drift is perpendicular to both fields
    efield = SVector.(gcastates.ex, gcastates.ey, gcastates.ez)
    bfield = SVector.(gcastates.bx, gcastates.by, gcastates.bz)
    drift = SVector.(
        gcastates.exbdrift_x, gcastates.exbdrift_y, gcastates.exbdrift_z
    )
    @test all(@. abs(drift ⋅ bfield) <= 1e-10 * norm(drift) * norm(bfield))
    @test all(@. abs(drift ⋅ efield) <= 1e-10 * norm(drift) * norm(efield))
end

@testset "Reproducibility" begin
    _, repeatedfname = runexperiment(fanspineparameters("repeated"; seed=params.seed))
    @test h5_getall(repeatedfname) == df
end

@testset "Predefined initial conditions" begin
    # Two particles with a magnetic moment large enough to make them
    # relativistic, and one thermal particle, all starting at the same point.
    u0 = [
        [1.63e7, ypos, -6.3e6, 1e6],
        [1.63e7, ypos, -6.3e6, -1e6],
        [1.63e7, ypos, -6.3e6, 1e4],
    ]
    mu = [4.1e-11, 4.1e-11, 4.1e-15]
    tspan = [(t0, tf) for _ in u0]
    odeparams = [
        GCAParams(
            charge=params.charge,
            mass=params.mass,
            magneticmoment=mu[i],
            electromagneticfield=electromagneticfield
        ) for i in eachindex(u0)
    ]
    prob = EnsembleProblem(
        tpprob.ensemble_prob.prob;
        prob_func=PredefinedICs(u0, tspan, odeparams)
    )
    sol = solve(
        prob, params.alg, EnsembleSerial();
        trajectories=length(u0),
        batch_size=length(u0),
        abstol=params.abstol,
        reltol=params.reltol,
        callback=tpprob.callbackset,
        maxiters=params.maxiters,
    )

    energy(s, t) = get_observable(s, :energy, t; EoM="GCA")

    @testset "Predefined states are used as given" begin
        for i in eachindex(u0)
            @test first(sol.u[i].u) == u0[i]
            @test first(sol.u[i].t) == t0
            @test sol.u[i].prob.p.magneticmoment == mu[i]
        end
    end

    @testset "Relativistic particles" begin
        # Both particles start out above the energy limit and must be killed
        # by the relativistic callback at the first step.
        for i in 1:2
            @test sol.u[i].prob.p.terminationcode == TerminationCode.Relativistic
            @test sol.u[i].retcode == SciMLBase.ReturnCode.Terminated
            @test last(sol.u[i].t) < tf
            @test energy(sol.u[i], first(sol.u[i].t)) > relativisticenergy
            @test energy(sol.u[i], last(sol.u[i].t)) > relativisticenergy
        end
        # The two particles are identical except for the sign of their parallel
        # velocity, so they mirror each other about their starting point.
        northbound = last(sol.u[1].u)
        southbound = last(sol.u[2].u)
        @test sign(northbound[3] - u0[1][3]) == -sign(southbound[3] - u0[2][3])
    end

    @testset "Thermal particle" begin
        # A thermal particle is not relativistic, and either survives to the
        # end of the simulation or leaves the domain.
        s = sol.u[3]
        @test s.prob.p.terminationcode != TerminationCode.Relativistic
        @test energy(s, last(s.t)) < relativisticenergy
        if s.prob.p.terminationcode == TerminationCode.NotTerminated
            @test last(s.t) == tf
            @test isinsidedomain(last(s.u)[1], last(s.u)[3])
        end
    end
end

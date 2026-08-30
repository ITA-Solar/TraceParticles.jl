using TraceParticles
using OrdinaryDiffEq: ODEProblem, remake
using SciMLBase: EnsembleProblem, solve
using Random: Xoshiro

outdir = mktempdir()
expname = "demo"
mkpath(joinpath(outdir, expname))

mass = electron_mass      # kg
charge = electron_charge  # C
B = 1e-3 # T
E = 1e2  # V/m
T_c = 2π * mass / abs(charge * B) # s
temperature = 1e6 # K

uniformfield(x, y, z, t) = [0.0, E, 0.0], [0.0, 0.0, B]

x0 = 0.5
y0 = 0.5
z0 = 0.5
t0 = 0.0
rng = Xoshiro(0)
vx0 = TraceParticles.maxwellianvelocitysample(rng, Float64, temperature, mass)
vy0 = TraceParticles.maxwellianvelocitysample(rng, Float64, temperature, mass)
vz0 = TraceParticles.maxwellianvelocitysample(rng, Float64, temperature, mass)
E0, B0 = uniformfield(x0, y0, z0, t0)
R0, vparal0, μ = get_guidingcentre(
    [x0, y0, z0],
    [vx0, vy0, vz0],
    B0, E0, charge, mass
)
u0 = [R0..., vparal0]

tspan = (0.0, 20_000T_c)
npart = 5_000

reduction = SaveBatchAsHDF5(
    expdir=outdir,
    expname=expname,
    batchsize = npart,
    nbatches=1,
    verbose=true,
)

baseprob = ODEProblem(
    guidingcentreapproximation!,
    u0,
    tspan,
    GCAParams(
        charge=charge,
        mass=mass,
        electromagneticfield=uniformfield,
        magneticmoment=μ
    ),
)
eprob = EnsembleProblem(baseprob;
    output_func=output_func_lightweight_gca,
    prob_func = (prob, ctx) -> remake(prob;
        u0 = [R0..., vparal0*randn(rng)],
        p = GCAParams(
            charge=prob.p.charge,
            mass=prob.p.mass,
            electromagneticfield=prob.p.electromagneticfield,
            magneticmoment=prob.p.magneticmoment*rand(rng)
        ),
    ),
    reduction=reduction
)
sim = solve(eprob;
    trajectories=npart,
)

filename = get_filename(reduction)
h5_getdataset(filename, [:xf, :yf, :zf, :vparalf, :magneticmoment])[1:10, :]

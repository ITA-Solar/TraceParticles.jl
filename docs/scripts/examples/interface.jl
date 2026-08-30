using OrdinaryDiffEq: Tsit5, Rosenbrock23
using TraceParticles
using JLD2 # hide


B = 1e-3 # hide
E = 1e2 # hide
gyroperiod = 2π * electron_mass / abs(electron_charge * B) # hide
fielddir = mktempdir() #hide
outdir = mktempdir() #hide
outdir = mktempdir() #hide
uniformfield(x,y,z,t) = [0.0, E, 0.0], [0.0, 0.0, B] #hide
temperature(x,y,z,t) = 1e6 # hide
density(x,y,z,t) = 1.0 # hide
JLD2.save(joinpath(fielddir, "emfield.itp.jld2"), "itp", uniformfield) #hide
JLD2.save(joinpath(fielddir, "temperature.itp.jld2"), "itp", temperature) #hide
JLD2.save(joinpath(fielddir, "particledensity.itp.jld2"), "itp", density) #hide

params = TraceParticlesParameters(
    emfield_file = joinpath(fielddir, "emfield.itp.jld2"),

    tf = 1.0,
    charge = electron_charge,
    mass = electron_mass,
    eom = hybridgcafo!,
    npart = 10,
    seed = 42,
    initialeomid = EoMID.FullOrbit,

    # Parameters for saving data
    expname = "demo",
    # Where to save the results
    datadir = outdir,
     # Solve parameters
    alg = Tsit5(),
    reltol = 1e-6,
    abstol = 1e-6,
    maxiters = 10_000,

    prob_func = SampleICsFromMHD(
        tg_file = joinpath(fielddir, "temperature.itp.jld2"),
        target_distr_file = joinpath(fielddir, "particledensity.itp.jld2"),
        xbounds = (0,1),
        ybounds = (0,1),
        zbounds = (0,1),
        t0bounds = (0.0, 0.1)
    ),

    # Termination criteria
    callback_specs = (
        OutOfBounds(
            x = (0,1),
            z = (0,1)
        ),
        Relativistic(fraction = 0.12)
    ),
)

prob = TraceParticlesProblem(params)
sim = solve(prob, params);
h5_getdataset(get_filename(prob), [:xf, :yf, :weight])

e0, ef = h5_getenrgies(get_filename(prob))
# Convert mean final energy from Joule to electron volt
sum(ef)/length(ef)*TracePartiles.J2eV 

params = setparameters(p; eom=guidingcentreapproximation!, tf=20_000T)
sim = solve(prob, params);

sim = solve(prob; alg=Rosenbrock23(), trajectories=2, reltol=1e-7)


## [Using the [high-level interface](@ref interface-api)](@id interface-ex)
TraceParticles.jl provides a [high level interface](@ref interface-api) to ease systematic running of ensemble simulations.
By defining a `TraceParticlesParameters` instance, e.g.
```@example interface-example
using OrdinaryDiffEq: Tsit5()
using TraceParticles

outdir = mktempdir()

params = TraceParticlesParameters(
    emfield_file = "emfield.itp.jld2",

    tf = 20T,
    charge = electron_charge,
    mass = electron_mass,
    eom = hybridgcafo!,
    npart = 10,
    seed = 42,
    initialeomid = EoMID.FullOrbit,

    # Parameters for saving data
    expname = "demo",
    # Where to save the results
    datadir = outdir,
     # Solve parameters
    alg = Tsit5(),
    reltol = 1e-6,
    abstol = 1e-6,
    maxiters = 10_000,

    prob_func = SampleICsFromMHD(
        tg_file = "temperature.itp.jld2",
        target_distr_file = "particledensity.itp.jld2"
        xbounds = (0,1),
        ybounds = (0,1),
        zbounds = (0,1),
        t0bounds = (t0, tf)
    ),

    # Termination criteria
    callback_specs = (
        OutOfBounds(
            x = (0,1),
            z = (0,1)
        ),
        Relativistic(fraction = 0.12)
    ),
)
```
and assemblying a `TraceParticlesProblem`,
```@example interface-example
prob = TraceParticlesProblem(params)
```
the simulation can be run with
```@example interface-example
sim = solve(prob, params);
h5_getdataset(get_filename(prob), [:xf, :yf, :weight])
```
The `solve(prob::TraceParticlesProblem, p::TraceParticleParameters`-method also computes and stores initial and final energies of the particles to file.
They are fetched with
```@example interface-example
e0, ef = h5_getenrgies(get_filename(prob))
# Convert mean final energy from Joule to electron volt
sum(ef)/length(ef)*TracePartiles.J2eV 
```
The particle ensemble can rerun with adjusted parameters:
```@example interface-example
params = setparameters(p; eom=guidingcentreapproximation!, tf=20_000T)
sim = solve(prob, params);
```
Solver options can be given directly, bypassing `params`:
```@example interface-example
sim = solve(prob; alg=Rosenbrock23(), trajectories=2, reltol=1e-7)
```

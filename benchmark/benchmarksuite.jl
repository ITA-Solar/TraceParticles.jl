# Main benchmark-suite. Contains all benchmarks.
using BenchmarkTools
using DifferentialEquations
using JLD2
using LinearAlgebra
using Random
using TraceParticles
using StaticArrays

suite = BenchmarkGroup()

#------- Benchmarks for physics.jl ---------------------------------------------
suite["physics.jl"] = BenchmarkGroup()

# gradBdrift and gca_drift_and_acceleration
rng = Xoshiro(1)
gradb = SMatrix{3,3,Float64}(rand(rng, 3, 3))
gradexb = SMatrix{3,3,Float64}(rand(rng, 3, 3))
gradB = SVector{3,Float64}(rand(rng, 3))
B = SVector{3,Float64}(rand(rng, 3))
E = SVector{3,Float64}(rand(rng, 3))
vparal = rand(rng)
q = rand(rng)
m = rand(rng)
mu = rand(rng)
b = B / norm(B)
B_inv = 1 / norm(B)
q_inv = 1 / q
suite["physics.jl"]["gradbdrift"] =
    @benchmarkable TraceParticles.gradbdrift(
        $b, $gradB, $mu, $B_inv, $q_inv
    )
suite["physics.jl"]["gca_drift_and_acceleration"] =
    @benchmarkable TraceParticles.gca_drift_and_acceleration(
        $gradb, $gradexb, $gradB, $B, $E, $vparal, $q, $m, $mu
    )

# fieldgradients
struct EMField1
    E0::Float64
    B0::Float64
    d::Float64
end
function (self::EMField1)(x, y, _)
    return SVector(0, 0, self.E0), self.B0 / self.d * SVector(y, x, 0)
end
struct EMField2
    itp::Vector{<:Any}
end
function (self::EMField2)(x, y, z)
    #domm
    ex = self.itp[1](x, y, z)
    ey = self.itp[2](x, y, z)
    ez = self.itp[3](x, y, z)
    bx = self.itp[4](x, y, z)
    by = self.itp[5](x, y, z)
    bz = self.itp[6](x, y, z)
    #efield = SVector{3,Float64}([self.itp[i](x, y, z) for i in 1:3])
    #bfield = SVector{3,Float64}([self.itp[i](x, y, z) for i in 4:6])
    return SVector(ex, ey, ez), SVector(bx, by, bz)
    #efield = [self.itp[i](x, y, z) for i in 1:3]
    #bfield = [self.itp[i](x, y, z) for i in 4:6]
    #return efield, bfield
end
R_analytical = [1.0, 0.0, 0.0]
SR_analytical = SVector{3,Float64}(R_analytical)
emfield_analytical = EMField1(0, 20, 0.02)
R_numerical = [1.6070314f7, 1e6, -6.3503845f6]
SR_numerical = SVector{3,Float64}(R_numerical)
emfield_itp = load_object(joinpath(Base.source_dir(), "cs_BE.jld2"))
emfield_numerical = EMField2(emfield_itp)
suite["physics.jl"]["fieldgradients_analytical"] =
    @benchmarkable TraceParticles.fieldgradients(
        $SR_analytical,
        $emfield_analytical,
    )
suite["physics.jl"]["fieldgradients_numerical"] =
    @benchmarkable TraceParticles.fieldgradients(
        $SR_numerical,
        $emfield_numerical,
    )

#------- Benchmarks for equations_of_motion.jl ---------------------------------
suite["equations_of_motion.jl"] = BenchmarkGroup()

# guidingcentreapproximation!
p_analytical = (
    charge=1.0,
    mass=1.0,
    magneticmoment=0.025,
    electromagneticfield=emfield_analytical,
)
p_numerical = (
    charge=-TraceParticles.e,
    mass=TraceParticles.m_e,
    magneticmoment=1e-13,
    electromagneticfield=emfield_numerical,
)
du_analytical = zeros(4)
du_numerical = zeros(4)
u_analytical0 = [R_analytical; 0.5]
u_numerical0 = [R_numerical; 2e6]
suite["equations_of_motion.jl"]["guidingcentreapproximation!_analytical"] =
    @benchmarkable guidingcentreapproximation!(
        $du_analytical, $u_analytical0, $p_analytical, 0.0
    )
suite["equations_of_motion.jl"]["guidingcentreapproximation!_numerical"] =
    @benchmarkable guidingcentreapproximation!(
        $du_numerical, $u_numerical0, $p_numerical, 0.0
    )

# lorentzforce!
p_numerical_proton = (
    charge=TraceParticles.e,
    mass=TraceParticles.m_p,
    magneticmoment=1e-13,
    electromagneticfield=emfield_numerical,
)
du_analytical_proton = zeros(6)
du_numerical_proton = zeros(6)
u_analytical_proton0 = [R_analytical; -1; 0.0; -0.5]
u_numerical_proton0 = [R_numerical; 2e6; 0.0; 0.0]
suite["equations_of_motion.jl"]["lorentzforce!_analytical"] =
    @benchmarkable lorentzforce!(
        $du_analytical_proton, $u_analytical_proton0, $p_analytical, 0.0
    )
suite["equations_of_motion.jl"]["lorentzforce!_numerical"] =
    @benchmarkable lorentzforce!(
        $du_numerical_proton, $u_numerical_proton0, $p_numerical_proton, 0.0
    )

#------- Benchmarks for an ODEProblem ------------------------------------------
algorithm = Rosenbrock23()
rtol = 3e-8
atol = 1e-0
tspan_analytical = (0, 0.2)
tspan_numerical = (0, 1)
tspan_numerical_proton = (0, 0.001)
prob_analytical = ODEProblem(
    guidingcentreapproximation!,
    u_analytical0,
    tspan_analytical,
    p_analytical
)
prob_numerical = ODEProblem(
    guidingcentreapproximation!,
    u_numerical0,
    tspan_numerical,
    p_numerical
)
prob_analytical_proton = ODEProblem(
    lorentzforce!,
    u_analytical_proton0,
    tspan_analytical,
    p_analytical,
)
prob_numerical_proton = ODEProblem(
    lorentzforce!,
    u_numerical_proton0,
    tspan_numerical_proton,
    p_numerical
)

# ODEProblem for the guiding centre approximation
suite["ODEProblems"]["gca_analytical"] =
    @benchmarkable solve($prob_analytical, $algorithm)
suite["ODEProblems"]["gca_numerical"] =
    @benchmarkable solve($prob_numerical, $algorithm; reltol=$rtol, abstol=$atol)
# ODEProblem for the Lorentz force
suite["ODEProblems"]["lorentzforce!_analytical"] =
    @benchmarkable solve($prob_analytical_proton, $algorithm)
suite["ODEProblems"]["lorentzforce!_numerical"] =
    @benchmarkable solve($prob_numerical_proton, $algorithm; reltol=$rtol, abstol=$atol)

#------- Benchmarks for problem_functions.jl -----------------------------------
xrange = emfield_numerical.itp[1].itp.itp.ranges[1]
zrange = emfield_numerical.itp[1].itp.itp.ranges[2]
xbounds = (xrange[1], xrange[end])
zbounds = (zrange[1], zrange[end])
tg_itp = load_object(joinpath(Base.source_dir(), "cs_tg_normfalse.jld2"))
ne_itp = load_object(joinpath(Base.source_dir(), "cs_ne_normfalse.jld2"))
maxval = maximum(ne_itp.itp.itp.itp.coefs)
probfunc = DposMBvel(
    rng, ne_itp, ne_itp, [xbounds, zbounds], maxval, tg_itp
)
suite["problem_functions.jl"]["DposMBvel"] =
    @benchmarkable probfunc($prob_numerical, 1, false)

#------- Benchmarks for callbacks.jl -------------------------------------------
mass = p_numerical.mass
relfraction = 0.02
gcatol = 0.001
out_cb = DiscreteCallback(
    OutOfDomainCondition_2Dxz(xbounds, zbounds),
    outofdomainaffect!
)
rel_cb = DiscreteCallback(
    RelativisticConditionGCA(; mass=mass, fraction=relfraction),
    relativisticaffect!
)
gca_cb = DiscreteCallback(
    GCABreakDownCondition(gcatol),
    gcabreakdownaffect!,
)
suite["callbacks.jl"]["OutOfDomainCondition_2Dxz"] =
    @benchmarkable out_cb.condition($u_numerical0, 1, $prob_numerical)
suite["callbacks.jl"]["RelativisticConditionGCA"] =
    @benchmarkable rel_cb.condition($u_numerical0, 1, $prob_numerical)
suite["callbacks.jl"]["GCABreakDownCondition"] =
    @benchmarkable gca_cb.condition($u_numerical0, 1, $prob_numerical)

#------- Benchmarks for output_functions.jl ------------------------------------
prob_numerical2 = ODEProblem(
    guidingcentreapproximation!,
    u_numerical0,
    tspan_numerical,
    (
        charge=p_numerical.charge,
        mass=p_numerical.mass,
        magneticmoment=p_numerical.magneticmoment,
        electromagneticfield=p_numerical.electromagneticfield,
        weight=1.0,
        acceptenceratio=1.0,
        userdata=TraceParticles.UserData("None"),
    )
)
sol = solve(prob_numerical2, algorithm; reltol=rtol, abstol=atol)
suite["output_functions.jl"]["output_func_max_lightweight"] =
    @benchmarkable output_func_max_lightweight($sol, 1)

#------- Benchmarks for ensemble problems --------------------------------------

ensemble_prob = EnsembleProblem(
    prob_numerical;
    prob_func=probfunc,
    output_func=output_func_max_lightweight,
    safetycopy=false
)
cbset = CallbackSet(out_cb, rel_cb, gca_cb)
npart = 1000
suite["Ensembles"]["numerical"] =
    @benchmarkable solve(
        $ensemble_prob, $algorithm;
        trajectories=$npart,
        reltol=$rtol,
        abstol=$atol,
        callback=$cbset,
        batch_size=$npart,
    )

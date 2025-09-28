using DifferentialEquations
using Test

using TraceParticles:
    OutOfBoundsCondition,
    MagneticGradientCondition,
    MagneticCurvatureCondition,
    ParallelElectricFieldCondition,
    RelativisticConditionGCA,
    outofboundsaffect!,
    magneticgradientaffect!,
    magneticcurvatureaffect!,
    parallelelectricfieldaffect!,
    relativisticaffect!,
    TerminationCode,
    guidingcentreapproximation!,
    GCAParams

if !isdefined(Main, :verbose)
    verbose = 4
end

@testset verbose = verbose ≥ 4 "callbacks.jl" begin
    xbounds = (0.0,1.0)
    zbounds=(-1.0,0.0)
    u1 = [0.0,1,0.0]
    u2 = [0.1,1, -0.1]
    u3 = [0.1,1, -1.1]
    u4 = [1.1,1,-0.5]
    u5 = ones(4)
    mass = 2.0
    fraction = 0.01
    tol = 0.001
    mass = 1.0
    charge = 2.0
    magneticmoment = 3.0
    t = 0.0
    emfield(x,y,z,t) = [0, x, x^2], [0, 0, x]
    tspan=(0,1)

    out_cb = DiscreteCallback(
        # Periodic boundary conditions in y
        OutOfBoundsCondition( (xbounds, zbounds), [1,3] ), 
        outofboundsaffect!;
    )
    rel_cb = DiscreteCallback(
        RelativisticConditionGCA(;mass=mass, fraction=fraction),
        relativisticaffect!;
    )
    mg_cb = DiscreteCallback(
        MagneticGradientCondition(tol),
        magneticgradientaffect!,
    )
    mc_cb = DiscreteCallback(
        MagneticCurvatureCondition(tol),
        magneticcurvatureaffect!,
    )
    pef_cb = DiscreteCallback(
        ParallelElectricFieldCondition(tol),
        parallelelectricfieldaffect!,
    )

    integrator = ODEProblem(
        guidingcentreapproximation!,
        u5,
        tspan,
        GCAParams(;
            mass=mass,
            charge=charge,
            magneticmoment=magneticmoment,
            electromagneticfield=emfield,
            )
        )
    @test !out_cb.condition(u1, t, integrator)
    @test !out_cb.condition(u2, t, integrator)
    @test out_cb.condition(u3, t, integrator)
    @test out_cb.condition(u4, t, integrator)

    @test !rel_cb.condition(u5, t, integrator)
    @test rel_cb.condition([1,1,1,1e8], t, integrator)
    @test rel_cb.condition([1e15,1,1,1], t, integrator)
    @test pef_cb.condition(u5, t, integrator)
    @test !pef_cb.condition([0.0005, 1.0, 1.0, 1.0, 1.0], t, integrator)
    @test mg_cb.condition(u5, t, integrator)
    @test !mg_cb.condition([200,1.0,1.0,1.0], t, integrator)
    @test !mc_cb.condition(u5, t, integrator)
    @test !mc_cb.condition([1e3,1,1,1], t, integrator)
end

using Test
using TraceParticles: 
    EoMID,
    init_probfunc,
    guidingcentreapproximation!,
    hybridgcafo!,
    HybridScheme,
    lorentzforce!,
    MHDSampler,
    SampleFullOrbit,
    SampleGCA,
    SampleHybrid,
    TraceParticlesParameters

if !isdefined(Main, :verbose)
    verbose = 5
end

function fixture_mhdsampler(; lb=0.0, ub=1.0, tmin=0.0, tmax=0.0)
    temperature(x, y, z, t) = 1e6 # K
    target(x, y, z, t) = 1e-10 # kg/m^3
    return MHDSampler(;
        temperature=temperature,
        target_distr=target,
        xmin=lb, ymin=lb, zmin=lb, tmin=tmin,
        xmax=ub, ymax=ub, zmax=ub, tmax=tmax,
        maxval=1e-10,
    )
end

@testset verbose = verbose ≥ 4 "problem_functions.jl" begin
    baseparams = (;
        tf = 1.0
    )
    makeparams(; kwargs...) = (;
        baseparams..., kwargs...
    )
    dummy_field(_,_,_,_) = zeros(3), zeros(3)

    @testset verbose = verbose ≥ 5 "init_probfunc" begin
        custom = (prob, ctx) -> prob
        @test init_probfunc(custom, makeparams(), dummy_field) === custom

        sampler = fixture_mhdsampler()

        # The equations of motion pick the problem function.
        fo = init_probfunc(sampler, makeparams(eom=lorentzforce!), dummy_field)
        @test fo isa SampleFullOrbit
        @test fo.tf == 1.0
        @test fo.sampler === sampler

        gca = init_probfunc(
            sampler, makeparams(eom=guidingcentreapproximation!),
            dummy_field,
        ) 
        @test gca isa SampleGCA
        @test fo.tf == 1.0
        @test fo.sampler === sampler

        hyb = init_probfunc(
            sampler,
            makeparams(
                eom=hybridgcafo!,
                initialeomid=EoMID.FullOrbit,
                save_switchinfo=false,
            ),
            dummy_field,
        )
        @test hyb.tf == 1.0
        @test hyb.sampler === sampler
        @test hyb isa SampleHybrid
        @test hyb.initialeomid == EoMID.FullOrbit
        # `save_switchinfo` is the single source of truth for the scheme.
        @test hyb.hybridscheme == HybridScheme.Default
        @test init_probfunc(
            sampler, makeparams(
                eom=hybridgcafo!,
                initialeomid=EoMID.FullOrbit,
                save_switchinfo=true
            ),
            dummy_field,
        ).hybridscheme == HybridScheme.SaveSwitches

        notaneom!(du, u, p, t) = nothing
        @test_throws ArgumentError init_probfunc(
            sampler, makeparams(eom=notaneom!), dummy_field
        )
    end
end

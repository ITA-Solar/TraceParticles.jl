using OrdinaryDiffEq: ODEProblem, solve, DiscreteCallback
using TraceParticles

"""
    linearisedharrisheet(x; η, d, g, b)
- `η`: Strength of component normal to the current sheet: η*b
- `d`: Current sheet width: |B/∇B| = d (with η = 0)
- `g`: The guide field strength: g*b
- `B0`: Characteristic field strength
)
A magnetic field configuration with a y-component that reverses at x=0.
Equivalent to the Harris current sheet in the regime x/d << 1.

                   ||| || | |  | | || |||
      y            ||| || | |  | | || |||
      ^            ||| || | |  | | || |||
      |            ||| || | |  | | || |||
      |/           ^^^ ^^ ^ ^  v v vv vvv
     -|-----> x    ||| || | |  | | || |||
     /             ||| || | |  | | || |||   η = 0
    /              ||| || | |  | | || |||
   z               ||| || | |  | | || |||


"""
function linearisedharrisheet(x; η, d, g, b)
    return [η, -x / d, g] * b
end
function xpoint(x, y; d, g, b)
    return [y/d, x/d, g] * b
end
function harrissheet(x; η, d, g, b)
    return [η, tanh(-x / d), g] * b
end


mass = proton_mass
charge = proton_charge

L0 = 1e4 # m
vel0 = [0.0, 0.0, 0.0] #478941.6577971818]
pos0 = [1e-6, 1e-10, 1e-10] * L0

b = 1e-2
η = 0.025 #0.025
d = 1e-4
g = 0.0

t_eject = π * mass / (charge * η * b)
velyτ = -0.5 * 3.0 * a / (η * b)

t0 = 0.0
tf = 1.1 * t_eject
tspan = (t0, tf)

v0 = 1e7
a = 1e-2 * v0 * b
Ez = -a

function currentsheet(x, _, _, _)
    return [0, 0, Ez], linearisedharrisheet(x; η=η, d=d, g=g, b=b)
end


#-------------------------------------------------------------------------------
# CREATE PROBLEM
# Set initial position and velocities

function run()
    prob = ODEProblem(
        lorentzforce!,
        [pos0; vel0],
        tspan,
        (charge=charge, mass=mass, electromagneticfield=currentsheet)
    )
    hprob = ODEProblem(
        hybridgcafo!,
        [pos0; vel0],
        tspan,
        HybridParamsWithDetection(
            charge=charge,
            mass=mass,
            electromagneticfield=currentsheet,
            getphase=(integrator) -> π / 2,
        )
    )
    R, vparal, mu = get_guidingcentre(
        pos0,
        vel0,
        currentsheet(pos0..., t0)[2],
        currentsheet(pos0..., t0)[1],
        charge,
        mass
    )
    gcaprob = ODEProblem(
        guidingcentreapproximation!,
        [R..., vparal],
        tspan,
        (
            charge=charge,
            mass=mass,
            electromagneticfield=currentsheet,
            magneticmoment=mu,
        )
    )
    hcb = DiscreteCallback(
        HybridSwitchCondition(1e-2, 1e-2, Inf, 0.9),
        TraceParticles.hybridswitchaffect_withdetection!
    )

    fosol = solve(prob)
    hsol = solve(hprob, callback=hcb)
    gcasol = solve(gcaprob)
    return fosol, hsol, gcasol
end

function plot(fosol, hsol, gcasol)

    velysimτ = fosol(t_eject)[5]
    hvelysimτ = hsol(t_eject)[5]
    @test isapprox(velysimτ, velyτ, atol=abs(5 * velyτ))
    @test isapprox(hvelysimτ, velyτ, atol=abs(5 * velyτ))

    t = range(t0, tf, 1000)
    #fenergy = [get_observable(fosol, :energy, ti, EoM="FO") for ti in t]
    #henergy = [get_observable(hsol, :energy, ti) for ti in t]

    fu = fosol(range(t0, tf, 10_000))
    hu = hsol(range(t0, tf, 10_000))
    gcau = gcasol(range(t0, tf, 1_000))

    # Plot
    fig = Figure()

    ax = Axis(fig[1,1]; aspect=DataAspect())
    lines!(ax, fu[1, :], fu[2, :]; label="Full orbit")
    lines!(ax, hu[1, :], hu[2, :]; label="Hybrid")
    axislegend(ax; position=:rt)
    ax.xlabel = "x"
    ax.ylabel = "y"

    ax = Axis(fig[1,2]; aspect=DataAspect())
    lines!(ax, fu[2, :], fu[3, :])
    lines!(ax, hu[2, :], hu[3, :])
    ax.xlabel = "y"
    ax.ylabel = "z"

    ax = Axis(fig[2,1]; aspect=DataAspect())
    lines!(ax, fu[1, :], fu[3, :])
    lines!(ax, hu[1, :], hu[3, :])
    ax.xlabel = "x"
    ax.ylabel = "z"

    return fig
end
# fosol, hsol, gcasol = run()
# fig = plot(fosol, hsol, gcasol)

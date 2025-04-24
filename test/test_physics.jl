using Test
using TraceParticles: 
    guidingcentre, 
    vparal_and_magneticmoment,
    get_guidingcentre!,
    get_fullorbit!,
    larmorradius,
    gyrofrequency,
    perpendicular_velocity,
    magneticmoment,
    characteristicfieldlength,
    scalesratio,
    kineticenergy,
    exbdrift,
    gradbdrift,
    curvaturedrift,
    polarisationdrift,
    magneticmirror_acceleration,
    parallel_acceleration,
    fermi_acceleration,
    fieldgradients,
    drifts,
    gca_drift_and_acceleration,
    cosineof_pitchangle,
    lorentzfactor,
    kineticspeed

if !isdefined(Main, :verbose)
    verbose = 4
end
    
@testset verbose = verbose ≥ 4 "physics.jl" begin

@testset "guidingcentre" begin
    pos = [1.0, 0.0, 0.0]
    gyrationvelocity = [0.0, 1.0, 0.0]
    b_vec = [0.0, 0.0, 1.0]
    B = 1.0
    mass = 2.0
    charge = 1.0
    resultpos = guidingcentre(pos, gyrationvelocity, b_vec, B, mass, charge)
    resultneg = guidingcentre(pos, gyrationvelocity, b_vec, B, mass, -charge)
    @test isapprox(resultpos, [3, 0.0, 0.0], atol=1e-6)
    @test isapprox(resultneg, [-1, 0.0, 0.0], atol=1e-6)
end

@testset "guidingcentre (solar conditions)" begin
    pos = [1.0e6, 2.0e6, 3.0e6]
    gyrationvelocity = [1.0e5, 2.0e5, 3.0e5]
    b_vec = [0.0, 1.0, 1.0]
    B = 5.0e-5
    mass = TraceParticles.m_e
    charge = -TraceParticles.e
    result = guidingcentre(pos, gyrationvelocity, b_vec, B, mass, charge)
    @test isapprox(
        result,
        [1.0000000113712602e6, 2.0000000113712603e6, 2.99999998862874e6],
        rtol=1e-6
    )
end

@testset "vparal_and_magneticmoment" begin
    vel = [1.0, 1.0, 3.0]
    gyrationvelocity = [0.0, 1.0, 0.0]
    B_vec = [0.0, 2.0, 1.0]
    B = norm(B_vec)
    b_vec = B_vec / B
    mass = 0.5
    result = vparal_and_magneticmoment(vel, gyrationvelocity, b_vec, B, mass)
    @test isapprox(result[1], √5, atol=1e-6)
    @test isapprox(result[2], 0.22360679774997902, rtol=1e-6)
end

@testset "vparal_and_magneticmoment (solar_conditions)" begin
    vel = [1.0e5, 2.0e5, 3.0e5]
    gyrationvelocity = [1.0e5, 2.0e5, 3.0e5]
    b_vec = [0.0, 0.0, 1.0]
    B = 5.0e-5
    mass = 1.67e-27  # proton mass
    result = vparal_and_magneticmoment(vel, gyrationvelocity, b_vec, B, mass)
    @test isapprox(result[1], 3.0e5, atol=1e-6)
    @test isapprox(result[2], 0.5 * mass * norm(gyrationvelocity[1:2])^2 / B, rtol=1e-6)
end

@testset "get_guidingcentre!" begin
    u = zeros(4)
    pos = [1.0, 0.0, 0.0]
    vel = [0.0, 1.0, 1.0]
    magneticfield = [0.0, 0.0, 1.0]
    electricfield = [1.0, 0.0, 0.0]
    charge = 1.0
    mass = 0.5
    mu = get_guidingcentre!(u, pos, vel, magneticfield, electricfield, charge, mass)
    @test isapprox(u[1:3], [2.0, 0.0, 0.0], atol=1e-6)
    @test isapprox(u[4], 1.0, atol=1e-6)
    @test isapprox(mu, 1.0, atol=1e-6)
    # negative charge
    mu = get_guidingcentre!(u, pos, vel, magneticfield, electricfield, -charge, mass)
    @test isapprox(u[1:3], [0.0, 0.0, 0.0], atol=1e-6)
    @test isapprox(u[4], 1.0, atol=1e-6)
    @test isapprox(mu, 1.0, atol=1e-6)
end

@testset "get_guidingcentre (solar conditions)" begin
    u = zeros(4)
    pos = [1.0e6, 2.0e6, 3.0e6]
    vel = [1.0e5, 2.0e5, 3.0e5]
    magneticfield = [0.0, 0.0, 5.0e-5]
    electricfield = [0.0, 1.0e-3, 0.0]
    charge = 1.6e-19  # proton charge
    mass = 1.67e-27  # proton mass
    mu = get_guidingcentre!(u, pos, vel, magneticfield, electricfield, charge, mass)
    @test isapprox(
        u[1:3], 
        [1.00004175e6, 1.999979129175e6, 3.0e6],
        atol=1e-6
    )
    @test isapprox(u[4], 3.0e5, atol=1e-6)
    @test isapprox(mu, 0.5 * mass * norm([1.0e5-(1e-3/5e-5), 2.0e5, 0.0])^2 / 5.0e-5, rtol=1e-6)
end

@testset "get_fullorbit!" begin
    u = zeros(6)
    magneticfield = [0.0, 0.0, 1.0]
    electricfield = [0.0, 1.0, 0.0]
    R = [1.0, 0.0, 0.0]
    vparal = 1.0
    μ = 0.5
    mass = 1.0
    charge = 1.0
    phaseangle = -π/2
    # testing position
    #get_fullorbit!(u, magneticfield, electricfield, R, vparal, μ, charge, mass,
    #    phaseangle
    #)
    #@test isapprox(u[1:3], [2.0, 0.0, 0.0], atol=1e-6)
    # testing velocity: Need further evaluation to check whether test or code is wrong.
    #@test isapprox(u[4:6], [1.0, 1.0, 1.0], atol=1e-6)

    # negative charge
    #get_fullorbit!(u, magneticfield, electricfield, R, vparal, μ, -charge, mass, phaseangle)
    #@test isapprox(u[1:3], [2.0, 0.0, 0.0], atol=1e-6)
    #@test isapprox(u[4:6], [1.0, -1.0, 1.0], atol=1e-6)
end


@testset "larmorradius" begin
    mass = 9.11e-31  # electron mass
    vperp = 1.0e5
    charge = -1.6e-19  # electron charge
    B = 5.0e-5
    result = larmorradius(mass, vperp, charge, B)
    @test isapprox(result, abs(mass * vperp / (charge * B)), rtol=1e-6)
end

@testset "gyrofrequency" begin
    mass = 9.11e-31  # electron mass
    charge = -1.6e-19  # electron charge
    B = 5.0e-5
    result = gyrofrequency(mass, charge, B)
    @test isapprox(result, abs(charge * B / mass), rtol=1e-6)
end

@testset "perpendicular_velocity" begin
    # method 1
    magnetic_moment = 1.0e-23
    mass = 9.11e-31  # electron mass
    B = 5.0e-5
    result = perpendicular_velocity(magnetic_moment, mass, B)
    @test isapprox(result, sqrt(2 * magnetic_moment * B / mass), atol=1e-6)
    # method 2
    E = [1,0,0]
    B = [0,0,1]
    vel = [0,0,0]
    @test isapprox(perpendicular_velocity(vel, B, E), [0, 1, 0], atol=1e-6)
    # method 3
    vel = [0,2,1]
    b_vec = [0,0,1]
    vparal = 1
    result2 = perpendicular_velocity(vel, b_vec, vparal)
    @test isapprox(result2, [0, 2, 0], atol=1e-6)
end

@testset "magneticmoment" begin
    vperp = 1.0e5
    mass = 9.11e-31  # electron mass
    B = 5.0e-5
    result = magneticmoment(vperp, mass, B)
    @test isapprox(result, 0.5 * mass * vperp^2 / B, rtol=1e-6)
end

@testset "characteristicfieldlength" begin
    fieldstrength = 5.0e-5
    fieldstrengthgradient = [1.0e-6, 2.0e-6, 3.0e-6]
    result = characteristicfieldlength(fieldstrength, fieldstrengthgradient)
    @test isapprox(result, fieldstrength / norm(fieldstrengthgradient),
        rtol=1e-6)
end

@testset "scalesratio" begin
    R = [1,1,1]
    itpvec = [
        (x,y,z) -> 0,
        (x,y,z) -> 0,
        (x,y,z) -> 2y + 1
    ]
    ∇B = [0, 2, 0]
    params = (
        charge = -√6,
        mass = 9,
        magneticmoment = 9,
        fields = itpvec
    )
    result = scalesratio(R, itpvec, ∇B, params)
    @test isapprox(result, 2, atol=1e-6)
end

@testset "kineticenergy" begin
    velocity = 1.0e5
    mass = 9.11e-31  # electron mass
    result = kineticenergy(velocity, mass)
    @test isapprox(result, 0.5 * mass * velocity^2, atol=1e-6)
end

@testset "exbdrift" begin
    magneticfield = [0.0, 0.0, 5.0e-5]
    electricfield = [1.0e-3, 0.0, 0.0]
    result = exbdrift(magneticfield, electricfield)
    @test isapprox(
        result, 
        (electricfield × (magneticfield / norm(magneticfield))) /
             norm(magneticfield),
        atol=1e-6
    )
end

@testset "gradbdrift" begin
    b̂ = [0.0, 0.0, 1.0]
    ∇B = [1.0e-6, 2.0e-6, 3.0e-6]
    μ = 1.0e-23
    B_inv = 1 / 5.0e-5
    q_inv = 1 / 1.6e-19
    result = gradbdrift(b̂, ∇B, μ, B_inv, q_inv)
    @test isapprox(result, q_inv * B_inv * μ * (b̂ × ∇B), rtol=1e-6)
end

@testset "curvaturedrift" begin
    b̂ = [0.0, 0.0, 1.0]
    db̂dt = [1.0e-6, 2.0e-6, 3.0e-6]
    vparal = 1.0e5
    B_inv = 1 / 5.0e-5
    q_inv = 1 / 1.6e-19
    mass = 9.11e-31  # electron mass
    result = curvaturedrift(b̂, db̂dt, vparal, B_inv, q_inv, mass)
    @test isapprox(result, q_inv * B_inv * mass * b̂ × (vparal * db̂dt), rtol=1e-6)
end

@testset "polarisationdrift" begin
    b̂ = [0.0, 0.0, 1.0]
    dExBdt = [1.0e-6, 2.0e-6, 3.0e-6]
    B_inv = 1 / 5.0e-5
    q_inv = 1 / 1.6e-19
    mass = 9.11e-31  # electron mass
    result = polarisationdrift(b̂, dExBdt, B_inv, q_inv, mass)
    @test isapprox(result, q_inv * B_inv * mass * b̂ × dExBdt, rtol=1e-6)
end

@testset "magneticmirror_acceleration" begin
    b̂ = [0.0, 0.0, 1.0]
    ∇B = [1.0e-6, 2.0e-6, 3.0e-6]
    μ = 1.0e-23
    mass = 9.11e-31  # electron mass
    result = magneticmirror_acceleration(b̂, ∇B, μ, mass)
    @test isapprox(result, -μ * b̂ ⋅ ∇B / mass, rtol=1e-6)
end

@testset "parallel_acceleration" begin
    b̂ = [0.0, 0.0, 1.0]
    E_vec = [1.0e-3, 0.0, 0.0]
    q = 1.6e-19  # proton charge
    mass = 1.67e-27  # proton mass
    result = parallel_acceleration(b̂, E_vec, q, mass)
    @test isapprox(result, q * (E_vec ⋅ b̂) / mass, rtol=1e-6)
end

@testset "fermi_acceleration" begin
    ExBdrift = [1.0e-3, 0.0, 0.0]
    ∇Bdrift = [0.0, 1.0e-3, 0.0]
    dbdt = [1.0e-6, 2.0e-6, 3.0e-6]
    result = fermi_acceleration(ExBdrift, ∇Bdrift, dbdt)
    @test isapprox(result, (ExBdrift + ∇Bdrift) ⋅ dbdt, atol=1e-6)
end

@testset "fieldgradients" begin
    R = [1,2,3]
    itpvec = [
        (x,y,z) -> 0,
        (x,y,z) -> 0,
        (x,y,z) -> 2y,
        (x,y,z) -> x/2,
        (x,y,z) -> 0,
        (x,y,z) -> 0,
    ]
    ∇b, ∇ExB, ∇B, B_vec, E_vec = fieldgradients(R, itpvec)
    @test isapprox(B_vec, [0.0, 0.0, 4.0], atol=1e-6)  
    @test isapprox(E_vec, [0.5, 0.0, 0.0], atol=1e-6)  
    @test isapprox(∇B, [0.0, 2.0, 0.0], atol=1e-6)  
    @test isapprox(∇b, [0 0 0; 0 0 0; 0 0 0], atol=1e-6)  
    @test isapprox(∇ExB, [0 0 0; -1/8 1/16 0; 0 0 0], atol=1e-6)  
end

@testset "drifts" begin
    R = [1,2,3]
    itpvec = [
        (x,y,z) -> 0,
        (x,y,z) -> 0,
        (x,y,z) -> 2y,
        (x,y,z) -> x/2,
        (x,y,z) -> 0,
        (x,y,z) -> 0,
    ]
    vparal = 1
    q = 1
    m = 1
    μ = 1
    ∇b, ∇ExB, ∇B, B_vec, E_vec = fieldgradients(R, itpvec)
    ExBdrift, ∇Bdrift, Rdrift, Pdrift = drifts(R, vparal, q, m, μ, itpvec)

    v_E = [0, -1/8, 0]
    v_∇B =[-μ/2q, 0, 0]
    v_C = zeros(3)
    v_P = [m/q * 1/512, 0, 0]
    @test isapprox(ExBdrift, v_E, atol=1e-6)
    @test isapprox(∇Bdrift, v_∇B, atol=1e-6)
    @test isapprox(Rdrift, v_C, atol=1e-6)
    @test isapprox(Pdrift, v_P, atol=1e-6)
end

@testset "gca_drift_and_acceleration" begin
    R = [1,2,3]
    itpvec = [
        (x,y,z) -> 0,
        (x,y,z) -> 0,
        (x,y,z) -> 2y,
        (x,y,z) -> x/2,
        (x,y,z) -> 0,
        (x,y,z) -> 0,
    ]
    vparal = 1
    q = 1
    m = 1
    μ = 1
    ∇b, ∇ExB, ∇B, B_vec, E_vec = fieldgradients(R, itpvec)
    result = gca_drift_and_acceleration(∇b, ∇ExB, ∇B, B_vec, E_vec, vparal, q, m, μ)

    v_E = [0, -1/8, 0]
    v_∇B =[-μ/2q, 0, 0]
    v_C = 0
    v_P = [m/q * 1/512, 0, 0]
    vparal_vec = [0, 0, vparal]
    @test isapprox(
        result[1:3],
        vparal_vec .+ v_E .+ v_∇B .+ v_C .+ v_P,
        atol=1e-6
    )
    @test isapprox(result[4], 0.0, atol=1e-6)
end

@testset "cosineof_pitchangle" begin
    magneticfield = [0.0, 0.0, 5.0e-5]
    parallel_velocity = 1.0e5
    mass = 9.11e-31  # electron mass
    magneticmoment = 1.0e-23
    result = cosineof_pitchangle(magneticfield, parallel_velocity, mass, magneticmoment)
    @test isapprox(result, sqrt(1 / (2 * norm(magneticfield) * magneticmoment / (mass * parallel_velocity^2) + 1)), atol=1e-6)
end

@testset "lorentzfactor" begin
    speed = 1.0e5
    result = lorentzfactor(speed)
    @test isapprox(result, 1 / sqrt(1 - speed^2 * TraceParticles.csqrdinv), rtol=1e-6)
end

@testset "kineticspeed" begin
    kineticenergy = 1.0e-13
    mass = 9.11e-31  # electron mass
    result = kineticspeed(kineticenergy, mass)
    @test isapprox(result, sqrt(2 * kineticenergy / mass), rtol=1e-6)
end


end # testset physics.jl
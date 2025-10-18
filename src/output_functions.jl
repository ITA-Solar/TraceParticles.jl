"""
    output_func_max_lightweight(sol, i)
Output function for ensemble simulations using `DifferentialEquations`.

Finds the maximum larmor radius, minimum field length and maximum scales ratio
of the particle's trajectory using the inherent interpolation function in the
particle solution. Also returns the initial and final state of the particle,
in addition to its charge, mass and magnetic moment.

We set the required "rerun" return argument to false.
"""
function output_func_max_lightweight(sol, i)
    ntimes = length(sol.t)
    maxrl, meanrl = find_max_larmorradius(sol, ntimes=ntimes)
    minlb, meanlb = find_min_fieldlength(sol, ntimes=ntimes)
    maxscalesratio, meanscalesratio = find_max_scalesratio(
        sol,
        ntimes=ntimes,
    )
    x0, y0, z0, vparal0 = first(sol)
    xf, yf, zf, vparalf = last(sol)
    # Construct and retunr the output tuple
    return (
        x0=x0, y0=y0, z0=z0, vparal0=vparal0,
        xf=xf, yf=yf, zf=zf, vparalf=vparalf,
        charge=sol.prob.p.charge,
        mass=sol.prob.p.mass,
        magneticmoment=sol.prob.p.magneticmoment,
        weight=sol.prob.p.weight,
        nrejections=sol.prob.p.nrejections,
        t0=first(sol.t),
        tf=last(sol.t),
        nt=length(sol.t),
        maxrl=maxrl.max,
        maxrl_x=maxrl.max_u[1],
        maxrl_y=maxrl.max_u[2],
        maxrl_z=maxrl.max_u[3],
        maxrl_vparal=maxrl.max_u[4],
        maxrl_t=maxrl.max_t,
        minlb=minlb.min,
        minlb_x=minlb.min_u[1],
        minlb_y=minlb.min_u[2],
        minlb_z=minlb.min_u[3],
        minlb_vparal=minlb.min_u[4],
        minlb_t=minlb.min_t,
        maxscalesratio=maxscalesratio.max,
        maxscalesratio_x=maxscalesratio.max_u[1],
        maxscalesratio_y=maxscalesratio.max_u[2],
        maxscalesratio_z=maxscalesratio.max_u[3],
        maxscalesratio_vparal=maxscalesratio.max_u[4],
        maxscalesratio_t=maxscalesratio.max_t,
        meanrl=meanrl,
        meanlb=meanlb,
        meanscalesratio=meanscalesratio,
        retcode=Int(sol.retcode),
        terminationcode=Int(sol.prob.p.terminationcode)
    ),
    false
end


"""
    output_func_lightweight(sol, i)
Output function for ensemble simulations using `DifferentialEquations`.

Returns only the initial and final state of the particle, in addition to its
charge, mass and magnetic moment. We set the required "rerun" return argument
to false.
"""
function output_func_lightweight(sol, i)
    x0, y0, z0, vparal0 = first(sol)
    xf, yf, zf, vparalf = last(sol)
    return (
        x0=x0, y0=y0, z0=z0, vparal0=vparal0,
        xf=xf, yf=yf, zf=zf, vparalf=vparalf,
        nt=length(sol.t),
        t0=first(sol.t),
        tf=last(sol.t),
        charge=sol.prob.p.charge,
        mass=sol.prob.p.mass,
        magneticmoment=sol.prob.p.magneticmoment,
        weight=sol.prob.p.weight,
        nrejections=sol.prob.p.nrejections,
        retcode=Int(sol.retcode),
        terminationcode=Int(sol.prob.p.terminationcode)
    ),
    false
end
function output_func_lightweight_hybrid(sol, i)
    x0, y0, z0, vx0, vy0, vz0 = first(sol)
    xf, yf, zf, vxf, vyf, vzf = last(sol)
    return (
        x0=x0, y0=y0, z0=z0, vx0=vx0, vy0=vy0, vz0=vz0,
        xf=xf, yf=yf, zf=zf, vxf=vxf, vyf=vyf, vzf=vzf,
        nt=length(sol.t),
        t0=first(sol.t),
        tf=last(sol.t),
        charge=sol.prob.p.charge,
        mass=sol.prob.p.mass,
        initialmagneticmoment=sol.prob.p.initialmagneticmoment,
        magneticmoment=sol.prob.p.magneticmoment,
        weight=sol.prob.p.weight,
        nrejections=sol.prob.p.nrejections,
        retcode=Int(sol.retcode),
        terminationcode=Int(sol.prob.p.terminationcode),
        nswitches=sol.prob.p.nswitches,
        eomid=sol.prob.p.eomid,
        initialeomid=sol.prob.p.initialeomid,
    ),
    false
end

function out2(sol, i)
    x0, y0, z0, vx0, vy0, vz0 = first(sol)
    xf, yf, zf, vxf, vyf, vzf = last(sol)
    t0, tf = first(sol.t), last(sol.t)
    keys = (
        :x0, :y0, :z0, :vx0, :vy0, :vz0,
        :xf, :yf, :zf, :vxf, :vyf, :vzf,
        :nt,
        :t0, :tf,
        :charge, :mass, :magneticmoment, :weight,
        :nrejections,
        :retcode, :terminationcode,
        :nswitches,
        :eomid, :initialeomid,
    )
    values = (
        x0, y0, z0, vx0, vy0, vz0,
        xf, yf, zf, vxf, vyf, vzf,
        length(sol.t),
        t0, tf,
        sol.prob.p.charge, sol.prob.p.mass, sol.prob.p.magneticmoment, sol.prob.p.weight,
        sol.prob.p.nrejections,
        Int(sol.retcode), Int(sol.prob.p.terminationcode),
        sol.prob.p.nswitches,
        sol.prob.p.eomid, sol.prob.p.initialeomid,
    )
    return NamedTuple{keys}(values),
    false
end

function output_func_preallocated(sol, i)
    x0, y0, z0, vx0, vy0, vz0 = first(sol)
    xf, yf, zf, vxf, vyf, vzf = last(sol)
    t0, tf = first(sol.t), last(sol.t)
    sol.prob.p.out[1, i] = x0
    sol.prob.p.out[2, i] = y0
    sol.prob.p.out[3, i] = z0
    sol.prob.p.out[4, i] = vx0
    sol.prob.p.out[5, i] = vy0
    sol.prob.p.out[6, i] = vz0
    sol.prob.p.out[7, i] = xf
    sol.prob.p.out[8, i] = yf
    sol.prob.p.out[9, i] = zf
    sol.prob.p.out[10, i] = vxf
    sol.prob.p.out[11, i] = vyf
    sol.prob.p.out[12, i] = vzf
    sol.prob.p.out[13, i] = length(sol.t)
    sol.prob.p.out[14, i] = t0
    sol.prob.p.out[15, i] = tf
    sol.prob.p.out[16, i] = sol.prob.p.charge
    sol.prob.p.out[17, i] = sol.prob.p.mass
    return nothing, false
end

struct OutputFuncSaveatHybrid{T1<:AbstractVector,T2<:NTuple{N,Symbol} where {N}}
    saveat::T1
    intermediate_saveats::T1
    outputkeys::T2

    function OutputFuncSaveatHybrid(saveat::T1) where {T1<:AbstractVector}
        N = length(saveat)
        Nintermediate = 1:N-2
        tkeys = Tuple(Symbol("t", i) for i in Nintermediate)
        xkeys = Tuple(Symbol("x", i) for i in Nintermediate)
        ykeys = Tuple(Symbol("y", i) for i in Nintermediate)
        zkeys = Tuple(Symbol("z", i) for i in Nintermediate)
        vxkeys = Tuple(Symbol("vx", i) for i in Nintermediate)
        vykeys = Tuple(Symbol("vy", i) for i in Nintermediate)
        vzkeys = Tuple(Symbol("vz", i) for i in Nintermediate)
        outputkeys = (
            :x0, xkeys..., :xf,
            :y0, ykeys..., :yf,
            :z0, zkeys..., :zf,
            :vx0, vxkeys..., :vxf,
            :vy0, vykeys..., :vyf,
            :vz0, vzkeys..., :vzf,
            :nt, :t0, tkeys..., :tf,
            :charge, :mass, :magneticmoment, :weight, :nrejections,
            :retcode, :terminationcode,
            :nswitches, :eomid, :initialeomid,
        )
        new{T1,typeof(outputkeys)}(saveat, saveat[2:end-1], outputkeys)
    end
end
function (self::OutputFuncSaveatHybrid)(sol, i)
    ntimes = length(sol.t)
    x0, y0, z0, vx0, vy0, vz0 = first(sol)
    xf, yf, zf, vxf, vyf, vzf = last(sol)
    t0, tf = first(sol.t), last(sol.t)
    mask = Tuple(
        findfirst(x -> x == t, sol.t) for t in self.intermediate_saveats
    )
    tvals = Tuple(
        isnothing(i) ? missing : sol.t[i] for i in mask
    )
    xvals = Tuple(
        isnothing(i) ? missing : sol[i][1] for i in mask
    )
    yvals = Tuple(
        isnothing(i) ? missing : sol[i][2] for i in mask
    )
    zvals = Tuple(
        isnothing(i) ? missing : sol[i][3] for i in mask
    )
    vxvals = Tuple(
        isnothing(i) ? missing : sol[i][4] for i in mask
    )
    vyvals = Tuple(
        isnothing(i) ? missing : sol[i][5] for i in mask
    )
    vzvals = Tuple(
        isnothing(i) ? missing : sol[i][6] for i in mask
    )
    outputvalues = (
        x0, xvals..., xf,
        y0, yvals..., yf,
        z0, zvals..., zf,
        vx0, vxvals..., vxf,
        vy0, vyvals..., vyf,
        vz0, vzvals..., vzf,
        ntimes, t0, tvals..., tf,
        sol.prob.p.charge, sol.prob.p.mass, sol.prob.p.magneticmoment, sol.prob.p.weight,
        sol.prob.p.nrejections,
        Int(sol.retcode), Int(sol.prob.p.terminationcode),
        sol.prob.p.nswitches,
        sol.prob.p.eomid, sol.prob.p.initialeomid,
    )
    return NamedTuple{self.outputkeys}(outputvalues), false
end

function out3(sol, i)
    x0, y0, z0, vx0, vy0, vz0 = first(sol)
    xf, yf, zf, vxf, vyf, vzf = last(sol)
    return (
        x0=x0, y0=y0, z0=z0, vx0=vx0, vy0=vy0, vz0=vz0,
        xf=xf, yf=yf, zf=zf, vxf=vxf, vyf=vyf, vzf=vzf,
        nt=length(sol.t),
        t0=first(sol.t),
        tf=last(sol.t),
        charge=sol.prob.p.charge,
        mass=sol.prob.p.mass,
        magneticmoment=sol.prob.p.magneticmoment,
        weight=sol.prob.p.weight,
        nrejections=sol.prob.p.nrejections,
        retcode=Int(sol.retcode),
        terminationcode=Int(sol.prob.p.terminationcode),
        nswitches=sol.prob.p.nswitches,
        eomid=sol.prob.p.eomid,
        initialeomid=sol.prob.p.initialeomid,
        x=sol[:, 1],
        y=sol[:, 2],
        z=sol[:, 3],
        vx=sol[:, 4],
        vy=sol[:, 5],
        vz=sol[:, 6],
        t=sol.t,
    ),
    false
end

"""
    find_max_larmorradius(sol; ntimes=5length(sol.t))
Find the maximum larmor radius of the particle's trajectory using the inherent
interpolation function in the particle solution. Also returns the mean larmor
radius. The solution is evaluated at `ntimes` points in time.

`sol` should return the solution at a given time `t` with the syntax `sol(t)`,
and the initial and final times must be accessible with `first(sol.t)` and
`last(sol.t)`, respectively.
"""
function find_max_larmorradius(sol; ntimes=5length(sol.t))
    f(t) = larmorradius(
        sol(first(t), idxs=1:3),
        t,
        sol.prob.p.magneticmoment,
        sol.prob.p.charge,
        sol.prob.p.mass,
        sol.prob.p.electromagneticfield,
    )
    if ntimes == 1
        time = sol.t[1]
        farray = larmorradius(
            sol[1][1:3],
            time,
            sol.prob.p.magneticmoment,
            sol.prob.p.charge,
            sol.prob.p.mass,
            sol.prob.p.electromagneticfield,
        )
        max = MaxValue(farray, sol[1], time)
        mean = sum(farray) / length(farray)
    else
        times = range(first(sol.t), last(sol.t), length=ntimes)
        farray = [f(t) for t in times]
        maxidx = argmax(farray)
        max = MaxValue(farray[maxidx], sol(times[maxidx]), times[maxidx])
        mean = sum(farray) / length(farray)
    end
    return max, mean
end


"""
    find_min_fieldlength(sol; ntimes=5length(sol.t))
Find the minimum field length of the particle's trajectory using the inherent
interpolation function in the particle solution. Also returns the mean field
length. The solution is evaluated at `ntimes` points in time.

`sol` should return the solution at a given time `t` with the syntax `sol(t)`,
and the initial and final times must be accessible with `first(sol.t)` and
`last(sol.t)`, respectively.
"""
function find_min_fieldlength(sol; ntimes=5length(sol.t))
    f(t) = characteristicfieldlength(
        sol(first(t), idxs=1:3), t,
        sol.prob.p.electromagneticfield,
    )
    if ntimes == 1
        time = sol.t[1]
        farray = characteristicfieldlength(
            sol[1][1:3],
            time,
            sol.prob.p.electromagneticfield,
        )
        min = MinValue(farray, sol[1], time)
        mean = farray
    else
        times = range(first(sol.t), last(sol.t), length=ntimes)
        farray = [f(t) for t in times]
        minidx = argmin(farray)
        min = MinValue(farray[minidx], sol(times[minidx]), times[minidx])
        mean = sum(farray) / length(farray)
    end
    return min, mean
end


"""
    find_max_scalesratio(sol; ntimes=5length(sol.t))
Find the maximum scales ratio of the particle's trajectory using the inherent
interpolation function in the particle solution. Also returns the mean scales
ratio. The solution is evaluated at `ntimes` points in time.

`sol` should return the solution at a given time `t` with the syntax `sol(t)`,
and the initial and final times must be accessible with `first(sol.t)` and
`last(sol.t)`, respectively.
"""
function find_max_scalesratio(sol; ntimes=5length(sol.t))
    f(t) = scalesratio(
        sol(first(t), idxs=1:3),
        t,
        sol.prob.p.mass,
        sol.prob.p.charge,
        sol.prob.p.magneticmoment,
        sol.prob.p.electromagneticfield,
    )
    if ntimes == 1
        time = sol.t[1]
        farray = scalesratio(
            time,
            sol[1][1:3],
            sol.prob.p.mass,
            sol.prob.p.charge,
            sol.prob.p.magneticmoment,
            sol.prob.p.electromagneticfield,
        )
        max = MaxValue(farray, sol[1], time)
        mean = farray
    else
        times = range(first(sol.t), last(sol.t), length=ntimes)
        farray = [f(t) for t in times]
        maxidx = argmax(farray)
        max = MaxValue(farray[maxidx], sol(times[maxidx]), times[maxidx])
        mean = sum(farray) / length(farray)
    end
    return max, mean
end

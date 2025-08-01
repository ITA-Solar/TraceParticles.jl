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
        retmsg=sol.prob.p.retcode,
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
        tf=last(sol.t),
        charge=sol.prob.p.charge,
        mass=sol.prob.p.mass,
        magneticmoment=sol.prob.p.magneticmoment,
        weight=sol.prob.p.weight,
        nrejections=sol.prob.p.nrejections,
        retcode=Int(sol.retcode),
        retmsg=sol.prob.p.retcode
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
        sol.prob.p.magneticmoment,
        sol.prob.p.charge,
        sol.prob.p.mass,
        sol.prob.p.electromagneticfield,
    )
    if ntimes == 1
        time = sol.t[1]
        farray = larmorradius(
            sol[1][1:3],
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
        sol(first(t), idxs=1:3),
        sol.prob.p.electromagneticfield,
    )
    if ntimes == 1
        time = sol.t[1]
        farray = characteristicfieldlength(
            sol[1][1:3],
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
        sol.prob.p.mass,
        sol.prob.p.charge,
        sol.prob.p.magneticmoment,
        sol.prob.p.electromagneticfield,
    )
    if ntimes == 1
        time = sol.t[1]
        farray = scalesratio(
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

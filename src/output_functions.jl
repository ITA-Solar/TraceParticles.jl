"""
    maxrl_minlb_maxscalesratio(sol, i)
Output function for ensemble simulations using `DifferentialEquations`.

Finds the maximum larmor radius, minimum field length and maximum scales ratio
of the particle's trajectory using the inherent interpolation function in the
particle solution. Also returns the full particle solution.
"""
function output_func_max_full(sol, i)
    ntimes = length(sol.t)
    maxrl = find_max_larmorradius(sol, ntimes=ntimes)
    minlb = find_min_fieldlength(sol, ntimes=ntimes)
    maxscalesratio = find_max_scalesratio(sol, ntimes=ntimes)
    return (
        sol=sol,
        larmorradius=maxrl,
        fieldlength=minlb,
        scalesratio=maxscalesratio
        ),
        false
end


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
    maxscalesratio, meanscalesratio = find_max_scalesratio(sol, ntimes=ntimes)
    x0, y0, z0, vparal0 = first(sol)
    xf, yf, zf, vparalf = last(sol)
    # Select parameters to save
    exclusionlist = [:fields, :temp, :userdata]
    param_to_save = filter(x -> !(x in exclusionlist), keys(sol.prob.p))
    # Construct and retunr the output tuple
    return (
        x0=x0, y0=y0, z0=z0, vparal0=vparal0,
        xf=xf, yf=yf, zf=zf, vparalf=vparalf,
        sol.prob.p[param_to_save]...,
        t0 = first(sol.t),
        tf = last(sol.t),
        nt = length(sol.t),
        maxrl = maxrl.max,
        maxrl_x = maxrl.max_u[1],
        maxrl_y = maxrl.max_u[2],
        maxrl_z = maxrl.max_u[3],
        maxrl_vparal = maxrl.max_u[4],
        maxrl_t = maxrl.max_t,
        minlb = minlb.min,
        minlb_x = minlb.min_u[1],
        minlb_y = minlb.min_u[2],
        minlb_z = minlb.min_u[3],
        minlb_vparal = minlb.min_u[4],
        minlb_t = minlb.min_t,
        maxscalesratio = maxscalesratio.max,
        maxscalesratio_x = maxscalesratio.max_u[1],
        maxscalesratio_y = maxscalesratio.max_u[2],
        maxscalesratio_z = maxscalesratio.max_u[3],
        maxscalesratio_vparal = maxscalesratio.max_u[4],
        maxscalesratio_t = maxscalesratio.max_t,
        meanrl = meanrl,
        meanlb = meanlb,
        meanscalesratio = meanscalesratio,
        retcode = string(sol.retcode),
        retmsg = sol.prob.p.userdata.retmsg,
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
    u0 = first(sol)
    uf = last(sol)
    return (
        u0=u0,
        uf=uf,
        charge = sol.prob.p.charge,
        mass = sol.prob.p.mass,
        magneticmoment=sol.prob.p.magneticmoment,
        ),
        false
end

function find_max_larmorradius(sol; ntimes=5length(sol.t))
    f(t) = larmorradius(sol(first(t))[1:2:3], sol.prob.p)
    if ntimes == 1
        time = sol.t[1]
        farray = larmorradius(sol[1][1:2:3], sol.prob.p)
        max = MaxValue(farray, sol[1], time)
        mean = sum(farray)/length(farray)
    else
        times = range(first(sol.t), last(sol.t), length=ntimes)
        farray = [f(t) for t in times]
        maxidx = argmax(farray)
        max = MaxValue(farray[maxidx], sol(times[maxidx]), times[maxidx])
        mean = sum(farray)/length(farray)
    end
    return max, mean
end


function find_min_fieldlength(sol; ntimes=5length(sol.t))
    f(t) = characteristicfieldlength(
        sol(first(t))[1:2:3], 
        sol.prob.p.fields[1:3]
        )
    if ntimes == 1
        time = sol.t[1]
        farray = characteristicfieldlength(
            sol[1][1:2:3], sol.prob.p.fields[1:3]
            )
        min = MinValue(farray, sol[1], time)
        mean = farray
    else
        times = range(first(sol.t), last(sol.t), length=ntimes)
        farray = [f(t) for t in times]
        minidx = argmin(farray)
        min = MinValue(farray[minidx], sol(times[minidx]), times[minidx])
        mean = sum(farray)/length(farray)
    end
    return min, mean
end


function find_max_scalesratio(sol; ntimes=5length(sol.t))
    frl(t) = larmorradius(sol(first(t))[1:2:3], sol.prob.p)
    flb(t) = characteristicfieldlength(
        sol(first(t))[1:2:3], 
        sol.prob.p.fields[1:3]
        )
    f(t) = frl(t)/flb(t)
    if ntimes == 1
        time = sol.t[1]
        farray = larmorradius(sol[1][1:2:3], sol.prob.p)/
            characteristicfieldlength(sol[1][1:2:3], sol.prob.p.fields[1:3])
        max = MaxValue(farray, sol[1], time)
        mean = farray
    else
        times = range(first(sol.t), last(sol.t), length=ntimes)
        farray = [f(t) for t in times]
        maxidx = argmax(farray)
        max = MaxValue(farray[maxidx], sol(times[maxidx]), times[maxidx])
        mean = sum(farray)/length(farray)
    end
    return max, mean
end


"""
    out_sol_maxrl_opt(sol, i)
Output function for ensemble simulations using `DifferentialEquations`. 

Finds the maximum larmor radius during the particle's trajectory using
`Optimization`. Also returns the full particle solution.
"""
function out_sol_maxrl_opt(sol, i)
    fminus(t,_) =  -larmorradius(sol(first(t))[1:2:3], sol.prob.p)
    optf = OptimizationFunction(fminus, Optimization.AutoForwardDiff())
    min_guess = (first(sol.t) + last(sol.t))/2
    optprob = OptimizationProblem(
        optf, [min_guess], lb = [first(sol.t)], ub = [last(sol.t)]
        )
    opt = Optimization.solve(optprob, NLopt.GN_ORIG_DIRECT_L())
    maxrl = MaxValue(-opt.minimum, sol(opt.u), opt.u)
    return (
        sol=sol,
        larmorradius=maxrl
        ),
        false
end

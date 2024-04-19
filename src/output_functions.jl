
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

function out_sol_maxrl(sol, i)
    f(t) =  larmorradius(sol(first(t))[1:2:3], sol.prob.p)
    ntimes = 2000
    times = range(first(sol.t), last(sol.t), length=ntimes)
    farray = [f(t) for t in times]
    maxidx = argmax(farray)
    maxrl = MaxValue(farray[maxidx], sol(times[maxidx]), times[maxidx])
    return (
        sol=sol,
        larmorradius=maxrl
        ),
        false
end


function out_sol_maxrllbscalesratio(sol, i)
    frl(t) = larmorradius(sol(first(t))[1:2:3], sol.prob.p)
    flb(t) = characteristicfieldlength(
        sol(first(t))[1:2:3], 
        sol.prob.p.fields[1:3]
        )
    fscalesratio(t) = frl(t)/flb(t)

    ntimes = 5length(sol.t)
    times = range(first(sol.t), last(sol.t), length=ntimes)

    frl_array = [frl(t) for t in times]
    flb_array = [flb(t) for t in times]
    fscalesratio_array = [fscalesratio(t) for t in times]

    maxidx_rl = argmax(frl_array)
    maxrl = MaxValue(
        frl_array[maxidx_rl],
        sol(times[maxidx_rl]),
        times[maxidx_rl]
        )

    minidx_lb = argmin(flb_array)
    minlb = MinValue(
        flb_array[minidx_lb],
        sol(times[minidx_lb]),
        times[minidx_lb]
        )

    maxidx_scalesratio = argmax(fscalesratio_array)
    maxscalesratio = MaxValue(
        fscalesratio_array[maxidx_scalesratio],
        sol(times[maxidx_scalesratio]),
        times[maxidx_scalesratio]
        )
    return (
        sol=sol,
        larmorradius=maxrl,
        fieldlength=minlb,
        scalesratio=maxscalesratio
        ),
        false
end

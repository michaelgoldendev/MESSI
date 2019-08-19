using NLopt

function myfunc(x::Vector)
    ret = 2.0*x[1]*x[1]*x[2] + x[2]*x[3] - x[3]
    println(sum(x),"\t",ret)
    return ret
end

function sumlessthan(x::Vector, maxsum::Float64=0.95)
    return x[12] + x[13] + x[14] - maxsum
end

opt = Opt(:LN_COBYLA, 3)
lower_bounds!(opt, ones(Float64,3)*1e-4)
upper_bounds!(opt, ones(Float64,3)*0.999) 
opt.xtol_rel = 1e-4
inequality_constraint!(opt, (x,g) -> sumlessthanone(x), 1e-8)
max_objective!(opt, ((x, grad) -> myfunc(x)))
(minf,minx,ret) = optimize(opt, [0.2,0.3,0.1])
numevals = opt.numevals # the number of function evaluations
println("got $minf at $minx after $numevals iterations (returned $ret)")
using GaussianProcesses
using PDMats
using NLopt



function getbounds(fixGU::Bool, fixLambdaWeight::Bool, unpairedmodel::Bool, siteCats::Int=3, lambdaCats::Int=5,initialparams::Array{Float64,1}=2.0*ones(Float64,15); coevolutionparamsonly::Bool=false)
  lower = ones(Float64, 15)*1e-4
  upper = ones(Float64, 15)*50.0
  lower[1] = 1.0
  lower[2] = 1.0
  lower[3] = 1.0

  lower[4] = 0.02
  lower[5] = 0.02
  lower[6] = 0.02
  lower[7] = 0.02
  lower[8] = 0.02
  lower[9] = 0.02


  lower[12] = 0.05
  lower[13] = 0.05
  lower[14] = 0.05
  lower[15] = 0.0001  
  if unpairedmodel
      lower[10] = 1.0
      lower[15] = 0.9999
  end
  if lambdaCats == 1
    lower[10] = 1.0
    lower[15] = 1e-4
  end
  upper[12] = 0.5
  upper[13] = 0.5
  upper[14] = 0.5
  upper[15] = 0.9999  
  if unpairedmodel
      upper[1] = 1.0
      upper[2] = 1.0
      upper[3] = 1.0
      upper[10] = 1.0
  end
  if lambdaCats == 1
    upper[10] = 1.0
    upper[15] = 1e-4
  end
  if fixGU
    upper[3] = 1.0
  end
  if fixLambdaWeight
      upper[15] = 1.0
  end
  #=
  println("A")
  initialparams = 2.0*ones(Float64,15)
  if initparams == nothing
      initialparams[1] = 1.5
      initialparams[2] = 1.5
      initialparams[3] = 1.5
      initialparams[12] = 0.25
      initialparams[13] = 0.25
      initialparams[14] = 0.25
      initialparams[15] = 0.5
      println("B")
  else
      initialparams = getparamsvector(initparams)
  end=#

  for z=1:length(initialparams) # correct any out of bounds values
    initialparams[z] = max(lower[z], initialparams[z])
    initialparams[z] = min(upper[z], initialparams[z])
  end

  if coevolutionparamsonly
    for z=4:length(initialparams)
      lower[z] = initialparams[z]
      upper[z] = initialparams[z]
    end
  end

  return lower,upper,initialparams
end

function randomsample(gp::GPBase, x, r=randn(size(x,2),1), A::Array{Float64,2}=zeros(size(x,2),1))
	nobs = size(x,2)
    n_sample = 1

    if gp.nobs == 0
        # Prior mean and covariance
        μ = mean(gp.mean, x);
        Σraw = cov(gp.kernel, x, x);
        Σraw = Matrix(Σraw)
        Σraw, chol = GaussianProcesses.make_posdef!(Σraw)
        Σ = PDMat(Σraw, chol)
    else
        # Posterior mean and covariance
        μ, Σraw = predict_f(gp, x; full_cov=true)
        Σraw = Matrix(Σraw)
        Σraw, chol = GaussianProcesses.make_posdef!(Σraw)
        Σ = PDMat(Σraw, chol)
    end
    return broadcast!(+, A, μ, unwhiten(Σ, r))
end

function likelihood(gp::GPBase, xin, r)
	x = zeros(Float64, length(xin),1)
	for z=1:length(xin)
		x[z,1] = xin[z]
	end
	ll = randomsample(gp, x, r)[1]
	println(xin,"\t",ll,"\t",r)
	return ll
end

cachefile = "results/rhinovirus_a/rhinovirus_a.cache.list.sitecats3.lambdacats1"

x = Array{Float64,1}[]
y = Float64[]
fin = open(cachefile,"r")
for line in readlines(fin)
    spl = split(line,"\t")
    if length(spl) == 2
        push!(y, parse(Float64, spl[1]))
        arrstr = split(spl[2][2:end-1], ", ")
        push!(x, Float64[parse(Float64,v) for v in arrstr])
    end
end
x = x[end-200:end]
y = y[end-200:end]

dim = length(x[1])

model = ElasticGPE(dim,mean = MeanConst(0.),         
                   kernel = SEArd(zeros(Float64,dim), 5.),
                   capacity = 3000)   
append!(model, hcat(x...), y)

println("start")
timeelapsed = @elapsed GaussianProcesses.optimize!(model, noise=false)
println("end")
println(timeelapsed)
x = rand(dim,1)

lower,upper,initial = getbounds(false, false, false, 3, 1)

initial = lower .+ 0.5*(upper.-lower)

r = randn(1,1)
localObjectiveFunction = ((param, grad) -> likelihood(model,param,r))        
opt = Opt(:LN_NELDERMEAD, 15)
lower_bounds!(opt, lower)
upper_bounds!(opt, upper)
maxeval!(opt, 1000)
max_objective!(opt, localObjectiveFunction)  
xtol_rel!(opt,1e-4)
(minf,minx,ret) = NLopt.optimize(opt, initial)
println(minx,"\t",minf)
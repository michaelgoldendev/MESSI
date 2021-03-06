module CommonUtils
using Distributions
using SHA
using Random

export eye
function eye(n::Int)
    #Diagonal{Float64}(I, n)
    mat = zeros(Float64,n,n)
    for i=1:n
        mat[i,i] = 1.0
    end
    return mat
end

function spearmanspvalue(x1::Array{Float64,1}, y1::Array{Float64,1})
    n = length(x1)
    r = corspearman(x1,y1)
    t = r*sqrt((n-2.0)/(1.0-(r*r)))
    p = cdf(TDist(n-2), -abs(t))

    stderr = 1.0 /sqrt(n - 3.0)
    delta = 1.96 * stderr
    lower = tanh(atanh(r) - delta)
    upper = tanh(atanh(r) + delta)
    return r,lower,upper,n,p
end

function spearmanspvalue(paired::Array{Int,1}, x::Array{Float64,1}, y::Array{Float64,1})

    x1 = Float64[]
    y1 = Float64[]
    for i=1:length(x)
        if paired[i] > i
            if !isnan(x[i]) && !isnan(y[i])
                push!(x1, x[i] + randn()*1e-10)
                push!(y1, y[i] + randn()*1e-10)
            end
        end
    end
    n = length(x1)
    r = corspearman(x1,y1)
    t = r*sqrt((n-2.0)/(1.0-(r*r)))
    p = cdf(TDist(n-2), -abs(t))

    stderr = 1.0 /sqrt(n - 3.0)
    delta = 1.96 * stderr
    lower = tanh(atanh(r) - delta)
    upper = tanh(atanh(r) + delta)
    return r,lower,upper,n,p
end

function softmapping(mapping::Array{Int,1})
    ret = zeros(Int, length(mapping))
    for i=1:length(mapping)
        ret[i] = mapping[i]
        if mapping[i] <= 0
            for j=1:length(mapping)
                if i-j >= 1 && mapping[i-j] > 0
                    ret[i] = mapping[i-j]
                    break
                elseif i+j <= length(mapping) && mapping[i+j] > 0
                    ret[i] = mapping[i+j];
                    break
                end
            end
        end
    end
    return ret
end


function getrevmapping(mapping::Array{Int,1})
    ret = zeros(Int, maximum(mapping))
    for i=1:length(mapping)
        if mapping[i] > 0
            ret[mapping[i]] = i
        end
    end
    return ret
end


function sha256base36(input::AbstractString)
    return string(parse(BigInt, bytes2hex(sha256(input)), base=16), base=36)
end

function canonicalise_sequence(seq::AbstractString)
    return replace(uppercase(replace(strip(seq), "-" => "")), "U" => "T")
end

export sampletruncated
function sampletruncated(rng::AbstractRNG, d::UnivariateDistribution, lower::Float64, upper::Float64)
    v = quantile(d,rand(rng))
    while !(lower <= v <= upper)
        v = quantile(d,rand(rng))
    end
    return v
end

export combinelogarrays
function combinelogarrays(v1::Array{Float64,1}, v2::Array{Float64,1}, w::Float64)::Array{Float64,1}
    logw1 = log(w)
    logw2 = log(1.0-w)
    len = length(v1)
    res = zeros(Float64, len)
    for i=1:len
        res[i] = logsumexp(logw1+v1[i], logw2+v2[i])
    end
    return res
end

export combinelogarrays
function combinelogarrays(v1::Array{Float64,2}, v2::Array{Float64,2}, w::Float64)::Array{Float64,2}
    logw1 = log(w)
    logw2 = log(1.0-w)
    d1 = size(v1,1)
    d2 = size(v2,2)
    res = zeros(Float64, d1, d2)
    for i=1:d1
        for j=1:d2
            res[i,j] = logsumexp(logw1+v1[i,j], logw2+v2[i,j])
        end
    end
    return res
end

export sample
function sample(rng::AbstractRNG, v::Array{Float64})::Int
    s = sum(v)
    n = length(v)
    r = rand(rng)
    cumsum = 0.0
    for i=1:n
        cumsum += v[i] / s
        if cumsum >= r
            return i
        end
    end

    return n
end

export logsumexp
function logsumexp(a::Float64, b::Float64)::Float64
    if b == -Inf
        return a
    elseif a == -Inf
        return b
    elseif a < b
        return b + log1p(exp(a-b))
    else
        return a + log1p(exp(b-a))
    end
end

export logsumexp
function logsumexp(v::Array{Float64,1}, s::Int, e::Int)
  m = -1e20
  for i=s+1:e
    if v[i] > m
      m = v[i]
    end
  end
  cutoff = m - 30.0

  ll = v[s]
  a = 0.0
  b = 0.0
  diff = 0.0
  for i=s+1:e
    b = v[i]
    if b > cutoff
      a = ll
      diff = a-b
      if a == -Inf || diff < -30.0

      elseif b == -Inf || diff > 30.0
        ll = a
      elseif a < b
        ll = b + log1p(exp(a-b))
      else
        ll = a + log1p(exp(b-a))
      end
    end
  end
  return ll
end

export absmat
function absmat(M::Array{Float64,2})::Array{Float64,2}
    dim1 = size(M,1)
    dim2 = size(M,2)
    for i=1:dim1
        for j=1:dim2
            if M[i,j] <= 0.0
                M[i,j] = 1e-50
            end
        end
    end
    return M
end

export summarisemcmc
function summarisemcmc(mcmclogfile)
    fin = open(mcmclogfile,"r")
    params = Dict{AbstractString,Array{Float64,1}}()
    header = true
    keys = AbstractString[]
    for line in readlines(fin)
        spl = split(strip(line))
        if header
            keys = spl
            for key in keys
                params[key] = Float64[]
            end
            header = false
        else
            for i=1:length(spl)
                p = get(params, keys[i], Float64[])
                push!(p, parse(Float64,spl[i]))
                params[keys[i]] = p
            end
        end
    end
    close(fin)
    for key in keys
        len = length(params[key])
        if len > 1
            data = params[key][div(len,2):end]
            #println(key,"\t",length(data),"\t", mean(data), "\t", std(data))
        end
    end
    return keys,params
end

using HypothesisTests
export mannwhitneyu
function mannwhitneyu(x,y)
    nx = length(x)
    ny = length(y)
    bigN = nx+ny
    n = ((nx * ny) / (bigN * (bigN - 1)))
    tieCorrectionFactor = 0.0
    d = ((bigN^3.0 - bigN) / 12.0 - tieCorrectionFactor)
    variance = sqrt(n*d)
    mann = MannWhitneyUTest(x, y)
    zscore =  (mann.U - (nx * ny / 2.0)) / variance
    return zscore, pvalue(mann)
end
end

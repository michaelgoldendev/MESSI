include("RNATools.jl")

#include("FancyPhylogenetics.jl")
#using FancyPhylogenetics
#include("../../juliamolev/src/MCMCLogging.jl")
#include("/mnt/usb-Seagate_BUP_Slim_BK_NA7QWPMK-0:0-part1/juliamolev/src/MCMCLogging.jl")
using MCMCLogging
#include("/mnt/usb-Seagate_BUP_Slim_BK_NA7QWPMK-0:0-part1/juliamolev/src/MolecularEvolution.jl")
#include("../../juliamolev/src/MolecularEvolution.jl")
using MolecularEvolution
#using CUDArt
using FastaIO
using Distributions
include("TimingsLogger.jl")
using LinearAlgebra
using CommonUtils
using Random
using SharedArrays

mutable struct FastCache
  size::Int
  linewidth::Int
  keys::Array{Int,1}
  cache::Array{Float64,2}
  logm::Array{Float64,1}

  function FastCache(size::Int, linewidth::Int)
    keys = zeros(Int,size)
    cache = zeros(Float64, size, linewidth)
    logm = zeros(Float64, size)
    new(size,linewidth,keys,cache, logm)
  end
end

function putarray(fastcache::FastCache, key::Int, v::Array{Float64,1})
  index = hash(key) % fastcache.size + 1
  fastcache.keys[index] = key
  for i=1:fastcache.linewidth
    fastcache.cache[index,i] = v[i]
  end
  return index
end

function getcacheindex(fastcache::FastCache, key::Int)
  index = hash(key) % fastcache.size + 1
  if fastcache.keys[index] == key
    return index
  end
  return 0
end

function getandkillcacheindex(fastcache::FastCache, key::Int)
  index = hash(key) % fastcache.size + 1
  fastcache.keys[index] = 0
  return index
end

#=
function getcolumnkey(paramindex::Int, numparams::Int, subcolumnrefs::Array{Int,2}, numcols::Int, columns::Array{Int,1})
  nodeindex = 1
  len = length(columns)
  key = 0
  for i=1:len
    key *= numcols
    key += (subcolumnrefs[nodeindex,columns[i]]-1)
  end
  key *= numparams
  key += paramindex
  return key+1
end=#
function getcolumnkey(subcolumnrefs::Array{Int,2}, numcols::Int, columns::Array{Int,1})
  nodeindex = 1
  len = length(columns)
  key = 0
  for i=1:len
    key *= numcols
    key += (subcolumnrefs[nodeindex,columns[i]]-1)
  end
  return key+1
end

mutable struct ColumnCache
  size::Int
  keys::Array{Int,1}
  cache::Array{Float64,1}

  function ColumnCache(size::Int)
    keys = zeros(Int,size)
    cache = zeros(Float64, size)
    new(size,keys,cache)
  end
end

function putvalue(cache::ColumnCache, key::Int, value::Float64)
  index = hash(key) % cache.size + 1
  cache.keys[index] = key
  cache.cache[index] = value
  return index
end

function getcacheindex(cache::ColumnCache, key::Int)
  index = hash(key) % cache.size + 1
  if cache.keys[index] == key
    return index
  end
  return 0
end



mutable struct PairedProbs
  probs::Array{Float64,1}

  PairedProbs(probs::Array{Float64,1}) = new(probs)
end

function getpairedsitestring(paired)
  pairedsites = ""
  unpairedsites = ""
  numpaired = 0
  numunpaired = 0
  for x=1:length(paired)
    y = paired[x]
    if y > x
      if length(pairedsites) > 0
        pairedsites = string(pairedsites,",")
      end
      pairedsites = string(pairedsites,x-1,",",y-1)
      numpaired += 1
    elseif y == 0
      if length(unpairedsites) > 0
        unpairedsites = string(unpairedsites,",")
      end
      unpairedsites = string(unpairedsites,x-1)
      numunpaired += 1
    end
  end
  return pairedsites, unpairedsites, numpaired, numunpaired
end


mutable struct PartialMatrix
  len::Int
  dlen::Int
  matrix::Array{Float64,2}
  defaultvalue::Float64

  function PartialMatrix(len::Int, dlen::Int, defaultvalue::Float64=-Inf)
    matrix = ones(Float64,len,dlen+1)*defaultvalue
    new(len, dlen+1, matrix, defaultvalue)
  end
end

function Base.getindex(A::PartialMatrix, i1::Real, i2::Real)
  j = i2-i1
  if abs(j) > A.dlen || i1 < 1 || i1 > A.len
    return A.defaultvalue
  end
  if j > 0
    return A.matrix[i1,j]
  elseif j < 0
    return A.matrix[i2,-j]
  end
  return A.defaultvalue
end

function Base.setindex!(A::PartialMatrix, x::Real, i1::Real, i2::Real)
  j = i2-i1
  if j > 0
    A.matrix[i1,j] = x
  elseif j < 0
    A.matrix[i2,-j] = x
  end
end

function gibbsampler(numGridPoints::Int, numSites::Int, conditionals::Array{Float64,2}, lambdas::Array{Float64,1}, alpha::Float64, iters::Int=1000)
  rng = MersenneTwister(719294791640174294)
  φ = zeros(Float64, numGridPoints)
  θ = ones(Float64, numGridPoints)/numGridPoints
  v = zeros(Float64, numGridPoints)
  θsum = zeros(Float64, numGridPoints)
  burnin::Int = div(iters, 5)

  gibbswriter = open("farce.log", "w")
  c::Float64 = 0.0
  s::Float64 = 0.0
  write(gibbswriter, string("iter\t",join([string("site",i) for i=1:numSites],"\t"), "\n"))

  tempalloc = zeros(Int, numSites)
  postprobs = zeros(Float64, numSites)
  totaliters = 0.0
  for iter=1:iters
    @simd for j=1:numGridPoints
      @inbounds φ[j] = alpha
    end

    for i=1:numSites
        s = 0.0
        @simd for j=1:numGridPoints
          @inbounds v[j] = θ[j]*conditionals[j,i]
          @inbounds s += v[j]
        end
        siteallocation = CommonUtils.sample(rng,v,s)
        tempalloc[i] = siteallocation
        @inbounds φ[siteallocation] += 1.0
    end
    φ[div(length(lambdas),2)+1] += numSites
    #println(lambdas[div(length(lambdas),2)+1])

    if iter % 20 == 1
      write(gibbswriter, string(iter-1,"\t",join([string(lambdas[tempalloc[i]]) for i=1:numSites],"\t"), "\n"))
    end


    rand!(Dirichlet(φ), θ)
    if iter > burnin
      for i=1:numSites
        if lambdas[tempalloc[i]] > 1.0
          postprobs[i] += 1
        end
      end
      totaliters += 1.0
      @simd for j=1:numGridPoints
        @inbounds θsum[j] += θ[j]
      end
      c += 1.0
    end
  end

  println(postprobs/totaliters)

  close(gibbswriter)
  return θsum/c
end

#include("InsideOutsideAlgorithm.jl")

function getpartitions(len::Int,numparts::Int)
   	numpartitions = numparts
   	numpairedsites = div(len*len-len,2)
   	partsize = div(numpairedsites,numpartitions)
   	ret = Any[]
   	spos = 1
   	tmp = 0
   	for i=1:len
   		tmp += len-i
   		if tmp > partsize || i == len
   			push!(ret, (spos,i,tmp))
   			spos = i+1
   			numpairedsites -= tmp
   			numpartitions -= 1
   			if numpartitions > 0
   				partsize = div(numpairedsites,numpartitions)
   			end
   			tmp = 0
   		end
   	end

   	return ret
end

function summarizemcmc(mcmclogfile)
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
      println(key,"\t",length(data),"\t", mean(data), "\t", std(data))
    end
  end
  len = length(params[keys[1]])
  println("lambdaGC > lambdaAT: ", mean([lambdaGC > lambdaAT ? 1.0 : 0.0 for (lambdaAT,lambdaGC) in zip(params["lambdaAT"][div(len,2):end],params["lambdaGC"][div(len,2):end])]))
  println("lambdaAT > lambdaGT: ", mean([lambdaAT > lambdaGT ? 1.0 : 0.0 for (lambdaGT,lambdaAT) in zip(params["lambdaGT"][div(len,2):end],params["lambdaAT"][div(len,2):end])]))
  println("lambdaGC > lambdaGT: ", mean([lambdaGC > lambdaGT ? 1.0 : 0.0 for (lambdaGT,lambdaGC) in zip(params["lambdaGT"][div(len,2):end],params["lambdaGC"][div(len,2):end])]))

  startindex = div(len,2)
  len2 = length(params["lambdaGC"][startindex:end])
  data = zeros(len2,3)
  for i=1:len2
    data[i,1] = params["lambdaGC"][startindex+i-1]
    data[i,2] = params["lambdaAT"][startindex+i-1]
    data[i,3] = params["lambdaGT"][startindex+i-1]
  end
  A = full(chol(cov(data)))
  return keys,params,A
end

function getmatrixstring(A)
  rows = AbstractString[]
  for i=1:size(A,1)
    push!(rows,string("{",join(AbstractString[string(v) for v in A[i,:]],","),"}"))
  end
  return string("{",join(rows,","),"}")
end

function getdata(fastafile)
  sequences = AbstractString[]
  names = AbstractString[]
  seqnametoindex = Dict{AbstractString,Int}()
  len = 0
  seqindex = 1
  FastaReader(fastafile) do fr
   for (desc, seq) in fr
     len = length(seq)
     push!(names,desc)
     push!(sequences, seq)
     seqnametoindex[desc] = seqindex
     seqindex += 1
   end
  end
  data = zeros(Float64,length(sequences),len,4)
  for s=1:length(sequences)
    seq = sequences[s]
    for j=1:len
      nuc = get(nucmapping,seq[j],0)
      if nuc > 0
        data[s,j,nuc] = 1.0
      else
        data[s,j,:] = 1.0
      end
    end
  end

  return data, seqnametoindex
end

function state(data::Array{Float64,3}, seq::Int, columns::Array{Int,1})
  if length(columns) == 2
    dinuc = zeros(Float64,16)
    for i=1:4
      for j=1:4
        dinuc[(i-1)*4+j] = data[seq,columns[1],i]*data[seq,columns[2],j]
      end
    end
    return dinuc
  end
end

kronsum(a::Vector, b::Vector) = vec([ a[i]+b[j] for j=1:length(b), i=1:length(a)])
kronprod(a::Vector, b::Vector) = vec([ a[i]*b[j] for j=1:length(b), i=1:length(a)])

function logstate(logdata::Array{Float64,3}, seq::Int, columns::Array{Int,1})
  v = vec(logdata[seq,columns[1],:])
  for c=2:length(columns)
    v = kronsum(v, vec(logdata[seq,columns[c],:]))
  end
  return v
end

function state(data::Array{Float64,3}, seq::Int, columns::Array{Int,1})
  v = vec(data[seq,columns[1],:])
  for c=2:length(columns)
    v = kronprod(v, vec(data[seq,columns[c],:]))
  end
  push!(v,0.0)
  return v
end

function getkey(numnodes::Int, nodeindex::Int, subcolumnrefs::Array{Int,2}, numcols::Int, columns::Array{Int,1})
  len = length(columns)
  key = (nodeindex-1)
  for i=1:len
    key *= numcols
    key += (subcolumnrefs[nodeindex,columns[i]]-1)
  end
  return key+1
end

function absmat(M::Array{Float64,2}, Q::Array{Float64,2}, branchlength::Float64)
  print = false
  dim1 = size(M,1)
  dim2 = size(M,2)
  for i=1:dim1
    for j=1:dim2
      if M[i,j] <= 0.0
        print = true
        M[i,j] = 1e-50
      end
    end
  end
  print = false
  if print
    println("branchlength",branchlength)
    println(M)
    println(Q)
  end
  return M
end

function gettransitionmatriceslist(branchlengths::Array{Float64,1}, Q::Array{Float64,2}, takelog::Bool=false)
  t = 1.0
  decomposition = eigen(Q)
  D, V = decomposition.values, decomposition.vectors
  Vi = inv(V)
  dim = size(Q,1)
  Plist = Array{Float64,2}[zeros(Float64,dim,dim)]
  for n=2:length(branchlengths)
    if takelog
      if branchlengths[n] == 0.0
        push!(Plist, log(eye(dim)))
      else
        push!(Plist, log(V*Diagonal(exp.(D*branchlengths[n]))*Vi))
      end
    else
      if branchlengths[n] == 0.0
        push!(Plist, eye(dim))
      else
        push!(Plist, absmat(real(V*Diagonal(exp.(D*branchlengths[n]))*Vi), Q, branchlengths[n]))
      end
    end
  end
  return Plist
end

mutable struct MuseModel
  alphabet::Int
  obsfreqs::Array{Float64,1}
  freqs::Array{Float64,1}
  logfreqs::Array{Float64,1}
  Q::Array{Float64,2}
  lambdaGC::Float64
  lambdaAT::Float64
  lambdaGT::Float64
  q1::Float64
  q2::Float64
  q3::Float64
  q4::Float64
  q5::Float64
  q6::Float64
  siteRate1::Float64
  siteRate2::Float64

  function MuseModel(obsfreqs::Array{Float64,1}, lambdaGC::Float64, lambdaAT::Float64, lambdaGT::Float64, q1::Float64, q2::Float64, q3::Float64, q4::Float64, q5::Float64, q6::Float64, siteRate1::Float64, siteRate2::Float64)
    freqs = calculatemusefreqs(obsfreqs,lambdaGC,lambdaAT,lambdaGT)
    logfreqs = log.(freqs)
    Q =  muse(lambdaGC, lambdaAT, lambdaGT, q1, q2, q3, q4, q5, q6, obsfreqs, siteRate1, siteRate2)
    new(16, obsfreqs, freqs, logfreqs, Q, lambdaGC, lambdaAT, lambdaGT, q1, q2, q3, q4, q5, q6, siteRate1, siteRate2)
  end
end

mutable struct ModelParameters
  freqs::Array{Float64,1}
  lambdaGC::Float64
  lambdaAT::Float64
  lambdaGT::Float64
  q1::Float64
  q2::Float64
  q3::Float64
  q4::Float64
  q5::Float64
  q6::Float64
  branchlengths::Array{Float64,1}
  lambdazeroweight::Float64
  lambdaratesGC::Array{Float64,1}
  lambdaweightsGC::Array{Float64,1}
  lambdaratesAT::Array{Float64,1}
  lambdaweightsAT::Array{Float64,1}
  lambdaratesGT::Array{Float64,1}
  lambdaweightsGT::Array{Float64,1}
  lambdaCats::Int
  lambdaGammaShapeGC::Float64
  lambdaGammaScaleGC::Float64
  lambdaGammaShapeAT::Float64
  lambdaGammaScaleAT::Float64
  lambdaGammaShapeGT::Float64
  lambdaGammaScaleGT::Float64
  pairedstates::Array{Int,1}
  siteRates::Array{Float64,1}
  siteWeights::Array{Float64,1}
  siteCats::Int
  siteGammaShape::Float64
  siteGammaScale::Float64
  states::Array{Int,1}
  logposterior::Float64
  loglikelihood::Float64
  logprior::Float64


  function ModelParameters(freqs::Array{Float64,1}, lambdaGC::Float64, lambdaAT::Float64, lambdaGT::Float64, q1::Float64, q2::Float64, q3::Float64, q4::Float64, q5::Float64, q6::Float64, nodelist::Array{TreeNode,1}, lambdazeroweight::Float64=0.5, lambdaratesGC::Array{Float64,1}=ones(Float64,1), lambdaweightsGC::Array{Float64,1}=ones(Float64,1), lambdaratesAT::Array{Float64,1}=ones(Float64,1), lambdaweightsAT::Array{Float64,1}=ones(Float64,1), lambdaratesGT::Array{Float64,1}=ones(Float64,1), lambdaweightsGT::Array{Float64,1}=ones(Float64,1), lambdaCats::Int=1, lambdaGammaScaleGC::Float64=1.0, lambdaGammaShapeGC::Float64=1.0, lambdaGammaScaleAT::Float64=1.0, lambdaGammaShapeAT::Float64=1.0, lambdaGammaScaleGT::Float64=1.0, lambdaGammaShapeGT::Float64=1.0, pairedstates::Array{Int,1}=ones(Int,1), siteWeights::Array{Float64,1}=ones(Float64,1), siteRates::Array{Float64,1}=ones(Float64,1), siteCats::Int=1, siteGammaShape::Float64=1.0, siteGammaScale::Float64=1.0, states::Array{Int,1}=ones(Int,1))
    branchlengths = [node.branchlength for node in nodelist]
    new(freqs, lambdaGC, lambdaAT, lambdaGT, q1, q2, q3, q4, q5, q6, branchlengths, lambdazeroweight, lambdaratesGC, lambdaweightsGC, lambdaratesAT, lambdaweightsAT, lambdaratesGT, lambdaweightsGT, lambdaCats, lambdaGammaShapeGC, lambdaGammaScaleGC, lambdaGammaShapeAT, lambdaGammaScaleAT, lambdaGammaShapeGT, lambdaGammaScaleGT, pairedstates, siteWeights, siteRates, siteCats, siteGammaShape, siteGammaScale, states,0.0,0.0,0.0)
  end
end

function getbounds(fixGU::Bool, fixLambdaWeight::Bool, unpairedmodel::Bool, siteCats::Int=3, lambdaCats::Int=5,initialparams::Array{Float64,1}=2.0*ones(Float64,15))
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
  upper[12] = 0.80
  upper[13] = 0.80
  upper[14] = 0.80
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

  return lower,upper,initialparams
end

function getparamsvector(params::ModelParameters)
    p = zeros(Float64,15)
    p[1] = params.lambdaGC
    p[2] = params.lambdaAT
    p[3] = params.lambdaGT
    p[4] = params.q1
    p[5] = params.q2
    p[6] = params.q3
    p[7] = params.q4
    p[8] = params.q5
    p[9] = params.q6
    p[10] = params.lambdaGammaShapeGC
    p[11] = params.siteGammaShape
    """
    p[12] = params.freqs[1]
    p[13] = params.freqs[2]/(1.0-params.freqs[1])
    p[14] = params.freqs[3]/(1.0-params.freqs[1]-params.freqs[2])
    """
    p[12] = params.freqs[1]
    p[13] = params.freqs[2]
    p[14] = params.freqs[3]

    p[15] = params.lambdazeroweight
    return p
end

function getparams(params::Array{Float64,1}, dataset::Dataset, siteCats::Int=3, lambdacats::Int=5, parameterisation::Int=1,fixGU::Bool=false,fixGCAU::Bool=false)
  #=
  f1 = params[12]
  f2 = (1.0-f1)*params[13]
  f3 = (1.0-f1-f2)*params[14]
  f4 = (1.0-f1-f2-f3)
  freqs = Float64[f1, f2, f3, f4]=#
  freqs = Float64[params[12], params[13], params[14], 1.0 - params[12] - params[13] - params[14]]
  currentparams = ModelParameters(freqs, params[1], params[2], params[3], params[4], params[5], params[6], params[7], params[8], params[9], getnodelist(dataset.root), params[15])
  currentparams.lambdaGC = params[1]
  currentparams.lambdaAT = params[2]
  currentparams.lambdaGT = params[3]
  currentparams.lambdaGammaShapeGC = params[10]
  currentparams.siteGammaShape = params[11]
  currentparams.siteGammaScale = 1.0 / currentparams.siteGammaShape
  currentparams.lambdaCats = lambdacats - 1

  #currentparams.lambdaweightsGC, currentparams.lambdaratesGC = discretizegamma2(currentparams.lambdazeroweight, currentparams.lambdaGammaShapeGC, currentparams.lambdaGammaScaleGC, currentparams.lambdaCats,parameterisation)
  #currentparams.lambdaweightsAT, currentparams.lambdaratesAT = discretizegamma2(currentparams.lambdazeroweight, currentparams.lambdaGammaShapeAT, currentparams.lambdaGammaScaleAT, currentparams.lambdaCats,parameterisation)
  #currentparams.lambdaweightsGT, currentparams.lambdaratesGT = discretizegamma2(currentparams.lambdazeroweight, currentparams.lambdaGammaShapeGT, currentparams.lambdaGammaScaleGT, currentparams.lambdaCats,parameterisation)
  currentparams.lambdaweightsGC, currentparams.lambdaratesGC, currentparams.lambdaweightsAT, currentparams.lambdaratesAT, currentparams.lambdaweightsGT, currentparams.lambdaratesGT = discretizegamma3(currentparams.lambdazeroweight, currentparams.lambdaGC, currentparams.lambdaAT, currentparams.lambdaGT, currentparams.lambdaGammaScaleGC, currentparams.lambdaCats, parameterisation)

  if fixGCAU
    currentparams.lambdaratesGC = copy(currentparams.lambdaratesAT)
  end
  if fixGU
    currentparams.lambdaratesGT = zeros(Float64, length(currentparams.lambdaratesGT))
  end
  #println(currentparams.lambdaweightsGC,"\t",currentparams.lambdaweightsAT, "\t", currentparams.lambdaweightsGT)
  #println(currentparams.lambdaratesGC,"\t",currentparams.lambdaratesAT, "\t", currentparams.lambdaratesGT)
  #println(currentparams.lambdaGammaShapeGC,"\t",currentparams.lambdaGammaScaleGC, "\t", currentparams.lambdaGammaShapeAT, "\t", currentparams.lambdaGammaScaleAT,"\t",currentparams.lambdaGammaShapeGT, "\t", currentparams.lambdaGammaScaleGT)
  #println(sum((currentparams.lambdaratesGC+1.0).*currentparams.lambdaweightsGC),"\t",sum((currentparams.lambdaratesAT+1.0).*currentparams.lambdaweightsAT), "\t", sum((currentparams.lambdaratesGT+1.0).*currentparams.lambdaweightsGT))
  currentparams.siteCats = siteCats
  currentparams.siteWeights = ones(Float64,currentparams.siteCats)/currentparams.siteCats
  currentparams.siteRates = discretizegamma(currentparams.siteGammaShape, currentparams.siteGammaScale, currentparams.siteCats)
  currentparams.states = Int[rand(1:siteCats) for i=1:dataset.numcols]
  if length(currentparams.lambdaweightsGC) == 1
    currentparams.lambdaratesGC[1] = currentparams.lambdaGC-1.0
    currentparams.lambdaratesAT[1] = currentparams.lambdaAT-1.0
    currentparams.lambdaratesGT[1] = currentparams.lambdaGT-1.0
    currentparams.lambdaweightsGC[1] = 1.0
    currentparams.lambdaweightsAT[1] = 1.0
    currentparams.lambdaweightsGT[1] = 1.0
  end
  #println(currentparams.lambdaratesGC,"\t",currentparams.lambdaweightsGC)
  #println(currentparams.lambdaratesAT,"\t",currentparams.lambdaweightsAT)
  #println(currentparams.lambdaratesGT,"\t",currentparams.lambdaweightsGT)
  return currentparams
end

#=
function getparams(params::Array{Float64,1}, dataset::Dataset, siteCats::Int=3, lambdacats::Int=5, parameterisation::Int=1,fixGU::Bool=false,fixGCAU::Bool=false)
  freqs = Float64[params[12],params[13],params[14], 1.0 - params[12] - params[13] - params[14]]
  currentparams = ModelParameters(freqs, params[1], params[2], params[3], params[4], params[5], params[6], params[7], params[8], params[9], getnodelist(dataset.root), params[15])
  currentparams.lambdaGammaShapeGC = params[10]
  currentparams.lambdaGammaScaleGC = params[16]
  currentparams.lambdaGammaShapeAT = params[17]
  currentparams.lambdaGammaScaleAT = params[18]
  currentparams.lambdaGammaShapeGT = params[19]
  currentparams.lambdaGammaScaleGT = params[20]
  currentparams.siteGammaShape = params[11]
  currentparams.siteGammaScale = 1.0 / currentparams.siteGammaShape
  currentparams.lambdaCats = lambdacats - 1

  #currentparams.lambdaweightsGC, currentparams.lambdaratesGC = discretizegamma2(currentparams.lambdazeroweight, currentparams.lambdaGammaShapeGC, currentparams.lambdaGammaScaleGC, currentparams.lambdaCats,parameterisation)
  #currentparams.lambdaweightsAT, currentparams.lambdaratesAT = discretizegamma2(currentparams.lambdazeroweight, currentparams.lambdaGammaShapeAT, currentparams.lambdaGammaScaleAT, currentparams.lambdaCats,parameterisation)
  #currentparams.lambdaweightsGT, currentparams.lambdaratesGT = discretizegamma2(currentparams.lambdazeroweight, currentparams.lambdaGammaShapeGT, currentparams.lambdaGammaScaleGT, currentparams.lambdaCats,parameterisation)
  currentparams.lambdaweightsGC, currentparams.lambdaratesGC, currentparams.lambdaweightsAT, currentparams.lambdaratesAT, currentparams.lambdaweightsGT, currentparams.lambdaratesGT = discretizegamma3(currentparams.lambdazeroweight, currentparams.lambdaGammaShapeGC, currentparams.lambdaGammaShapeAT, currentparams.lambdaGammaShapeGT, currentparams.lambdaGammaScaleGC, currentparams.lambdaCats, parameterisation)

  if fixGCAU
    currentparams.lambdaratesGC = copy(currentparams.lambdaratesAT)
  end
  if fixGU
    currentparams.lambdaratesGT = zeros(Float64, length(currentparams.lambdaratesGT))
  end
  #println(currentparams.lambdaweightsGC,"\t",currentparams.lambdaweightsAT, "\t", currentparams.lambdaweightsGT)
  #println(currentparams.lambdaratesGC,"\t",currentparams.lambdaratesAT, "\t", currentparams.lambdaratesGT)
  #println(currentparams.lambdaGammaShapeGC,"\t",currentparams.lambdaGammaScaleGC, "\t", currentparams.lambdaGammaShapeAT, "\t", currentparams.lambdaGammaScaleAT,"\t",currentparams.lambdaGammaShapeGT, "\t", currentparams.lambdaGammaScaleGT)
  #println(sum((currentparams.lambdaratesGC+1.0).*currentparams.lambdaweightsGC),"\t",sum((currentparams.lambdaratesAT+1.0).*currentparams.lambdaweightsAT), "\t", sum((currentparams.lambdaratesGT+1.0).*currentparams.lambdaweightsGT))
  currentparams.siteCats = siteCats
  currentparams.siteWeights = ones(Float64,currentparams.siteCats)/currentparams.siteCats
  currentparams.siteRates = discretizegamma(currentparams.siteGammaShape, currentparams.siteGammaScale, currentparams.siteCats)
  currentparams.states = Int[rand(1:siteCats) for i=1:dataset.numcols]

  return currentparams
end=#




mutable struct MuseSpecificParameters
  lambdacat::Int
  prob::Float64
  logprob::Float64
  lambdarate::Float64
  lambdaGC::Float64
  lambdaAT::Float64
  lambdaGT::Float64

  function MuseSpecificParameters(lambdacat::Int, prob::Float64, lambdarate::Float64, lambdaGC::Float64, lambdaAT::Float64, lambdaGT::Float64)
    return new(lambdacat,prob,log(prob),lambdarate,lambdaGC,lambdaAT,lambdaGT)
  end
end

function getmusespecificparamsarray(params::ModelParameters)
  museparams = MuseSpecificParameters[]
  lambdacat = 1
  for (lambdaweight,lambdarate) in zip(params.lambdaweightsGC,params.lambdaratesGC)
    prob = lambdaweight
    #println("lambdaweight", lambdaweight)
    musespecificparam = MuseSpecificParameters(lambdacat, prob, lambdarate, params.lambdaGC, params.lambdaAT, params.lambdaGT)
    push!(museparams, musespecificparam)
    lambdacat += 1
  end
  return museparams
end

function getGC(params::ModelParameters, musespecificparams::MuseSpecificParameters)
  #println("getGC ", params.lambdaratesGC[musespecificparams.lambdacat]+1.0)
  return params.lambdaratesGC[musespecificparams.lambdacat]+1.0
end

function getAT(params::ModelParameters, musespecificparams::MuseSpecificParameters)
  #println("getAT ", params.lambdaratesAT[musespecificparams.lambdacat]+1.0)
  return params.lambdaratesAT[musespecificparams.lambdacat]+1.0
end

function getGT(params::ModelParameters, musespecificparams::MuseSpecificParameters)
  #println("getGT ", params.lambdaratesGT[musespecificparams.lambdacat]+1.0)
  return params.lambdaratesGT[musespecificparams.lambdacat]+1.0
end


function computeunpairedlikelihoodscat(dataset::Dataset, params::ModelParameters, states::Array{Int,1}, sites::Array{Int,1}, rateCat::Int, fastcache::FastCache=FastCache(100000,4))
  alphabet = 4
  root = dataset.root
  nodelist = getnodelist(root)
  revnodelist = reverse(nodelist)
  siteloglikelihoods = ones(Float64, dataset.numcols)*-Inf
  rate = params.siteRates[rateCat]
  freqs = params.freqs
  Q = gtr(params.q1, params.q2, params.q3, params.q4, params.q5, params.q6, freqs)
  transprobs = gettransitionmatriceslist(params.branchlengths, Q*rate)

  cache = Array{Float64,1}[zeros(Float64,alphabet) for i=1:dataset.numnodes]
  subcolumnrefs = dataset.subcolumnrefs
  loglikcache = Dict{Tuple,Float64}()
  logm = zeros(Float64, dataset.numnodes)
  for col in sites
    if states[col] == rateCat
      for z=1:dataset.numnodes
        cache[z][1] = -Inf
      end
      fill!(logm, 0.0)
      v = felsensteinstack(dataset, nodelist, cache, logm, fastcache, subcolumnrefs, transprobs, root, dataset.data, Int[col], alphabet)
      total = 0.0
      for a=1:alphabet
        total += freqs[a] * v[a]
      end
      total = log(total) + logm[1]
      siteloglikelihoods[col] = total
    end
  end
  return siteloglikelihoods
end

function computeunpairedlikelihoodscats(dataset::Dataset, params::ModelParameters, states::Array{Int,1}, sites::Array{Int,1})
  refs = []
  for rateCat=1:params.siteCats
    #ref = @spawn computeunpairedlikelihoodscat(dataset, params, states, sites, rateCat)
    ref = computeunpairedlikelihoodscat(dataset, params, states, sites, rateCat)
    push!(refs,ref)
  end

  siteloglikelihoods = ones(Float64, dataset.numcols)*-Inf
  for rateCat=1:params.siteCats
    #temp = @fetch(refs[rateCat])
    temp = refs[rateCat]
    for col=1:dataset.numcols
      if temp[col] != -Inf
        siteloglikelihoods[col] = temp[col]
      end
    end
  end
  return siteloglikelihoods
end


function computedpairedlikelihoods(dataset::Dataset, params::ModelParameters, paired::Array{Int,1}, parallel::Bool=true)
  if parallel
    refs = []
    ll = -Inf
    for lambdacat=1:length(params.lambdarates)
      lambdarate = params.lambdarates[lambdacat]
      #ref = @spawn computepairedll(dataset, params, paired, lambdarate, lambdacat, 1.0, 1.0)
      ref = computepairedll(dataset, params, paired, lambdarate, lambdacat, 1.0, 1.0)
      push!(refs, ref)
    end
    likelihoods = Float64[]
    for i=1:length(refs)
      #temp = fetch(refs[i])
      temp = refs[i]
      loglambdaweights = log(params.lambdaweightsGC[i])
      if i == 1
        likelihoods = temp + loglambdaweights
      else
        for j=1:length(temp)
          likelihoods[j] = CommonUtils.logsumexp(likelihoods[j], loglambdaweights+temp[j])
        end
      end
    end
    return sum(likelihoods)
  else
    likelihoods = Float64[]
    for i=1:length(params.lambdarates)
      lambdarate = params.lambdarates[i]
      temp = computepairedll(dataset, params, paired, lambdarate, i)
      loglambdaweights = log(params.lambdaweightsGC[i])
      if i == 1
        likelihoods = loglambdaweights+temp
      else
        for j=1:length(temp)
          likelihoods[j] = CommonUtils.logsumexp(likelihoods[j], loglambdaweights+temp[j])
        end
      end
    end
    return sum(likelihoods)
  end
end

function samplelambdastates(rng::AbstractRNG, dataset::Dataset, params::ModelParameters, states::Array{Int,1}, paired::Array{Int,1}, parallel::Bool=true)
  refs = []
  for siteCat1=1:params.siteCats
    for siteCat2=1:params.siteCats
      maskedpaired = copy(paired)
      for i=1:length(maskedpaired)
        if maskedpaired[i] > i
          if states[i] == siteCat1 && states[maskedpaired[i]] == siteCat2

          else
            maskedpaired[maskedpaired[i]] = 0
            maskedpaired[i] = 0
          end
        end
      end
      for lambdacat=1:length(params.lambdarates)
        lambdarate = params.lambdarates[lambdacat]
        siteRate1 = params.siteRates[siteCat1]
        siteRate2 = params.siteRates[siteCat2]
        ref = computepairedll(dataset, params, maskedpaired, lambdarate, lambdacat, siteRate1, siteRate2, siteCat1, siteCat2)
        push!(refs, ref)
      end
    end
  end

  loglikelihoodcats = zeros(Float64,dataset.numcols, params.lambdaCats)
  loglikelihoods = zeros(Float64,dataset.numcols)
  index = 1
  for j=1:params.siteCats
    for k=1:params.siteCats
      for lambdacat=1:length(params.lambdarates)
        loglambdaweight = log(params.lambdaweightsGC[lambdacat])
        #temp = fetch(refs[index])
        temp = refs[index]
        pairindex = 1
        for i=1:dataset.numcols
          #sample = zeros()
          if paired[i] > i && states[i] == j && states[paired[i]] == k
            loglikelihoodcats[i, lambdacat] = loglambdaweight+temp[pairindex]
            if lambdacat == 1
              loglikelihoods[i] = loglambdaweight+temp[pairindex]
            else
              loglikelihoods[i] = CommonUtils.logsumexp(loglikelihoods[i], loglambdaweight+temp[pairindex])
            end
            pairindex += 1
          end
        end
        index += 1
      end
    end
  end

  sampledlambdastates = zeros(Int, dataset.numcols)
  for i=1:dataset.numcols
    if paired[i] > i
      loglikelihoodcats[i,:] = exp(loglikelihoodcats[i,:]-loglikelihoods[i])
      sampledlambdastates[i] = CommonUtils.sample(rng, loglikelihoodcats[i,:])
    end
  end

  return sampledlambdastates
end

function computepairedllhelper(dataset::Dataset, params::ModelParameters, states::Array{Int,1}, paired::Array{Int,1}, lambdarate::Float64, lambdacat::Int, siteRate1::Float64, siteRate2::Float64, siteCat1::Int, siteCat2::Int)
  maskedpaired = copy(paired)
  for i=1:length(maskedpaired)
    if maskedpaired[i] > i
      if states[i] == siteCat1 && states[maskedpaired[i]] == siteCat2

      else
        maskedpaired[maskedpaired[i]] = 0
        maskedpaired[i] = 0
      end
    end
  end
  return computepairedll(dataset, params, maskedpaired, lambdarate, lambdacat, siteRate1, siteRate2, siteCat1, siteCat2)
end

function computedpairedlikelihoodscatshelper(dataset::Dataset, params::ModelParameters, states::Array{Int,1}, paired::Array{Int,1}, siteCat1::Int, siteCat2::Int, keeplikelihoods::Bool=false)
  loglikelihoods = zeros(Float64, dataset.numcols)
  siteRate1 = params.siteRates[siteCat1]
  siteRate2 = params.siteRates[siteCat2]
  maskedpaired = copy(paired)
  for i=1:length(maskedpaired)
    if maskedpaired[i] > i
      if states[i] == siteCat1 && states[maskedpaired[i]] == siteCat2

      else
        maskedpaired[maskedpaired[i]] = 0
        maskedpaired[i] = 0
      end
    end
  end

  #=
  refs = []
  for lambdacat=1:length(params.lambdarates)
    lambdarate = params.lambdarates[lambdacat]
    #temp = computepairedll(dataset, params, maskedpaired, lambdarate, lambdacat, siteRate1, siteRate2, siteCat1, siteCat2)
    ref = @spawn computepairedll(dataset, params, maskedpaired, lambdarate, lambdacat, siteRate1, siteRate2, siteCat1, siteCat2)
    push!(refs,ref)
  end=#

  ret = Array{Float64,1}[]
  for lambdacat=1:length(params.lambdarates)
    loglambdaweight = log(params.lambdaweightsGC[lambdacat])
    lambdarate = params.lambdarates[lambdacat]
    temp = computepairedll(dataset, params, maskedpaired, lambdarate, lambdacat, siteRate1, siteRate2, siteCat1, siteCat2)
    #temp = fetch(refs[lambdacat])
    if keeplikelihoods
      push!(ret, loglambdaweight+temp)
    end
    pairindex = 1
    for i=1:dataset.numcols
      if paired[i] > i && states[i] == siteCat1 && states[paired[i]] == siteCat2
        if lambdacat == 1
          loglikelihoods[i] = loglambdaweight+temp[pairindex]
        else
          loglikelihoods[i] = CommonUtils.logsumexp(loglikelihoods[i], loglambdaweight+temp[pairindex])
        end
        pairindex += 1
      end
    end
  end
  return loglikelihoods, ret
end

function calculatecachedlikelihoodparallel(dataset::Dataset, params::ModelParameters, musespecificparam::MuseSpecificParameters, paired::Array{Int,1}, siteCat1::Int, siteCat2::Int, pairs::Array{Tuple{Int,Int},1}, cache::FastCache)
    temp = computepairedll(dataset, params, musespecificparam, paired, musespecificparam.lambdacat,  params.siteRates[siteCat1], params.siteRates[siteCat2], siteCat1, siteCat2, pairs, cache)
    return temp, cache
end

function calculatecachedlikelihood(dataset::Dataset, params::ModelParameters, museparams::Array{MuseSpecificParameters,1}, states::Array{Int,1}, paired::Array{Int,1}, cache::Array{FastCache,1}, unpairedcache::Array{FastCache,1}, integratesiterates::Bool, unpairedcolumncache::ColumnCache=ColumnCache(100000), pairedcolumncache::ColumnCache=ColumnCache(1000000))
    if integratesiterates
      unpaired = ones(Float64, dataset.numcols)*-Inf
      pairedloglikelihoods = ones(Float64, dataset.numcols)*-Inf
      cols = Int[]
      for i=1:dataset.numcols
        key = getcolumnkey(dataset.subcolumnrefs, dataset.numcols, Int[i])
        index = getcacheindex(unpairedcolumncache, key)
        if index == 0
          push!(cols, i)
        else
          unpaired[i] = unpairedcolumncache.cache[index]
        end
      end

      pairs = Tuple{Int,Int}[]
      for i=1:length(paired)
        if paired[i] > i
          key = getcolumnkey(dataset.subcolumnrefs, dataset.numcols, Int[i,paired[i]])
          index = getcacheindex(pairedcolumncache, key)
          if index == 0
            push!(pairs, (i,paired[i]))
          else
            pairedloglikelihoods[i] = pairedcolumncache.cache[index]
          end
        end
      end

      unpairedconditionals = zeros(Float64, dataset.numcols, params.siteCats)
      for rateCat=1:params.siteCats
        temp = computeunpairedlikelihoodscat(dataset, params, Int[rateCat for i=1:dataset.numcols], cols, rateCat, unpairedcache[rateCat])
        weight = log(params.siteWeights[rateCat])
        for i=1:dataset.numcols
          unpairedconditionals[i,rateCat] = weight + temp[i]
        end
      end


      ps = Tuple{Int,Int,MuseSpecificParameters}[]
      for siteCat1=1:params.siteCats
        for siteCat2=1:params.siteCats
          for musespecificparam in museparams
            if musespecificparam.lambdarate != 0.0
              push!(ps, (siteCat1, siteCat2, musespecificparam))
            end
          end
        end
      end

      #=
      refs = []
      cacheindex = 1
      cacheindex2 = 1
      active = 0
      fullconditionals = zeros(Float64, dataset.numcols, params.siteCats, params.siteCats, length(museparams))
      for (siteCat1, siteCat2, musespecificparam) in ps

        if active < 4
          fastcache = cache[cacheindex]
          ref = @spawn calculatecachedlikelihoodparallel(dataset, params, musespecificparam, paired, siteCat1, siteCat2, pairs, fastcache)
          cacheindex += 1
          push!(refs, ref)
          active += 1
        end

        for ref in refs
            temp, fastcache = fetch(ref)
            cache[cacheindex2] = fastcache
            weight = log(params.siteWeights[siteCat1]) + log(params.siteWeights[siteCat2]) +  musespecificparam.logprob
            for i=1:length(temp)
              fullconditionals[pairs[i][1], siteCat1, siteCat2, musespecificparam.lambdacat] = weight + temp[i]
            end
            cacheindex2 += 1
        end
        refs = []
        active = 0
      end=#

      cacheindex = 1
      fullconditionals = zeros(Float64, dataset.numcols, params.siteCats, params.siteCats, length(museparams))
      for (siteCat1, siteCat2, musespecificparam) in ps
        weight = log(params.siteWeights[siteCat1]) + log(params.siteWeights[siteCat2]) +  musespecificparam.logprob
        #temp = computepairedll(dataset, params, musespecificparam, paired, musespecificparam.lambdacat,  params.siteRates[siteCat1], params.siteRates[siteCat2], siteCat1, siteCat2, pairs, FastCache(125000,16))
        temp = computepairedll(dataset, params, musespecificparam, paired, musespecificparam.lambdacat,  params.siteRates[siteCat1], params.siteRates[siteCat2], siteCat1, siteCat2, pairs, cache[cacheindex])
        for i=1:length(temp)
          fullconditionals[pairs[i][1], siteCat1, siteCat2, musespecificparam.lambdacat] = weight + temp[i]
        end
        cacheindex += 1
      end



      for i=1:dataset.numcols
          if unpaired[i] == -Inf
            for siteCat=1:params.siteCats
              if siteCat == 1
               unpaired[i] = unpairedconditionals[i, siteCat]
              else
               unpaired[i] = CommonUtils.logsumexp(unpaired[i], unpairedconditionals[i, siteCat])
              end
            end
            key = getcolumnkey(dataset.subcolumnrefs, dataset.numcols, Int[i])
            putvalue(unpairedcolumncache, key, unpaired[i])
          end
      end

      for pair in pairs
        i = pair[1]
        if pairedloglikelihoods[i] == -Inf
          for musespecificparam in museparams
            lambdacat = musespecificparam.lambdacat
            if musespecificparam.lambdarate == 0.0
              pairedloglikelihoods[i] = CommonUtils.logsumexp(pairedloglikelihoods[i], musespecificparam.logprob  + unpaired[i] + unpaired[paired[i]])
            end
          end

          for siteCat1=1:params.siteCats
            for siteCat2=1:params.siteCats
              for musespecificparam in museparams
                lambdacat = musespecificparam.lambdacat
                if musespecificparam.lambdarate == 0.0
                  #ll = musespecificparam.logprob + unpairedconditionals[i, siteCat1] + unpairedconditionals[paired[i],siteCat2]
                  #pairedloglikelihoods[i] = CommonUtils.logsumexp(pairedloglikelihoods[i], ll)
                else
                  pairedloglikelihoods[i] = CommonUtils.logsumexp(pairedloglikelihoods[i], fullconditionals[i, siteCat1, siteCat2, lambdacat])
                end
              end
            end
          end
          key = getcolumnkey(dataset.subcolumnrefs, dataset.numcols, Int[i,paired[i]])
          putvalue(pairedcolumncache, key, pairedloglikelihoods[i])
        end
      end


      unpairedll = 0.0
      pairedll = 0.0
      for i=1:dataset.numcols
        if paired[i] == 0
          unpairedll += unpaired[i]
        elseif paired[i] > i
          pairedll += pairedloglikelihoods[i]
        end
      end

      return unpairedll + pairedll
    else
      unpaired = zeros(Float64, dataset.numcols)
      for rateCat=1:params.siteCats
        temp = computeunpairedlikelihoodscat(dataset, params, states, Int[i for i=1:dataset.numcols], rateCat, unpairedcache[rateCat])
        for i=1:dataset.numcols
          if states[i] == rateCat
            unpaired[i] = temp[i]
          end
        end
      end

      ps = Tuple{Int,Int,MuseSpecificParameters}[]
      for siteCat1=1:params.siteCats
        for siteCat2=1:params.siteCats
          for musespecificparam in museparams
            if musespecificparam.lambdarate != 0.0
              push!(ps, (siteCat1, siteCat2, musespecificparam))
            end
          end
        end
      end

      v = zeros(Float64, dataset.numcols, length(params.lambdarates))
      cacheindex = 1
      for (siteCat1, siteCat2, musespecificparam) in ps
        pairs = Tuple{Int,Int}[]
        for i=1:length(paired)
          if paired[i] > i && states[i] == siteCat1 && states[paired[i]] == siteCat2
            push!(pairs, (i,paired[i]))
          end
        end
        temp = computepairedll(dataset, params, musespecificparam, paired, musespecificparam.lambdacat,  params.siteRates[siteCat1], params.siteRates[siteCat2], siteCat1, siteCat2, pairs, cache[cacheindex])
        for i=1:length(temp)
          v[pairs[i][1],musespecificparam.lambdacat] = musespecificparam.logprob + temp[i]
        end
        cacheindex += 1
      end

      loglikelihoods = zeros(Float64, dataset.numcols)
      for i=1:dataset.numcols
        if paired[i] > i
          for musespecificparam in museparams
            if musespecificparam.lambdarate == 0.0
              loglikelihoods[i] = musespecificparam.logprob  + unpaired[i] + unpaired[paired[i]]
            else
              loglikelihoods[i] = CommonUtils.logsumexp(loglikelihoods[i], v[i,musespecificparam.lambdacat])
            end
          end
        end
      end

      for i=1:dataset.numcols
        if paired[i] != 0
          unpaired[i] = 0.0
        end
      end

      total = sum(unpaired) + sum(loglikelihoods)
      return total
    end
end

function computeunpairedlikelihoods(dataset::Dataset, params::ModelParameters, states::Array{Int,1})
  unpaired = SharedArray{Float64}(dataset.numcols)
  #@sync @parallel for rateCat=1:params.siteCats
  for rateCat=1:params.siteCats
    temp = computeunpairedlikelihoodscat(dataset, params, states, Int[i for i=1:dataset.numcols], rateCat)
    #weight = log(params.siteWeights[rateCat])
    for i=1:dataset.numcols
      if states[i] == rateCat
        unpaired[i] = temp[i]
      end
    end
  end
  unpaired2 = zeros(Float64, dataset.numcols)
  for i=1:dataset.numcols
    unpaired2[i] = unpaired[i]
  end
  return unpaired2
end

function computeunpairedlikelihoods(dataset::Dataset, params::ModelParameters)
  unpairedconditionals = SharedArray{Float64}(dataset.numcols, params.siteCats)
  #@sync @parallel for rateCat=1:params.siteCats
  for rateCat=1:params.siteCats
    temp = computeunpairedlikelihoodscat(dataset, params, Int[rateCat for i=1:dataset.numcols], Int[i for i=1:dataset.numcols], rateCat)
    weight = log(params.siteWeights[rateCat])
    for i=1:dataset.numcols
      unpairedconditionals[i,rateCat] = weight + temp[i]
    end
  end

  unpaired = zeros(Float64, dataset.numcols)
  for siteCat=1:params.siteCats
    for i=1:dataset.numcols
      if siteCat == 1
        unpaired[i] = unpairedconditionals[i, siteCat]
      else
        unpaired[i] = CommonUtils.logsumexp(unpaired[i], unpairedconditionals[i, siteCat])
      end
    end
  end
  return unpaired
end

function computeunpairedlikelihoods(dataset::Dataset, params::ModelParameters, museparams::Array{MuseSpecificParameters,1}, states::Array{Int,1}, paired::Array{Int,1}, B::Float64=1.0, conditionalRates::Bool=true, rng::AbstractRNG=MersenneTwister(291019011901), likelihoodsOnly::Bool=false)
  unpairedconditionals = SharedArray{Float64}(dataset.numcols, params.siteCats)
  #@sync @parallel for rateCat=1:params.siteCats
  for rateCat=1:params.siteCats
    temp = computeunpairedlikelihoodscat(dataset, params, Int[rateCat for i=1:dataset.numcols], Int[i for i=1:dataset.numcols], rateCat)
    weight = log(params.siteWeights[rateCat])
    for i=1:dataset.numcols
      unpairedconditionals[i,rateCat] = weight + temp[i]
    end
  end

  unpairedcolumncache = ColumnCache(10000)
  unpaired = zeros(Float64, dataset.numcols)
  for i=1:dataset.numcols
    if paired[i] == 0
      key = getcolumnkey(dataset.subcolumnrefs, dataset.numcols, Int[i])
      index = getcacheindex(unpairedcolumncache, key)
      if index == 0
        for siteCat=1:params.siteCats
          if siteCat == 1
            unpaired[i] = unpairedconditionals[i, siteCat]
          else
            unpaired[i] = CommonUtils.logsumexp(unpaired[i], unpairedconditionals[i, siteCat])
          end
        end
        putvalue(unpairedcolumncache, key, unpaired[i])
      else
        unpaired[i] = unpairedcolumncache.cache[index]
      end
    end
  end

  return unpaired
end

function computeuncachedlikelihood(dataset::Dataset, params::ModelParameters, museparams::Array{MuseSpecificParameters,1}, states::Array{Int,1}, paired::Array{Int,1}, B::Float64=1.0, conditionalRates::Bool=true, rng::AbstractRNG=MersenneTwister(291019011901), likelihoodsOnly::Bool=false)
  if conditionalRates
    unpaired = SharedArray{Float64}(dataset.numcols)
    #@sync @parallel for rateCat=1:params.siteCats
    for rateCat=1:params.siteCats
      temp = computeunpairedlikelihoodscat(dataset, params, states, Int[i for i=1:dataset.numcols], rateCat)
      #weight = log(params.siteWeights[rateCat])
      for i=1:dataset.numcols
        if states[i] == rateCat
          unpaired[i] = temp[i]
        end
      end
    end

    ps = Tuple{Int,Int,MuseSpecificParameters}[]
    for siteCat1=1:params.siteCats
      for siteCat2=1:params.siteCats
        for musespecificparam in museparams
          if musespecificparam.lambdarate != 0.0
            push!(ps, (siteCat1, siteCat2, musespecificparam))
          end
        end
      end
    end
    v = SharedArray{Float64}(dataset.numcols, length(museparams))
    #@sync @parallel for (siteCat1, siteCat2, musespecificparam) in ps
    for (siteCat1, siteCat2, musespecificparam) in ps
      pairs = Tuple{Int,Int}[]
      for i=1:length(paired)
        if paired[i] > i && states[i] == siteCat1 && states[paired[i]] == siteCat2
          push!(pairs, (i,paired[i]))
        end
      end
      loglambdaweight = musespecificparam.logprob
      temp = computepairedll(dataset, params, musespecificparam, paired, musespecificparam.lambdacat,  params.siteRates[siteCat1], params.siteRates[siteCat2], siteCat1, siteCat2, pairs,FastCache(125000,16))
      for i=1:length(temp)
        v[pairs[i][1],musespecificparam.lambdacat] = loglambdaweight + temp[i]
      end
    end

    loglikelihoods = zeros(Float64, dataset.numcols)
    for i=1:dataset.numcols
      if paired[i] > i
        for musespecificparam in museparams
          if musespecificparam.lambdarate == 0.0
            loglikelihoods[i] = musespecificparam.logprob + unpaired[i] + unpaired[paired[i]]
          else
            loglikelihoods[i] = CommonUtils.logsumexp(loglikelihoods[i], v[i,musespecificparam.lambdacat])
          end
        end
      end
    end

    for i=1:dataset.numcols
      if paired[i] != 0
        unpaired[i] = 0.0
      end
    end

    return unpaired,loglikelihoods
  else
    unpairedconditionals = SharedArray{Float64}(dataset.numcols, params.siteCats)
    #@sync @parallel for rateCat=1:params.siteCats
    for rateCat=1:params.siteCats
      temp = computeunpairedlikelihoodscat(dataset, params, Int[rateCat for i=1:dataset.numcols], Int[i for i=1:dataset.numcols], rateCat)
      weight = log(params.siteWeights[rateCat])
      for i=1:dataset.numcols
        unpairedconditionals[i,rateCat] = weight + temp[i]
      end
    end


    ps = Tuple{Int,Int,MuseSpecificParameters}[]
    for siteCat1=1:params.siteCats
      for siteCat2=1:params.siteCats
        for musespecificparam in museparams
          if musespecificparam.lambdarate != 0.0
            push!(ps, (siteCat1, siteCat2, musespecificparam))
          end
        end
      end
    end

    fullconditionals = SharedArray{Float64}(dataset.numcols, params.siteCats, params.siteCats, length(museparams))
    #@sync @parallel for (siteCat1, siteCat2, musespecificparam) in ps
    for (siteCat1, siteCat2, musespecificparam) in ps
      pairs = Tuple{Int,Int}[]
      for i=1:length(paired)
        if paired[i] > i
          push!(pairs, (i,paired[i]))
        end
      end
      weight = log(params.siteWeights[siteCat1]) + log(params.siteWeights[siteCat2]) +  musespecificparam.logprob
      temp = computepairedll(dataset, params, musespecificparam, paired, musespecificparam.lambdacat,  params.siteRates[siteCat1], params.siteRates[siteCat2], siteCat1, siteCat2, pairs)
      for i=1:length(temp)
        fullconditionals[pairs[i][1], siteCat1, siteCat2, musespecificparam.lambdacat] = weight + temp[i]
      end
    end

    unpairedcolumncache = ColumnCache(100000)
    pairedcolumncache = ColumnCache(100000)


    pairconditionals = ones(Float64, dataset.numcols, params.siteCats*params.siteCats)*-Inf
    museconditionals = ones(Float64, dataset.numcols, length(museparams))*-Inf
    pairedloglikelihoods = ones(Float64, dataset.numcols)*-Inf
    for i=1:dataset.numcols
      if paired[i] > i
        key = getcolumnkey(dataset.subcolumnrefs, dataset.numcols, Int[i,paired[i]])
        index = getcacheindex(pairedcolumncache, key)
        if index == 0 ||  !likelihoodsOnly
          for siteCat1=1:params.siteCats
            for siteCat2=1:params.siteCats
              rateindex = (siteCat1-1)*params.siteCats + siteCat2
              for musespecificparam in museparams
                lambdacat = musespecificparam.lambdacat
                if musespecificparam.lambdarate == 0.0
                  ll = musespecificparam.logprob + unpairedconditionals[i, siteCat1] + unpairedconditionals[paired[i],siteCat2]
                  pairedloglikelihoods[i] = CommonUtils.logsumexp(pairedloglikelihoods[i], ll)
                  if !likelihoodsOnly
                    pairconditionals[i, rateindex] = CommonUtils.logsumexp(pairconditionals[i, rateindex], ll)
                    museconditionals[i,lambdacat] = CommonUtils.logsumexp(museconditionals[i,lambdacat], ll)
                  end
                else
                  pairedloglikelihoods[i] = CommonUtils.logsumexp(pairedloglikelihoods[i], fullconditionals[i, siteCat1, siteCat2, lambdacat])
                  if !likelihoodsOnly
                    pairconditionals[i, rateindex] = CommonUtils.logsumexp(pairconditionals[i, rateindex], fullconditionals[i, siteCat1, siteCat2, lambdacat])
                    museconditionals[i,lambdacat] = CommonUtils.logsumexp(museconditionals[i,lambdacat], fullconditionals[i, siteCat1, siteCat2, lambdacat])
                  end
                end
              end
            end
          end
          putvalue(pairedcolumncache, key, pairedloglikelihoods[i])
        else
          pairedloglikelihoods[i] = pairedcolumncache.cache[index]
        end
      end
    end

    #index = getcacheindex(unpairedcolumncache, i)
    #=
   unpaired = zeros(Float64, dataset.numcols)
   for siteCat=1:params.siteCats
     for i=1:dataset.numcols
       if paired[i] == 0
         if siteCat == 1
           unpaired[i] = unpairedconditionals[i, siteCat]
         else
           unpaired[i] = CommonUtils.logsumexp(unpaired[i], unpairedconditionals[i, siteCat])
         end
       end
     end
    end

    if !likelihoodsOnly
      for i=1:dataset.numcols
        if paired[i] == 0
          unpairedconditionals[i,:] = exp(unpairedconditionals[i,:] - unpaired[i])
        end
      end

      for i=1:dataset.numcols
        if paired[i] > i
          for siteCat1=1:params.siteCats
            for siteCat2=1:params.siteCats
              pairconditionals[i, (siteCat1-1)*params.siteCats + siteCat2] = exp(pairconditionals[i, (siteCat1-1)*params.siteCats + siteCat2] - pairedloglikelihoods[i])
            end
          end
        end
      end
    end=#

    unpaired = zeros(Float64, dataset.numcols)
    for i=1:dataset.numcols
      if paired[i] == 0
        key = getcolumnkey(dataset.subcolumnrefs, dataset.numcols, Int[i])
        index = getcacheindex(unpairedcolumncache, key)
        if index == 0
          for siteCat=1:params.siteCats
            if siteCat == 1
              unpaired[i] = unpairedconditionals[i, siteCat]
            else
              unpaired[i] = CommonUtils.logsumexp(unpaired[i], unpairedconditionals[i, siteCat])
            end
          end
          putvalue(unpairedcolumncache, key, unpaired[i])
        else
          unpaired[i] = unpairedcolumncache.cache[index]
        end
      end
    end

   if !likelihoodsOnly
     for i=1:dataset.numcols
       if paired[i] == 0
         unpairedconditionals[i,:] = exp.(unpairedconditionals[i,:] .- unpaired[i])
       end
     end

     for i=1:dataset.numcols
       if paired[i] > i
         for siteCat1=1:params.siteCats
           for siteCat2=1:params.siteCats
             pairconditionals[i, (siteCat1-1)*params.siteCats + siteCat2] = exp.(pairconditionals[i, (siteCat1-1)*params.siteCats + siteCat2] .- pairedloglikelihoods[i])
           end
         end
       end
     end
   end


    sampledstates = zeros(Int, dataset.numcols)
    if !likelihoodsOnly
      for i=1:dataset.numcols
        if paired[i] > i
          r = CommonUtils.sample(rng, pairconditionals[i, :].^B)
          r1 = div(r-1, params.siteCats)
          r2 = (r-1) % params.siteCats
          sampledstates[i] = r1 + 1
          sampledstates[paired[i]] = r2 + 1
        elseif paired[i] == 0
          sampledstates[i] = CommonUtils.sample(rng, unpairedconditionals[i, :].^B)
        end
      end
    end

    return unpaired,pairedloglikelihoods,sampledstates,museconditionals
  end
end

function computedpairedlikelihoodscats(dataset::Dataset, params::ModelParameters, states::Array{Int,1}, paired::Array{Int,1}, parallel::Bool=true, keeplikelihoods::Bool=false)
  indices = Tuple{Int,Int}[]
  for siteCat1=1:params.siteCats
    for siteCat2=1:params.siteCats
      push!(indices, (siteCat1, siteCat2))
    end
  end
  temps = pmap(index -> computedpairedlikelihoodscatshelper(dataset, params, states, paired, index[1], index[2], keeplikelihoods), indices)

  index = 1
  loglikelihoods = zeros(Float64,dataset.numcols)
  ret = nothing
  if keeplikelihoods
    ret = zeros(Float64,dataset.numcols,dataset.numcols, params.lambdaCats)
  end
  for siteCat1=1:params.siteCats
    for siteCat2=1:params.siteCats

      if keeplikelihoods
        temp2 = temps[index][2]
        pairindex = 1
        for i=1:dataset.numcols
          if paired[i] > i && states[i] == siteCat1 && states[paired[i]] == siteCat2
            for j=1:params.lambdaCats
              ret[i,j] = temp2[j][pairindex]
            end
            pairindex += 1
          end
        end
      end

      temp = temps[index][1]
      for i=1:dataset.numcols
        if paired[i] > i && states[i] == siteCat1 && states[paired[i]] == siteCat2
          loglikelihoods[i] = temp[i]
        end
      end
      index += 1
    end
  end
  return loglikelihoods, ret
end

function discretizegamma(shape::Float64, scale::Float64, numcategories::Int)
  if numcategories == 1
    return Float64[1.0]
  else
    catwidth = 1.0 / numcategories
    vs = Float64[catwidth/2.0 + (i-1.0)*catwidth for i=1:numcategories]
    gammadist = Gamma(shape, scale)
    return Float64[quantile(gammadist, v) for v in vs]
  end
end

function discretizegamma2(lambdazeroweight::Float64, p1::Float64, p2::Float64, numcategories::Int, parameterisation::Int=1)
  shape = p1
  scale = p2
  if parameterisation == 1
    gammamean = p1
    gammavar = p2
    shape = gammamean/gammavar
    scale = gammavar
  elseif parameterisation == 2
    gammamean = p1
    gammaprecision = 1.0 / p2
    shape = gammamean/gammaprecision
    scale = gammaprecision
  end

  lambdarates = Float64[0.0]
  append!(lambdarates, discretizegamma(shape, scale, numcategories))
  weights = ones(Float64, numcategories+1)
  weights /= numcategories
  weights *= (1.0-lambdazeroweight)
  weights[1] = lambdazeroweight
  return weights, lambdarates
end

function discretizegamma3(lambdazeroweight::Float64, meanGC::Float64, meanAT::Float64, meanGT::Float64, gammavar::Float64, numcategories::Int, parameterisation::Int=0)

  rates =  Float64[0.0]
  append!(rates, discretizegamma(gammavar, 1.0/gammavar, numcategories))

  lambdaratesGC = rates*(meanGC-1.0)
  lambdaratesAT = rates*(meanAT-1.0)
  lambdaratesGT = rates*(meanGT-1.0)
  #println("lambdaratesGC", lambdaratesGC)
  #println("lambdaratesAT", lambdaratesAT)
  #println("lambdaratesGT", lambdaratesGT)

  weights = ones(Float64, numcategories+1)
  weights /= numcategories
  weights *= (1.0-lambdazeroweight)
  weights[1] = lambdazeroweight
  return weights, lambdaratesGC, weights, lambdaratesAT, weights, lambdaratesGT
end

function samplelambdaweights(rng::AbstractRNG, dataset::Dataset, params::ModelParameters, paired::Array{Int,1}, B::Float64=1.0, c::Float64=0.5)
  refs = []
  ll = -Inf
  for lambdacat=1:length(params.lambdarates)
    lambdarate = params.lambdarates[lambdacat]
    #ref = @spawn computepairedll(dataset, params, paired, lambdarate, lambdacat)
    ref = computepairedll(dataset, params, paired, lambdarate, lambdacat)
    push!(refs, ref)
  end

  likelihoods = Array{Float64,1}[]
  for ref in refs
    #push!(likelihoods, fetch(ref))
    push!(likelihoods, ref)
  end

  loglambdaweights = log(params.lambdaweightsGC)
  numratecats = length(params.lambdarates)
  numpairs = length(likelihoods[1])
  v = zeros(Float64, numratecats)
  alloc = ones(Float64, numratecats)*c
  for i=1:numpairs
    total = -Inf
    for r=1:numratecats
      v[r] = (loglambdaweights[r] + likelihoods[r][i])*B
      total = CommonUtils.logsumexp(total, v[r])
    end
    v = exp.(v.-total)
    alloc[CommonUtils.sample(rng, v)] += 1.0
  end
  return rand(Dirichlet(B*alloc-B+1.0))
end

include("FelsensteinCUDA.jl")


function computepairedll(dataset::Dataset, params::ModelParameters, musespecificparams::MuseSpecificParameters, paired::Array{Int,1}, lambdacat::Int, siteRate1::Float64, siteRate2::Float64, siteCat1::Int, siteCat2::Int, pairs::Array{Tuple{Int,Int}}=Tuple{Int,Int}[], fastcache::FastCache=FastCache(125000,16))
    musemodel = MuseModel(params.freqs, getGC(params,musespecificparams), getAT(params,musespecificparams), getGT(params,musespecificparams), params.q1, params.q2, params.q3, params.q4, params.q5, params.q6, siteRate1, siteRate2)
    transprobs = gettransitionmatriceslist(params.branchlengths, musemodel.Q)
    #pairs = Tuple{Int,Int}[]
    if length(pairs) == 0
      for i=1:length(paired)
        if paired[i] > i
          push!(pairs,(i,paired[i]))
        end
      end
    end
    return coevolution(dataset, musemodel.freqs, transprobs, pairs, fastcache)
end

function coevolution(dataset::Dataset, freqs::Array{Float64,1}, transprobs::Array{Array{Float64,2},1}, pairs::Array{Tuple{Int,Int},1}, fastcache::FastCache)
  root = dataset.root
  alphabet = 16
  cache = Array{Float64,1}[zeros(Float64,alphabet) for i=1:dataset.numnodes]
  loglikcache = Dict{Tuple,Float64}()
  nodelist = getnodelist(root)

  likelihoods = Float64[]
  for (i,j) in pairs
    columns = Int[i,j]
    key = (0, getkey(dataset.numnodes, root.data.nodeindex, dataset.subcolumnrefs, dataset.numcols, columns))
    total = get(loglikcache,key,0.0)
    if total == 0.0
      for z=1:dataset.numnodes
        cache[z][1] = -Inf
      end
      logm = zeros(Float64, dataset.numnodes)
      v = felsensteinstack(dataset, nodelist, cache, logm, fastcache, dataset.subcolumnrefs, transprobs, root, dataset.data, columns, alphabet)
      for a=1:alphabet
        total += freqs[a] * v[a]
      end
      total = log(total) + logm[1]
      loglikcache[key] = total
    end
    push!(likelihoods,total)
  end

  return likelihoods
end

function getposteriorsiterates(singleprobs::Array{Float64,1}, pairprobs::Array{Float64,2}, dataset::Dataset, params::ModelParameters, unpaired::Bool, parallel::Bool=true, keepmatrices=false, integratesiterates::Bool=false, maxbasepairdistance::Int=1000000)
  siterateconditionals = ones(Float64,params.siteCats,dataset.numcols)*-Inf
  refs = []
  for rateCat=1:params.siteCats
    #ref = @spawn computeunpairedlikelihoodscat(dataset, params, Int[rateCat for i=1:dataset.numcols], Int[i for i=1:dataset.numcols], rateCat)
    ref = computeunpairedlikelihoodscat(dataset, params, Int[rateCat for i=1:dataset.numcols], Int[i for i=1:dataset.numcols], rateCat)
    push!(refs,ref)
  end

  for rateCat=1:params.siteCats
    #temp = @fetch(refs[rateCat])
    temp = refs[rateCat]
    for col=1:dataset.numcols
      if temp[col] != -Inf
        singleprob = 0.0
        if !unpaired
          singleprob = log(singleprobs[col])
        end
        siterateconditionals[rateCat,col] =  CommonUtils.logsumexp(siterateconditionals[rateCat, col], log(params.siteWeights[rateCat]) + singleprob + temp[col])
      end
    end
  end

  if !unpaired
    museparams = getmusespecificparamsarray(params)
    refs = []
    for j=1:params.siteCats
      for k=1:params.siteCats
        for musespecificparam in museparams
          #ref = @spawn coevolutionall(dataset, params, musespecificparam, params.states, true, params.siteRates[j], params.siteRates[k], j, k, 1, min(dataset.numcols, maxbasepairdistance), parallel)
          ref = coevolutionall(dataset, params, musespecificparam, params.states, true, params.siteRates[j], params.siteRates[k], j, k, 1, min(dataset.numcols, maxbasepairdistance), parallel)
          push!(refs,ref)
        end
      end
    end
    ret = PartialMatrix[]
    pairedcolumncache = ColumnCache(1000000)
    for ref in refs
        push!(ret, ref)
      #push!(ret, fetch(ref))
    end

    columns = Int[0,0]
    for x=1:dataset.numcols
      ystart = max(1, x-maxbasepairdistance)
      yend = min(x+maxbasepairdistance,dataset.numcols)
      for y=x+1:yend
        if x != y
          columns[1] = x
          columns[2] = y
          logpairprobsxy =  log(pairprobs[x,y])
          i = 1
          for j=1:params.siteCats
            logj = log(params.siteWeights[j])
            for k=1:params.siteCats
              logk = log(params.siteWeights[k])
              for musespecificparam in museparams
                loglambdaweight = musespecificparam.logprob + logj + logk + logpairprobsxy
                siterateconditionals[j,x] = CommonUtils.logsumexp(siterateconditionals[j,x], loglambdaweight+ret[i][x,y])
                siterateconditionals[k,y] = CommonUtils.logsumexp(siterateconditionals[k,y], loglambdaweight+ret[i][x,y])
                i += 1
              end
            end
          end
        end
      end
    end
  end
  return siterateconditionals
end

function coevolutionall(dataset::Dataset, params::ModelParameters, parallel::Bool=true, keepmatrices=false, integratesiterates::Bool=false, maxbasepairdistance::Int=1000000, usecuda::Bool=true)
  if integratesiterates && usecuda
    museparams = getmusespecificparamsarray(params)
    mat, ret = felsensteincuda(dataset, params, museparams,maxbasepairdistance,keepmatrices)
    GC.gc()
    return mat,ret
  else
    #tic()
    museparams = getmusespecificparamsarray(params)
    refs = []
    for j=1:params.siteCats
      for k=1:params.siteCats
        for musespecificparam in museparams
          #ref = @spawn coevolutionall(dataset, params, musespecificparam, params.states, integratesiterates, params.siteRates[j], params.siteRates[k], j, k, 1, min(dataset.numcols, maxbasepairdistance), parallel)
          ref = coevolutionall(dataset, params, musespecificparam, params.states, integratesiterates, params.siteRates[j], params.siteRates[k], j, k, 1, min(dataset.numcols, maxbasepairdistance), parallel)
          push!(refs,ref)
        end
      end
    end


    #key = getcolumnkey(dataset.subcolumnrefs, dataset.numcols, Int[x,y])*params.siteCats*params.siteCats + j*params.siteCats + k
    #index = getcacheindex(pairedcolumncache, key)
    coevolutionll = ones(Float64, dataset.numcols, dataset.numcols)*-Inf
    ret = PartialMatrix[]
    if integratesiterates
      pairedcolumncache = ColumnCache(1000000)
      for ref in refs
        #push!(ret, fetch(ref))
        push!(ret, ref)
      end
      #elapsed = #toc()
      #println("Computing matrices", elapsed)
      #tic()
      columns = Int[0,0]
      for x=1:dataset.numcols
        yend = min(x+maxbasepairdistance,dataset.numcols)
        for y=x+1:yend
          columns[1] = x
          columns[2] = y
          key = getcolumnkey(dataset.subcolumnrefs, dataset.numcols, columns)
          index = getcacheindex(pairedcolumncache, key)
          if index == 0
            i = 1
            for j=1:params.siteCats
              logj = log(params.siteWeights[j])
              for k=1:params.siteCats
                logk = log(params.siteWeights[k])
                for musespecificparam in museparams
                  loglambdaweight = musespecificparam.logprob + logj + logk
                  coevolutionll[x,y] = CommonUtils.logsumexp(coevolutionll[x,y], loglambdaweight+ret[i][x,y])
                  i += 1
                end
              end
            end
            coevolutionll[y,x] = coevolutionll[x,y]
            putvalue(pairedcolumncache, key, coevolutionll[x,y])
          else
            coevolutionll[x,y] = pairedcolumncache.cache[index]
            coevolutionll[y,x] = coevolutionll[x,y]
          end
        end
      end
      #elapsed = #toc()
      #println("Computing final", elapsed)
    else
      pairedcolumncache = ColumnCache(1000000)
      for ref in refs
        push!(ret, fetch(ref))
      end
      columns = Int[0,0]
      for x=1:dataset.numcols
        yend = min(x+maxbasepairdistance,dataset.numcols)
        for y=x+1:yend
          columns[1] = x
          columns[2] = y
          j = params.states[x]
          k = params.states[y]
          key = (getcolumnkey(dataset.subcolumnrefs, dataset.numcols, columns)-1)*params.siteCats*params.siteCats + (j-1)*params.siteCats + k
          index = getcacheindex(pairedcolumncache, key)
          if index == 0
            i = (j-1)*params.siteCats*length(museparams) + (k-1)*length(museparams)
            for musespecificparam in museparams
              coevolutionll[x,y] =  CommonUtils.logsumexp(coevolutionll[x,y], musespecificparam.logprob+ret[i+musespecificparam.lambdacat][x,y])
            end
            coevolutionll[y,x] = coevolutionll[x,y]
            putvalue(pairedcolumncache, key, coevolutionll[x,y])
          else
            coevolutionll[x,y] = pairedcolumncache.cache[index]
            coevolutionll[y,x] = coevolutionll[x,y]
          end
        end
      end
    end

    #=
    ret2 = Array{Float64,2}[]
    if keepmatrices
      ret2 = Array{Float64,2}[ones(Float64, dataset.numcols, dataset.numcols)*-Inf for dummy in museparams]
      index = 1
      for j=1:params.siteCats
        for k=1:params.siteCats
          lambdaindex = 1
          for musespecificparam in museparams
            for x=1:dataset.numcols
              for y=1:dataset.numcols
                ret2[lambdaindex][x,y] = max(ret2[lambdaindex][x,y], ret[index][x,y])
              end
            end
            index += 1
            lambdaindex += 1
          end
        end
      end
    end=#
    ret2 = PartialMatrix[]
    if keepmatrices
      ret2 = ret
    end

    return coevolutionll,ret2
  end
end

function coevolutionall(dataset::Dataset, params::ModelParameters, musespecificparams::MuseSpecificParameters, states::Array{Int,1}, integratesiterates::Bool, siteRate1::Float64, siteRate2::Float64, siteCat1::Int, siteCat2::Int, numparts::Int,maxbasepairdistance::Int, parallel::Bool=true, evolutionmode::Int=0)

  if evolutionmode == 1
    musemodel = MuseModel(params.freqs, getGC(params,musespecificparams), getAT(params,musespecificparams), getGT(params,musespecificparams), params.q1, params.q2, params.q3, params.q4, params.q5, params.q6, siteRate1, siteRate2)
    freqs = musemodel.freqs
    coevolutionll = PartialMatrix(dataset.numcols,maxbasepairdistance)
    for x=1:dataset.numcols
      yend = min(x+maxbasepairdistance,dataset.numcols)
      for y=x+1:yend
        for seq=1:dataset.numseqs
          coevolutionll[x,y] += log(state(dataset.data, seq, Int[x,y]).*freqs)
        end
      end
    end
    return coevolutionll
  elseif evolutionmode == 2
    musemodel = MuseModel(params.freqs, getGC(params,musespecificparams), getAT(params,musespecificparams), getGT(params,musespecificparams), params.q1, params.q2, params.q3, params.q4, params.q5, params.q6, siteRate1, siteRate2)
    freqs = musemodel.freqs
    coevolutionll = PartialMatrix(dataset.numcols,maxbasepairdistance)
    temp1 = computeunpairedlikelihoodscat(dataset, params, Int[siteCat1 for i=1:dataset.numcols], Int[i for i=1:dataset.numcols], siteCat1)
    temp2 = computeunpairedlikelihoodscat(dataset, params, Int[siteCat2 for i=1:dataset.numcols], Int[i for i=1:dataset.numcols], siteCat2)
    for x=1:dataset.numcols
      yend = min(x+maxbasepairdistance,dataset.numcols)
      for y=x+1:yend
        coevolutionll[x,y] = temp1[x] + temp2[y]
      end
    end
    return coevolutionll
  end


  if musespecificparams.lambdarate == 0.0
    coevolutionll = PartialMatrix(dataset.numcols,maxbasepairdistance)
    temp1 = computeunpairedlikelihoodscat(dataset, params, Int[siteCat1 for i=1:dataset.numcols], Int[i for i=1:dataset.numcols], siteCat1)
    temp2 = computeunpairedlikelihoodscat(dataset, params, Int[siteCat2 for i=1:dataset.numcols], Int[i for i=1:dataset.numcols], siteCat2)
    for x=1:dataset.numcols
      yend = min(x+maxbasepairdistance,dataset.numcols)
      for y=x+1:yend
        coevolutionll[x,y] = temp1[x] + temp2[y]
      end
    end
    return coevolutionll
  else
    musemodel = MuseModel(params.freqs, getGC(params,musespecificparams), getAT(params,musespecificparams), getGT(params,musespecificparams), params.q1, params.q2, params.q3, params.q4, params.q5, params.q6, siteRate1, siteRate2)


    #musemodel = MuseModel(params.freqs, exp(log(params.lambdaGC)*lambdarate), exp(log(params.lambdaAT)*lambdarate), exp(log(params.lambdaGT)*lambdarate), params.q1, params.q2, params.q3, params.q4, params.q5, params.q6, siteRate1, siteRate2)
    freqs = musemodel.freqs
    Q = musemodel.Q
    alphabet = 16

    root = dataset.root
    transprobs = gettransitionmatriceslist(params.branchlengths, Q)
    partitions = getpartitions(dataset.numcols,numparts)

    return felsensteinall(dataset, params, freqs, transprobs, alphabet, states, integratesiterates, siteCat1, siteCat2, 1, dataset.numcols, maxbasepairdistance)
  end

  #=
  musemodel = MuseModel(params.freqs, getGC(params,musespecificparams), getAT(params,musespecificparams), getGT(params,musespecificparams), params.q1, params.q2, params.q3, params.q4, params.q5, params.q6, siteRate1, siteRate2)
  #musemodel = MuseModel(params.freqs, exp(log(params.lambdaGC)*lambdarate), exp(log(params.lambdaAT)*lambdarate), exp(log(params.lambdaGT)*lambdarate), params.q1, params.q2, params.q3, params.q4, params.q5, params.q6, siteRate1, siteRate2)
  freqs = musemodel.freqs
  Q = musemodel.Q
  alphabet = 16

  root = dataset.root
  transprobs = gettransitionmatriceslist(params.branchlengths, Q)
  partitions = getpartitions(dataset.numcols,numparts)

  refs = Any[]
  for part in partitions
    ref = @spawn felsensteinall(dataset, params, freqs, transprobs, alphabet, states, integratesiterates, siteCat1, siteCat2, part[1], part[2],maxbasepairdistance)
    push!(refs,ref)
  end
  coevolutionll = ones(Float64, dataset.numcols, dataset.numcols)*-Inf
  for p=1:length(partitions)
    coevolutionpart = fetch(refs[p])
    for i=partitions[p][1]:partitions[p][2]
      for j=i+1:dataset.numcols
        v = coevolutionpart[i-partitions[p][1]+1,j]
        if v != -Inf
          coevolutionll[i,j] = v
          coevolutionll[j,i] = v
        end
      end
    end
  end

  return coevolutionll=#
  #=
  if parallel
    refs = Any[]
    for part in partitions
      #println("A")
      ref = @spawn felsensteinall(dataset, params, freqs, transprobs, alphabet, siteCat1, siteCat2, part[1], part[2])
      push!(refs,ref)
    end
    coevolutionll = ones(Float64, dataset.numcols, dataset.numcols)*-Inf
    for p=1:length(partitions)
      #println("B",p)
      coevolutionpart = fetch(refs[p])
      #println("W",coevolutionpart)
      for i=partitions[p][1]:partitions[p][2]
        for j=i+1:dataset.numcols
          v = coevolutionpart[i-partitions[p][1]+1,j]
          if v != -Inf
            coevolutionll[i,j] = v
            coevolutionll[j,i] = v
          end
        end
      end
    end
    #println("D")
    return coevolutionll
  else
    coevolutionll = ones(Float64, dataset.numcols, dataset.numcols)*-Inf
    for p=1:length(partitions)
      part = partitions[p]
      coevolutionpart = felsensteinall(dataset, params, freqs, transprobs, alphabet, siteCat1, siteCat2, part[1], part[2])
      for i=partitions[p][1]:partitions[p][2]
        for j=i+1:dataset.numcols
          coevolutionll[i,j] = coevolutionpart[i-partitions[p][1]+1,j]
          coevolutionll[j,i] = coevolutionll[i,j]
        end
      end
    end
    return coevolutionll
  end
  =#

end

function felsensteinall(dataset::Dataset, params::ModelParameters, freqs::Array{Float64,1}, transprobs::Array{Array{Float64,2},1}, alphabet::Int, states::Array{Int,1}, integratesiterates::Bool=false, siteCat1::Int=1, siteCat2::Int=1, startindex::Int=1, endindex::Int=0, maxbasepairdistance::Int=1000000)
  if endindex == 0
    endindex = dataset.numcols
  end
  root = dataset.root
  nodelist = getnodelist(root)
  revnodelist = reverse(nodelist)
  subcolumnrefs = dataset.subcolumnrefs
  pairedcolumncache = ColumnCache(1000000)
  fastcache = FastCache(1000000,alphabet)



  cache = Array{Float64,1}[zeros(Float64,alphabet) for i=1:dataset.numnodes]
  count = 0
  data = dataset.data
  coevolutionll = PartialMatrix(dataset.numcols,maxbasepairdistance)
  #coevolutionll = ones(Float64, endindex-startindex+1, dataset.numcols)*-Inf
  columns = Int[0,0]
  logm = zeros(Float64, dataset.numnodes)
  for i=startindex:endindex
    if integratesiterates || states[i] == siteCat1
      jend = min(i+maxbasepairdistance, dataset.numcols)
      for j=i+1:jend
        if  integratesiterates || states[j] == siteCat2
          if abs(i-j) <= maxbasepairdistance
            columns[1] = i
            columns[2] = j
            key = getcolumnkey(dataset.subcolumnrefs, dataset.numcols, columns)
            index = getcacheindex(pairedcolumncache, key)
            total = 0.0
            if index == 0
              for z=1:dataset.numnodes
                cache[z][1] = -Inf
              end
              fill!(logm,0.0)
              v = felsensteinstack(dataset, nodelist, cache, logm, fastcache, subcolumnrefs, transprobs, root, data, columns, alphabet)
              for a=1:alphabet
                total += freqs[a] * v[a]
              end
              total = log(total) + logm[1]
              putvalue(pairedcolumncache, key, total)
            else
              total = pairedcolumncache.cache[index]
            end
            coevolutionll[i-startindex+1,j] = total
          end
        end
        count += 1
      end
    end
  end

  return coevolutionll
end

function felsensteinterative(dataset::Dataset, nodelist::Array{TreeNode,1}, cache::Array{Float64,2}, subcolumnloglikcache::Array{Dict{Int,Array{Float64,1}},1}, subcolumnrefs::Array{Int,2}, logtransprobs::Array{Float64,3}, logdata::Array{Float64,3}, columns::Array{Int,1})
  for node in nodelist
    #println(node.name, "\t", node.data.nodeindex,"\t",node.branchlength)
    if isleafnode(node)
      #v = logstate(logdata, node.data.seqindex,columns)
      #subcolumnloglikcache[node.data.nodeindex][key] = v
      cache[node.data.nodeindex,:] = logstate(logdata, node.data.seqindex, columns)
      #return v
    else
      key = (subcolumnrefs[node.data.nodeindex,columns[1]]-1)*dataset.numcols + subcolumnrefs[node.data.nodeindex,columns[2]]
      if haskey(subcolumnloglikcache[node.data.nodeindex],key)
        cache[node.data.nodeindex,:] = subcolumnloglikcache[node.data.nodeindex][key]
        #return subcolumnloglikcache[node.data.nodeindex][key]
      else
        v = zeros(Float64,16)
        for a=1:16
          bsum = -Inf
          csum = -Inf
          for d=1:16
            #bstore[d] = logtransprobs[node.children[1].data.nodeindex,a,d] + felsenstein(cache, logtransprobs, node.children[1], data, columns, d)
            #cstore[d] = logtransprobs[node.children[2].data.nodeindex,a,d] + felsenstein(cache, logtransprobs, node.children[2], data, columns, d)
            bsum = CommonUtils.logsumexp(bsum, logtransprobs[node.children[1].data.nodeindex,a,d] + cache[node.children[1].data.nodeindex,d])
            csum = CommonUtils.logsumexp(csum, logtransprobs[node.children[2].data.nodeindex,a,d] + cache[node.children[2].data.nodeindex,d])
          end
          v[a] = bsum+csum
        end
        subcolumnloglikcache[node.data.nodeindex][key] = v
        cache[node.data.nodeindex,:] = v
      end
    end
  end
  return cache[1,:]
end

function felsenstein(dataset::Dataset, subcolumnloglikcache::Array{Dict{Int,Array{Float64,1}},1}, subcolumnrefs::Array{Int,2}, logtransprobs::Array{Float64,3}, node::TreeNode, logdata::Array{Float64,3}, columns::Array{Int,1}, alphabet::Int)
  if isleafnode(node)
    return logstate(logdata, node.data.seqindex,columns)
  else
    key = getkey(dataset.numnodes, node.data.nodeindex, subcolumnrefs, dataset.numcols, columns)
    #key = (subcolumnrefs[node.data.nodeindex,columns[1]]-1)*dataset.numcols + subcolumnrefs[node.data.nodeindex,columns[2]]
    if haskey(subcolumnloglikcache[node.data.nodeindex],key)
      return subcolumnloglikcache[node.data.nodeindex][key]
    else
      v1 = felsenstein(dataset, subcolumnloglikcache, subcolumnrefs, logtransprobs, node.children[1], logdata, columns, alphabet)
      v2 = felsenstein(dataset, subcolumnloglikcache, subcolumnrefs, logtransprobs, node.children[2], logdata, columns, alphabet)
      v = zeros(Float64, alphabet)
      for a=1:alphabet
        bsum = -Inf
        csum = -Inf
        for d=1:alphabet
          bsum = CommonUtils.logsumexp(bsum,logtransprobs[node.children[1].data.nodeindex,a,d] + v1[d])
          csum = CommonUtils.logsumexp(csum,logtransprobs[node.children[2].data.nodeindex,a,d] + v2[d])
        end
        v[a] = bsum+csum
      end
      subcolumnloglikcache[node.data.nodeindex][key] = v
      return v
    end
  end
end


function felsensteinfast(dataset::Dataset, cache::Array{Float64,2}, fastcache::FastCache, subcolumnrefs::Array{Int,2}, transprobs::Array{Array{Float64,2},1}, node::TreeNode, data::Array{Float64,3}, columns::Array{Int,1}, alphabet::Int)
  nodeindex = node.data.nodeindex
  if isleafnode(node)
    return state(data, node.data.seqindex,columns)
  elseif cache[nodeindex,1] != 0.0
    return cache[nodeindex,:]
  else
    key = getkey(dataset.numnodes, nodeindex, subcolumnrefs, dataset.numcols, columns)
    index = getcacheindex(fastcache,key)

    if index != 0
      for z=1:alphabet+1
        cache[nodeindex,z] = fastcache.cache[index,z]
      end
      return fastcache.cache[index,:]
    else
      v1 = felsensteinfast(dataset, cache, fastcache, subcolumnrefs, transprobs, node.children[1], data, columns, alphabet)
      v2 = felsensteinfast(dataset, cache, fastcache, subcolumnrefs, transprobs, node.children[2], data, columns, alphabet)

      leftchildindex = node.children[1].data.nodeindex
      rightchildindex = node.children[2].data.nodeindex
      #transmatleft = transprobs[leftchildindex,:,:]
      #transmatright = transprobs[rightchildindex,:,:]
      #=
      m = 0.0
      for a=1:alphabet
        bsum = 0.0
        csum = 0.0
        transvecleft = transprobs[leftchildindex,a,:]
        transvecright = transprobs[rightchildindex,a,:]
        #bsum = dot(transprobs[leftchildindex,a,:], v1[1:alphabet])
        #csum = dot(transprobs[rightchildindex,a,:], v2[1:alphabet])

        for d=1:alphabet
          bsum += transvecleft[d]*v1[d]
          csum += transvecright[d]*v2[d]
          #bsum += transprobs[leftchildindex,a,d]*v1[d]
          #csum += transprobs[rightchildindex,a,d]*v2[d]
        end
        cache[nodeindex,a] = bsum*csum
        if cache[nodeindex,a] > m
          m = cache[nodeindex,a]
        end
      end=#


      cache[nodeindex,1:alphabet] = (transprobs[leftchildindex]*v1[1:alphabet]).*(transprobs[rightchildindex]*v2[1:alphabet])
      m = 0.0
      for a=1:alphabet
        if cache[nodeindex,a] > m
          m = cache[nodeindex,a]
        end
      end
      if m < 1e-20
        for a=1:alphabet
          cache[nodeindex,a] /= m
        end
        cache[nodeindex,alphabet+1] = log(m)+v1[alphabet+1]+v2[alphabet+1]
      else
        cache[nodeindex,alphabet+1] = v1[alphabet+1]+v2[alphabet+1]
      end

      index = getandkillcacheindex(fastcache, key)
      fastcache.keys[index] = key
      for z=1:alphabet+1
        fastcache.cache[index,z] = cache[nodeindex,z]
      end
      return cache[nodeindex,:]
    end
  end
end

function felsensteinstack(dataset::Dataset, nodelist::Array{TreeNode,1}, cache::Array{Array{Float64,1},1}, logm::Array{Float64,1}, fastcache::FastCache, subcolumnrefs::Array{Int,2}, transprobs::Array{Array{Float64,2},1}, node::TreeNode, data::Array{Float64,3}, columns::Array{Int,1}, alphabet::Int)
  stack = Int[1]
  while length(stack) > 0
    nodeindex = stack[end]
    node = nodelist[nodeindex]
    if cache[nodeindex][1] != -Inf
      pop!(stack)
    elseif isleafnode(node)
      v = state(data, node.data.seqindex,columns)
      for z=1:alphabet
        cache[nodeindex][z] = v[z]
      end
      pop!(stack)
    else
      key = getkey(dataset.numnodes, node.data.nodeindex, subcolumnrefs, dataset.numcols, columns)
      index = getcacheindex(fastcache,key)
      if index != 0
        for z=1:alphabet
          cache[nodeindex][z] = fastcache.cache[index,z]
        end
        logm[nodeindex] = fastcache.logm[index]
        pop!(stack)
      else
        leftchildindex = node.children[1].data.nodeindex
        rightchildindex = node.children[2].data.nodeindex
        cont = true
        if cache[leftchildindex][1] == -Inf
          push!(stack, leftchildindex)
          cont = false
        end
        if cache[rightchildindex][1] == -Inf
          push!(stack, rightchildindex)
          cont = false
        end

        if cont
          cache[nodeindex] = (transprobs[leftchildindex]*cache[leftchildindex]).*(transprobs[rightchildindex]*cache[rightchildindex])
          #=
          m = 0.0
          for a=1:alphabet
            if cache[nodeindex,a] > m
              m = cache[nodeindex,a]
            end
          end=#
          m = maximum(cache[nodeindex])
          if m < 1e-20
            cache[nodeindex] /= m
            #=
            for a=1:alphabet
              cache[nodeindex,a] /= m
            end=#
            logm[nodeindex] = log(m) + logm[leftchildindex] + logm[rightchildindex]
            #cache[nodeindex,alphabet+1] = log(m)+cache[leftchildindex,alphabet+1]+cache[rightchildindex,alphabet+1]
          else
            #cache[nodeindex,alphabet+1] = cache[leftchildindex,alphabet+1]+cache[rightchildindex,alphabet+1]
            logm[nodeindex] = logm[leftchildindex] + logm[rightchildindex]
          end

          index = getandkillcacheindex(fastcache, key)
          fastcache.keys[index] = key
          for z=1:alphabet
            fastcache.cache[index,z] = cache[nodeindex][z]
          end
          fastcache.logm[index] = logm[nodeindex]
          pop!(stack)
        end

      end
    end
  end

  return cache[1]
end


function felsensteinfastiterative(dataset::Dataset, revnodelist::Array{TreeNode,1}, cache::Array{Float64,2}, fastcache::FastCache, subcolumnrefs::Array{Int,2}, transprobs::Array{Array{Float64,2},1}, node::TreeNode, data::Array{Float64,3}, columns::Array{Int,1}, alphabet::Int)
  for node in revnodelist
    nodeindex = node.data.nodeindex
    if isleafnode(node)
      cache[nodeindex,:] = state(data, node.data.seqindex,columns)
    elseif cache[nodeindex,1] != 0.0

    else
      key = getkey(dataset.numnodes, node.data.nodeindex, subcolumnrefs, dataset.numcols, columns)
      index = getcacheindex(fastcache,key)

      if index > 0
        for z=1:alphabet+1
          cache[nodeindex,z] = fastcache.cache[index,z]
        end
      else
        leftchildindex = node.children[1].data.nodeindex
        rightchildindex = node.children[2].data.nodeindex
        #v1 =
        #v2 = cache[rightchildindex,:]
        #lefttransprobs = transprobs[leftchildindex]
        #righttransprobs = transprobs[rightchildindex]
        #=
        m = 0.0
        for a=1:alphabet
          bsum = 0.0
          csum = 0.0
          for d=1:alphabet
            bsum += lefttransprobs[a,d]*v1[d]
            csum += righttransprobs[a,d]*v2[d]
          end
          cache[nodeindex,a] = bsum*csum
          if cache[nodeindex,a] > m
            m = cache[nodeindex,a]
          end
        end=#

        cache[nodeindex,1:alphabet] = (transprobs[leftchildindex]*(cache[leftchildindex,1:alphabet]')).*(transprobs[rightchildindex]*(cache[rightchildindex,1:alphabet]'))
        m = 0.0
        for a=1:alphabet
          if cache[nodeindex,a] > m
            m = cache[nodeindex,a]
          end
        end
        if m < 1e-20
          for a=1:alphabet
            cache[nodeindex,a] /= m
          end
          cache[nodeindex,alphabet+1] = log(m)+cache[leftchildindex,alphabet+1]+cache[rightchildindex,alphabet+1]
        else
          cache[nodeindex,alphabet+1] = cache[leftchildindex,alphabet+1]+cache[rightchildindex,alphabet+1]
        end

        index = getandkillcacheindex(fastcache, key)
        fastcache.keys[index] = key
        for z=1:alphabet+1
          fastcache.cache[index,z] = cache[nodeindex,z]
        end
      end
    end
  end

  return cache[1,:]
end

function computestructureintegratedlikelihood(dataset::Dataset, params::ModelParameters, maxbasepairdistance::Int, samplebranchlengths::Bool, fixLambdaGT::Bool, integratesiterates::Bool, usecuda::Bool=true)
  tempGT = params.lambdaGT
  if fixLambdaGT
    params.lambdaGT = 1.0
  end

  #tic()
  unpairedlogprobs = nothing
  if integratesiterates
    unpairedlogprobs = computeunpairedlikelihoods(dataset, params)
  else
    unpairedlogprobs = computeunpairedlikelihoods(dataset, params, params.states)
  end
  pairedlogprobs, ret = coevolutionall(dataset,params,true,false,integratesiterates,maxbasepairdistance, usecuda)
  #elapsed = #toc()
  #println("Coevolution ", elapsed)
  #tic()
  Z = computeinsideKH99(unpairedlogprobs, pairedlogprobs, 1.0, true, usecuda)
  #elapsed = #toc()
  #println("Inside ", elapsed)
  params.logprior = computeprior(params, samplebranchlengths, fixLambdaGT)
  params.lambdaGT = tempGT
  return params.logprior + Z
end

function computetotallikelihood(rng::AbstractRNG, dataset::Dataset,params::ModelParameters, paired::Array{Int,1},samplebranchlengths::Bool,integratesiterates::Bool,integratestructure::Bool, M::Float64, useprior::Bool=true,maxbasepairdistance::Int=500, usecuda::Bool=true)
  if integratestructure
    if M == 1.0
      params.logposterior = computestructureintegratedlikelihood(dataset, params, maxbasepairdistance, samplebranchlengths, false, integratesiterates, usecuda)
      params.loglikelihood = params.logposterior - params.logprior
      return params.logposterior
    elseif M == 0.0
      params.logposterior = computestructureintegratedlikelihood(dataset, params, maxbasepairdistance, samplebranchlengths, true, integratesiterates, usecuda)
      params.loglikelihood = params.logposterior - params.logprior
      return params.logposterior
    else
      q0 = computestructureintegratedlikelihood(dataset, params, maxbasepairdistance, samplebranchlengths, true, integratesiterates, usecuda)
      q1 = computestructureintegratedlikelihood(dataset, params, maxbasepairdistance, samplebranchlengths, false, integratesiterates, usecuda)
      return q0*(1.0 - M) + q1*M
    end
  else
    if M == 1.0
      return computemodellikelihood(rng, dataset, params, paired, samplebranchlengths, integratesiterates, useprior,false, usecuda)
    elseif M == 0.0
      tempGT = params.lambdaGT
      params.lambdaGT = 1.0
      ll = computemodellikelihood(rng, dataset, params, paired, samplebranchlengths, integratesiterates, useprior,true, usecuda)
      params.lambdaGT = tempGT
      return ll
    else
      q1 = computemodellikelihood(rng, dataset, params, paired, samplebranchlengths, integratesiterates, useprior, false, usecuda)
      tempGT = params.lambdaGT
      params.lambdaGT = 1.0
      q0 = computemodellikelihood(rng, dataset, params, paired, samplebranchlengths, integratesiterates, useprior, true, usecuda)
      params.lambdaGT = tempGT
      return q0*(1.0 - M) + q1*M
    end
  end
end

function getmodelswitchlikelihoods(rng::AbstractRNG, dataset::Dataset,params::ModelParameters, paired::Array{Int,1},samplebranchlengths::Bool,integratesiterates::Bool,integratestructure::Bool, M::Float64, useprior::Bool=true, maxbasepairdistance::Int=500, usecuda::Bool=true)
    if integratestructure
      q0 = computestructureintegratedlikelihood(dataset, params, maxbasepairdistance, samplebranchlengths, true,integratesiterates,usecuda)
      q1 = computestructureintegratedlikelihood(dataset, params, maxbasepairdistance, samplebranchlengths, false,integratesiterates,usecuda)
      CommonUtils.logsumexp(log(1.0 - M)+q0, log(M)+q1), q0, q1, q1-q0
    else
      q1 = computemodellikelihood(rng, dataset, params, paired, samplebranchlengths, integratesiterates, useprior, false,usecuda)
      tempGT = params.lambdaGT
      params.lambdaGT = 1.0
      q0 = computemodellikelihood(rng, dataset, params, paired, samplebranchlengths, integratesiterates, useprior, true,usecuda)
      params.lambdaGT = tempGT
      return CommonUtils.logsumexp(log(1.0 - M)+q0, log(M)+q1), q0, q1, q1-q0
    end
end


gammaprior = Gamma(2.0, 0.5)
gtprior = Normal(0.0,4.0)
betaprior = Beta(2.0,2.0)
function computeprior(params::ModelParameters, samplebranchlengths::Bool, fixLambdaGT::Bool)
  priorll = 0.0

  priorll += -(params.lambdaGC-1.0)*0.1
  priorll += -(params.lambdaAT-1.0)*0.1
  if !fixLambdaGT
    priorll += -(params.lambdaGT-1.0)*0.1
  end
  #=
  if !fixLambdaGT
    priorll += logpdf(gtprior, log(params.lambdaGT))
  end=#
  #priorll += logpdf(gammaprior, params.lambdaGC-1.0)
  #priorll += logpdf(gammaprior, params.lambdaAT-1.0)
  #priorll += logpdf(gammaprior, params.lambdaGT-1.0)
  priorll += -params.q1*0.1
  priorll += -params.q2*0.1
  priorll += -params.q3*0.1
  priorll += -params.q4*0.1
  priorll += -params.q5*0.1
  priorll += -params.q6*0.1
  if samplebranchlengths
    for b=2:length(params.branchlengths)
      priorll += -params.branchlengths[b]*0.2
    end
  end
  priorll += -0.1*params.lambdaGammaShapeGC
  priorll += -0.1*params.lambdaGammaScaleGC
  priorll += -0.1*params.lambdaGammaShapeAT
  priorll += -0.1*params.lambdaGammaScaleAT
  priorll += -0.1*params.lambdaGammaShapeGT
  priorll += -0.1*params.lambdaGammaScaleGT
  priorll += -0.1*params.siteGammaShape
  priorll += logpdf(betaprior, params.lambdazeroweight)
  return priorll
end


function computemodellikelihood(rng::AbstractRNG, dataset::Dataset,params::ModelParameters, paired::Array{Int,1},samplebranchlengths::Bool,integratesiterates::Bool, useprior::Bool=true, fixLambdaGT::Bool=false, usecuda::Bool=true)

  priorll = 0.0
  if useprior
     priorll = computeprior(params, samplebranchlengths, fixLambdaGT)
  end
  #=
  for z=1:length()
    Gamma(2.0, 0.5)
  end=#
  #=
  for lambdarate in params.lambdarates
    #priorll += -(lambdarate-1.0)*1.0
    write(posteriorlog, "iter\t")
  end
  priorll += logpdf(Dirichlet(ones(Float64,length(params.lambdaweights))*c), params.lambdaweights)
  =#

  ##tic()
  #pairedllref = coevolution(dataset, params, paired)
  #pairedllref = @spawn coevolution(dataset, params, paired)

  sites = Int[]
  for i=1:length(paired)
    if paired[i] == 0
      push!(sites,i)
    end
  end
  #unpairedlogprobsref = computelikelihoods(dataset,params,sites)
  #unpairedlogprobsref = @spawn computelikelihoods(dataset,params,sites)
  #elapsed = #toc()
  #println("singlell=", elapsed)
  #pairedll = fetch(pairedllref)
  #unpairedlogprobs = fetch(unpairedlogprobsref)
  #pairedll =  computedpairedlikelihoods(mcmcCache, dataset, params, paired)
  #unpairedlogprobs = computeunpairedlikelihoods(dataset,params,sites)
  museparams = getmusespecificparamsarray(params)

  unpairedll = 0.0
  pairedll = 0.0
  if integratesiterates

    if usecuda
      #tic()
      unpairedll = sum(computeunpairedlikelihoods(dataset, params, museparams, params.states, paired, 1.0, false, rng, true))
      pairedll = felsensteincudapaired(dataset, paired, params, museparams)
      #GC.gc()
      #elapsed = toq()
      #logtime("structurell_cuda", elapsed, dataset.numcols)
    else
      #tic()
      unpairedloglikelihoods,pairedloglikelihoods,params.states,museconditionals = computeuncachedlikelihood(dataset, params, museparams, params.states, paired, 1.0, false, rng, true)
      for i=1:dataset.numcols
        if paired[i] == 0
          unpairedll += unpairedloglikelihoods[i]
        elseif paired[i] > i
          pairedll +=  pairedloglikelihoods[i]
        end
      end
      #elapsed = toq()
      #logtime("structurell_cpu", elapsed, dataset.numcols)
    end

  else
    unpairedloglikelihoods,pairedloglikelihoods = computeuncachedlikelihood(dataset, params, museparams, params.states, paired)
    for i=1:dataset.numcols
      if paired[i] == 0
        unpairedll += unpairedloglikelihoods[i]
      elseif paired[i] > i
        pairedll +=  pairedloglikelihoods[i]
      end
    end
  end
  ##println(pairedloglikelihoods)
  #println("£££££",unpairedll,"\t",pairedll)

  return priorll+unpairedll+pairedll
end

function samplestatespairedhelper(rng::AbstractRNG, dataset::Dataset, params::ModelParameters, paired::Array{Int,1}, siteCat1::Int, siteCat2::Int, parallel::Bool=true)
  states = zeros(Int,dataset.numcols)
  for i=1:length(paired)
    if paired[i] > i
      states[i] = siteCat1
      states[paired[i]] = siteCat2
    end
  end
  return computedpairedlikelihoodscats(dataset, params, states, paired)[1]
end

function samplestatessinglehelper(rng::AbstractRNG, dataset::Dataset, params::ModelParameters, sites::Array{Int,1}, siteCat::Int, parallel::Bool=true)
  states = ones(Int,dataset.numcols)*siteCat
  return computeunpairedlikelihoodscats(dataset, params, states, sites)
end

function samplestates(rng::AbstractRNG, dataset::Dataset, params::ModelParameters, paired::Array{Int,1}, B::Float64=1.0, parallel::Bool=true)
  museparams = getmusespecificparamsarray(params)
  unpairedloglikelihoods,pairedloglikelihoods,params.states,museconditionals = computeuncachedlikelihood(dataset, params, museparams, params.states, paired, B, false, rng)
  #=
  pairedchoice = zeros(Float64, dataset.numcols, params.siteCats*params.siteCats)
  singlechoice = zeros(Float64, dataset.numcols, params.siteCats)
  sites = Int[]
  for i=1:length(paired)
    if paired[i] == 0
      push!(sites,i)
    end
  end

  singlerefs = []
  pairrefs = []
  for siteCat1=1:params.siteCats
    for siteCat2=1:params.siteCats
      ref = @spawn samplestatespairedhelper(rng, dataset, params, paired, siteCat1, siteCat2, parallel)
      push!(pairrefs, ref)
    end
    ref =  @spawn samplestatessinglehelper(rng, dataset, params, sites, siteCat1, parallel)
    push!(singlerefs,ref)
  end

  for siteCat1=1:params.siteCats
    logunpairedweight = log(params.siteWeights[siteCat1])
    for siteCat2=1:params.siteCats
      logpairedweight = log(params.siteWeights[siteCat1]) + log(params.siteWeights[siteCat2])
      index = (siteCat1-1)*params.siteCats+siteCat2
      pairedloglikelihoods = fetch(pairrefs[index])
      for i=1:length(paired)
        if paired[i] > i
          pairedchoice[i,index] = (logpairedweight + pairedloglikelihoods[i])*B
        end
      end
    end

    unpairedloglikelihoods = fetch(singlerefs[siteCat1])
    for i=1:length(paired)
      if paired[i] == 0
        singlechoice[i,siteCat1] = (logunpairedweight+unpairedloglikelihoods[i])*B
      end
    end
  end

  for i=1:length(paired)
    if paired[i] > i
      total = -Inf
      for siteCat1=1:params.siteCats
        for siteCat2=1:params.siteCats
          index = (siteCat1-1)*params.siteCats+siteCat2
          total = CommonUtils.logsumexp(total, pairedchoice[i,index])
        end
      end
      pairedchoice[i,:] = exp(pairedchoice[i,:]-total)
      index = CommonUtils.sample(rng,pairedchoice[i,:])
      params.states[i] = div(index-1,params.siteCats) + 1
      params.states[paired[i]] = ((index-1) % params.siteCats) + 1
    else
      total = -Inf
      for siteCat1=1:params.siteCats
        total = CommonUtils.logsumexp(total, singlechoice[i,siteCat1])
      end
      singlechoice[i,:] = exp(singlechoice[i,:]-total)
      params.states[i] = CommonUtils.sample(rng,singlechoice[i,:])
    end
  end=#

end

include("FelsensteinCUDA.jl")

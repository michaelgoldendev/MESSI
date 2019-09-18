include("RNATools.jl")

MODE_FIXED_STRUCTURE = 0
MODE_IO_SAMPLE = 1
MODE_VIENNA_SAMPLE = 2

using JSON
using CommonUtils
#using HDF5

function getparamsvector2(params::ModelParameters)
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
    p[12] = params.freqs[1]
    p[13] = params.freqs[2]
    p[14] = params.freqs[3]
    p[15] = params.lambdazeroweight
    #=
    p[16] = params.lambdaGammaScaleGC
    p[17] = params.lambdaGammaShapeAT
    p[18] = params.lambdaGammaScaleAT
    p[19] = params.lambdaGammaShapeGT
    p[20] = params.lambdaGammaScaleGT=#
    #p[15] = params.freqs[4]
    return p
end

maximumZ2 = -Inf
function savemaximum2(Z::Float64, maxparams::Array{Float64,1}, maxfile, override::Bool=false)
  global maximumZ2
  maximumparams = maxparams

  if isfile(maxfile)
    jsondict = JSON.parse(open(maxfile))
    if Z > jsondict["Z"] || override
      if Z > maximumZ2  || override
        maximumZ2 = Z
        jsondict["Z"] = maximumZ2
        jsondict["maxparams"] = maxparams
        maximumparams = jsondict["maxparams"]
        out = open(maxfile,"w")
        ret = replace(JSON.json(jsondict),",\"" => ",\n\"")
        ret = replace(ret, "],[" => "],\n[")
        ret = replace(ret, "{" => "{\n")
        ret = replace(ret, "}" => "\n}")
        write(out,ret)
        close(out)
      end
    end
  elseif Z > maximumZ2 || override
    maximumZ2 = Z
    jsondict = Dict()
    jsondict["Z"] = maximumZ2
    jsondict["maxparams"] = maxparams
    maximumparams = jsondict["maxparams"]
    out = open(maxfile,"w")
    ret = replace(JSON.json(jsondict),",\"" => ",\n\"")
    ret = replace(ret, "],[" => "],\n[")
    ret = replace(ret, "{" => "{\n")
    ret = replace(ret, "}" => "\n}")
    write(out,ret)
    close(out)
  end

  return maximumZ2, maximumparams
end

function calculatemean(mcmclogfile, numcols::Int)
  fin = open(mcmclogfile,"r")
  header = true
  numlines = 1
  for line in readlines(fin)
    if header
      spl = split(line)
      header = false
    else
      spl = split(line)
      if length(spl) > 1
        numlines += 1
      end
    end
  end
  close(fin)

  fin = open(mcmclogfile,"r")
  header = true
  count = 0.0
  lineindex = 1
  totals = zeros(Float64,numcols)
  for line in readlines(fin)
    if header
      header = false
    else
      if lineindex > numlines/2
        spl = split(line)
        if length(spl) > 1
          for j=2:length(spl)
            totals[j-1] += parse(Float64, spl[j])
          end
          count += 1.0
        end
      end
      lineindex += 1
    end
  end
  close(fin)
  println(count,"\t",numlines)
  return totals / max(1.0, count)
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
  if len > 2
    println("lambdaGC > lambdaAT: ", mean([lambdaGC > lambdaAT ? 1.0 : 0.0 for (lambdaAT,lambdaGC) in zip(params["lambdaAT"][div(len,2):end],params["lambdaGC"][div(len,2):end])]))
    println("lambdaAT > lambdaGT: ", mean([lambdaAT > lambdaGT ? 1.0 : 0.0 for (lambdaGT,lambdaAT) in zip(params["lambdaGT"][div(len,2):end],params["lambdaAT"][div(len,2):end])]))
    println("lambdaGC > lambdaGT: ", mean([lambdaGC > lambdaGT ? 1.0 : 0.0 for (lambdaGT,lambdaGC) in zip(params["lambdaGT"][div(len,2):end],params["lambdaGC"][div(len,2):end])]))
  end
  #=
  startindex = div(len,2)
  len2 = length(params["lambdaGC"][startindex:end])
  data = zeros(len2,3)
  for i=1:len2
    data[i,1] = params["lambdaGC"][startindex+i-1]
    data[i,2] = params["lambdaAT"][startindex+i-1]
    data[i,3] = params["lambdaGT"][startindex+i-1]
  end
  return keys,params,A=#
end



spawnedindices = Int[]
spawned = Any[]
function getnext(spawnindex::Int, rng::AbstractRNG, dataset::Dataset, paired::Array{Int,1}, currentparams::ModelParameters, moveweights::Array{Float64,1}, samplebranchlengths, lambdaMv::MvNormal, QMv::MvNormal)
  #println("S", spawnedindices,"\t", spawnindex)
  while length(spawnedindices) > 0 && spawnedindices[1] != spawnindex
    #println("A")
    popfirst!(spawnedindices)
    popfirst!(spawned)
    #println("B")
  end

  #println("R", spawnedindices)
  if length(spawnedindices) > 0
    if isready(spawned[1])
      popfirst!(spawnedindices)
      return fetch(popfirst!(spawned))
    end
  end
    #println("E", spawnedindices)
  while length(spawnedindices) < 4
    rng2 =  MersenneTwister(rand(1:typemax(Int)))
    #ref = @spawn proposeparams(rng2, dataset::Dataset, paired::Array{Int,1}, currentparams, moveweights,samplebranchlengths, lambdaMv, QMv)
    ref = proposeparams(rng2, dataset::Dataset, paired::Array{Int,1}, currentparams, moveweights,samplebranchlengths, lambdaMv, QMv)
    push!(spawnedindices,spawnindex)
    push!(spawned,ref)
  end
  #println("Z",spawnedindices,"\t",spawnindex, "\t", length(spawned))
  popfirst!(spawnedindices)
  #return fetch(popfirst!(spawned))
  return popfirst!(spawned)

end

function getmatrix(matrixfile, m::Int)
    fin = open(matrixfile,"r")
    count = 0
    matrixlines = []
    for line in readlines(fin)
      if line[1] == '>'
        count += 1
      elseif m == count
        spl1 = split(strip(strip(line), ['[',']']), ";")
        for line2 in spl1
          push!(matrixlines, [parse(Float64,v) for v in split(strip(line2))])
        end
      end
    end
    close(fin)

    if m > 0 && m <= count
      dim1 = length(matrixlines)
      dim2 = length(matrixlines[1])
      matrix = zeros(Float64, dim1, dim2)
      for i=1:dim1
        for j=1:dim2
          matrix[i,j] = matrixlines[i][j]
        end
      end
      return matrix,count
    end
    return nothing,count
end


include("StructureVisualisation.jl")

mutable struct Chain
  iter::Int
  id::Int
  B::Float64
  rng::AbstractRNG
  currentll::Float64
  proposedll::Float64
  currentparams::ModelParameters
  proposedparams::ModelParameters
  paired::Array{Int,1}
  inside::Array{Float64,3}
  pairedlogprobs::Array{Float64,2}
  unpairedlogprobs::Array{Float64,1}
  logger::AcceptanceLogger
  timings::Array{Tuple{Int,Any},1}
  lasthash::UInt64
  hashindex::Int
  maxinside::Array{Float64,3}
  maxpairedlogprobs::Array{Float64,2}
  maxunpairedlogprobs::Array{Float64,1}
  maxparams::ModelParameters
  pairedlogprior::Float64

  function Chain(id::Int, B::Float64, rng::AbstractRNG, currentll::Float64, proposedll::Float64, currentparams::ModelParameters, paired::Array{Int,1}, inside::Array{Float64,3}, pairedlogprobs::Array{Float64,2}, unpairedlogprobs::Array{Float64,1}, maxinside::Array{Float64,3}, maxpairedlogprobs::Array{Float64,2}, maxunpairedlogprobs::Array{Float64,1}, maxparams::ModelParameters)
    timings = Tuple{Int,Any}[(-1,time())]
    new(0, id, B, rng, currentll, proposedll, currentparams, deepcopy(currentparams), paired, inside, pairedlogprobs, unpairedlogprobs, AcceptanceLogger(), timings, 0, 0, maxinside, maxpairedlogprobs, maxunpairedlogprobs, maxparams, -Inf)
  end
end

function onedarray(a, t::Type=Float64)
  ret = zeros(Float64, length(a))
  for i=1:size(ret,1)
    ret[i] = a[i]
  end
  return ret
end

function twodarray(matrix, t::Type=Float64)
  ret = zeros(t, length(matrix), length(matrix[1]))
  for i=1:size(ret,1)
    for j=1:size(ret,2)
      ret[i,j] = matrix[i][j]
    end
  end
  return ret
end

function sampletruncated(rng::AbstractRNG, mu::Float64, sigma::Float64, left::Float64, right::Float64)
  s = mu + randn(rng)*sigma
  while !(left < s < right)
    s = mu + randn(rng)*sigma
  end
  return s
end

function isvalid(params::ModelParameters)
  if params.lambdazeroweight < 0.0 || params.lambdazeroweight > 1.0 || params.q1 < 0.0 || params.q2 < 0.0 || params.q3 < 0.0 || params.q4 < 0.0 || params.q5 < 0.0 || params.q6 < 0.0  || params.siteGammaShape < 0.0 || params.lambdaGammaShapeGC < 1.0 || params.lambdaGammaScaleGC < 0.0  || params.lambdaGammaShapeAT < 1.0 || params.lambdaGammaScaleAT < 0.0  || params.lambdaGammaShapeGT < 1.0 || params.lambdaGammaScaleGT < 0.0
  #if chain.proposedparams.lambdaGC < 1.0 || chain.proposedparams.lambdaAT < 1.0 || chain.proposedparams.lambdaGT < 1.0 || chain.proposedparams.q1 < 0.0 || chain.proposedparams.q2 < 0.0 || chain.proposedparams.q3 < 0.0 || chain.proposedparams.q4 < 0.0 || chain.proposedparams.q5 < 0.0 || chain.proposedparams.q6 < 0.0
    return false
  end

  for b=1:length(params.branchlengths)
    if params.branchlengths[b] < 0.0
      return false
    end
  end

  return true
end
#
function gettuningvector(params::ModelParameters)
    p = zeros(Float64,14)
    p[1] = params.q1
    p[2] = params.q2
    p[3] = params.q3
    p[4] = params.q4
    p[5] = params.q5
    p[6] = params.q6
    p[7] = params.siteGammaShape
    p[8] = params.lambdazeroweight

    p[9] = params.lambdaGC
    p[10] = params.lambdaGammaScaleGC
    p[11] = params.lambdaAT
    p[12] = params.lambdaGammaScaleGC
    p[13] = params.lambdaGT
    p[14] = params.lambdaGammaScaleGC

    return p
end



function settuningparams(params::Array{Float64,1}, currentparams::ModelParameters)
  currentparams.q1 = params[1]
  currentparams.q2 = params[2]
  currentparams.q3 = params[3]
  currentparams.q4 = params[4]
  currentparams.q5 = params[5]
  currentparams.q6 = params[6]
  currentparams.siteGammaShape = params[7]
  currentparams.siteGammaScale = 1.0 / currentparams.siteGammaShape
  currentparams.lambdazeroweight = params[8]
  currentparams.lambdaGC = params[9]
  currentparams.lambdaGammaShapeGC = params[9]
  currentparams.lambdaGammaScaleGC = params[10]

  currentparams.lambdaAT = params[11]
  currentparams.lambdaGammaShapeAT = params[11]
  currentparams.lambdaGammaScaleAT = params[10]

  currentparams.lambdaGT = params[13]
  currentparams.lambdaGammaShapeGT = params[13]
  currentparams.lambdaGammaScaleGT = params[10]
  #=
  currentparams.lambdaGammaShapeGC = params[9]
  currentparams.lambdaGammaScaleGC = params[10]
  currentparams.lambdaGammaShapeAT = params[11]
  currentparams.lambdaGammaScaleAT = params[12]
  currentparams.lambdaGammaShapeGT = params[13]
  currentparams.lambdaGammaScaleGT = params[14]
  =#

  if isvalid(currentparams)
    #currentparams.lambdaweightsGC, currentparams.lambdaratesGC = discretizegamma2(currentparams.lambdazeroweight, currentparams.lambdaGammaShapeGC, currentparams.lambdaGammaScaleGC, currentparams.lambdaCats)
    #currentparams.lambdaweightsAT, currentparams.lambdaratesAT = discretizegamma2(currentparams.lambdazeroweight, currentparams.lambdaGammaShapeAT, currentparams.lambdaGammaScaleAT, currentparams.lambdaCats)
    #currentparams.lambdaweightsGT, currentparams.lambdaratesGT = discretizegamma2(currentparams.lambdazeroweight, currentparams.lambdaGammaShapeGT, currentparams.lambdaGammaScaleGT, currentparams.lambdaCats)
    currentparams.lambdaweightsGC, currentparams.lambdaratesGC, currentparams.lambdaweightsAT, currentparams.lambdaratesAT, currentparams.lambdaweightsGT, currentparams.lambdaratesGT = discretizegamma3(currentparams.lambdazeroweight, currentparams.lambdaGammaShapeGC, currentparams.lambdaGammaShapeAT, currentparams.lambdaGammaShapeGT, currentparams.lambdaGammaScaleGC, currentparams.lambdaCats)
    currentparams.siteRates = discretizegamma(currentparams.siteGammaShape, currentparams.siteGammaScale, currentparams.siteCats)
  else
    #println("INVALID", currentparams)
  end
  return currentparams
end

function maskarray(arr::Array{Float64,1}, mask::Array{Int,1})
  ret = copy(arr)
  for i=1:length(mask)
    if mask[i] == 0
      ret[i] = 0.0
    end
  end
  return ret
end

function hcathelper(tuningvectors::Array{Array{Float64,1},1})
  len = length(tuningvectors[1])
  mat = zeros(Float64, len, length(tuningvectors))
  for i=1:length(tuningvectors)
    for j=1:len
      mat[j,i] = tuningvectors[i][j]
    end
  end
  return mat
end

function gettunedparams(current::Array{Float64,1}, tuningvectors::Array{Array{Float64,1},1}, maskin::Array{Int,1}, covarScale::Float64=1.0, maxTuningIter::Int=1000, covmatin=nothing, C::Float64=0.005)
  covmat = covmatin
  mask = copy(maskin)
  d = length(current)
  d2 = sum(mask)
  e = eye(d)
  r = rand()
  if r < 0.25
    e *= 0.00025
  elseif r < 0.5
    e *= 0.0025
  elseif r < 0.75
    e *= 0.0075
  else
    e *= 0.025
  end
  len = length(tuningvectors)
  if  len < 16*d
    e[8,8] *= 0.25
    #eyemat = (e)/Float64(d)
    #eyemat = [0.0714252 0.0380393 0.0179428 0.00612939 -0.000207974 0.00307991 0.00304808 0.00698059 0.00248841 -0.0121903 -0.000179019 -0.000714703; 0.0380393 0.0515734 0.0154519 0.000574956 -0.00108632 0.00374097 0.00308298 0.00227931 0.00283558 -0.0102061 0.000335268 -1.20491e-5; 0.0179428 0.0154519 0.00900712 -0.000671783 0.00013602 0.000293199 0.00128999 -0.00506506 0.000271419 -0.00439561 0.000107318 -4.93604e-5; 0.00612939 0.000574956 -0.000671783 0.00808231 0.000861339 0.00178185 0.0026731 0.0148671 0.00275768 -0.00047444 -0.000393101 -5.4529e-5; -0.000207974 -0.00108632 0.00013602 0.000861339 0.008077 -0.000420196 0.00267563 -0.00131275 0.000997578 0.00012277 -0.000291644 -0.000223464; 0.00307991 0.00374097 0.000293199 0.00178185 -0.000420196 0.0063449 0.00125321 0.015334 0.00211563 -0.00085548 -0.000120409 9.03542e-6; 0.00304808 0.00308298 0.00128999 0.0026731 0.00267563 0.00125321 0.011305 0.0131275 0.00069293 -0.00088149 -0.000123356 -0.000140382; 0.00698059 0.00227931 -0.00506506 0.0148671 -0.00131275 0.015334 0.0131275 0.0895556 0.00737556 -0.000463992 -0.000569329 0.000248975; 0.00248841 0.00283558 0.000271419 0.00275768 0.000997578 0.00211563 0.00069293 0.00737556 0.00471388 -0.00104256 -0.000271204 0.000110928; -0.0121903 -0.0102061 -0.00439561 -0.00047444 0.00012277 -0.00085548 -0.00088149 -0.000463992 -0.00104256 0.00360694 -2.36933e-5 6.73943e-5; -0.000179019 0.000335268 0.000107318 -0.000393101 -0.000291644 -0.000120409 -0.000123356 -0.000569329 -0.000271204 -2.36933e-5 7.34892e-5 -2.89015e-6; -0.000714703 -1.20491e-5 -4.93604e-5 -5.4529e-5 -0.000223464 9.03542e-6 -0.000140382 0.000248975 0.000110928 6.73943e-5 -2.89015e-6 8.56547e-5]
    #eyemat *= 0.5
    #covmat = [3.40274 1.64505 0.682729 -0.0227855 -0.0349619 -0.112405 0.170246 -0.0216856 0.0489216 -3.63159 0.0133547 0.00763219; 1.64505 1.38686 0.456436 0.0023172 0.00225685 0.0399517 -0.00897246 -0.0434067 0.0360939 -2.32204 0.00443568 0.00248812; 0.682729 0.456436 0.425997 0.00479458 -0.016015 0.0222947 0.0238418 -0.0882016 0.0414733 -0.498733 -0.00214618 -0.0033361; -0.0227855 0.0023172 0.00479458 0.064213 0.00900069 0.0171433 0.00116769 0.026957 0.00574038 -0.140947 -0.00378226 -0.00227964; -0.0349619 0.00225685 -0.016015 0.00900069 0.0451184 0.0121727 -1.20564e-5 0.0123196 0.00367476 -0.0702724 -0.00294692 -0.000957184; -0.112405 0.0399517 0.0222947 0.0171433 0.0121727 0.0610466 -0.0231408 0.00768127 0.0159818 0.111935 -0.00328597 0.000510824; 0.170246 -0.00897246 0.0238418 0.00116769 -1.20564e-5 -0.0231408 0.0774362 0.022108 0.0079497 -0.0209017 -0.000604395 -0.00163211; -0.0216856 -0.0434067 -0.0882016 0.026957 0.0123196 0.00768127 0.022108 0.245942 0.0301434 0.110812 -0.00576267 -0.00361273; 0.0489216 0.0360939 0.0414733 0.00574038 0.00367476 0.0159818 0.0079497 0.0301434 0.0690555 -0.199385 -0.00391895 -0.00186482; -3.63159 -2.32204 -0.498733 -0.140947 -0.0702724 0.111935 -0.0209017 0.110812 -0.199385 80.8044 -0.00435359 0.0433615; 0.0133547 0.00443568 -0.00214618 -0.00378226 -0.00294692 -0.00328597 -0.000604395 -0.00576267 -0.00391895 -0.00435359 0.00130361 0.000320698; 0.00763219 0.00248812 -0.0033361 -0.00227964 -0.000957184 0.000510824 -0.00163211 -0.00361273 -0.00186482 0.0433615 0.000320698 0.00287422]

    if len % 2 == 0
      covmat = e
      d2 = sum(mask)
      covmat = (2.38*2.38*covmat)/d2/2.0
      return current + maskarray(rand(MvNormal(covmat)),mask), covmat
    else
      covmat = e
      covmat = (2.38*2.38*covmat)/d/2.0
      return current + rand(MvNormal(covmat)), covmat
    end
  else
    e = eye(d)*0.0001
    eyemat = (e)/d2

    reg = 0.001*eye(d)/len

    startindex = max(1, div(length(tuningvectors),2))
    endindex = length(tuningvectors)
    mat = hcathelper(tuningvectors[startindex:endindex])
    mat = mat'

    covmat = (2.38*2.38*cov(mat))/d2
    #=
    if length(tuningvectors) % 10 == 0
      println(cov(mat))
    end=#
    #end

    r1 = maskarray(rand(MvNormal(eyemat)), mask)
    r2 = maskarray(rand(MvNormal((covmat.+reg)*covarScale)), mask)
    return C*(current+r1) + (1.0-C)*(current + r2), covmat
  end
end


function getbranchtunedparams(current::Array{Float64,1}, tuningvectors::Array{Array{Float64,1},1}, maskin::Array{Int,1}, covarScale::Float64=1.0, maxTuningIter::Int=1000, C::Float64=0.05)
  mask = copy(maskin)
  d = length(current)
  e = eye(d)*0.0001
  len = length(tuningvectors)
  if  len < 4*d
    eyemat = (e)/Float64(d)
    return current + maskarray(rand(MvNormal(eyemat)),mask)
  else
    temptuning = tuningvectors

    d2 = sum(mask)
    eyemat = (e)/d2

    startindex = max(3, div(length(temptuning),5))
    endindex = length(temptuning)
    mat = hcathelper(temptuning[startindex:endindex])
    mat = mat'

      reg = 0.0005*eye(d)/len
      covmat = (2.38*2.38*cov(mat))/d2
      if length(temptuning) % 10 == 0
        #println(cov(mat))
      end
    #end

    r1 = maskarray(rand(MvNormal(eyemat)), mask)
    r2 = maskarray(rand(MvNormal((covmat.+reg)*covarScale)), mask)

    return C*(current+r1) + (1.0-C)*(current + r2)
  end
end

mutable struct MCMCoptions
  mode::Int
  M::Float64
  samplebranchlengths::Bool
  integratesiterates::Bool
  integratestructure::Bool
  maxbasepairdistance::Int
  usecuda::Bool


  function MCMCoptions(mode::Int, M::Float64, sample::Bool, integratesiterates::Bool, integratestructure::Bool, maxbasepairdistance::Int, usecuda::Bool)
    new(mode, M, sample, integratesiterates, integratestructure, maxbasepairdistance, usecuda)
  end
end

function estimateinformationentropyhelper(rng::AbstractRNG, dataset::Dataset, params::ModelParameters, inside::Array{Float64,3}, Z::Float64, pairedlogprobs::Array{Float64,2}, unpairedlogprobs::Array{Float64,1}, calcpriorentropy::Bool)
  museparams = getmusespecificparamsarray(params)
  len = params.siteCats*params.siteCats*length(params.lambdaratesGC)
  cache = FastCache[FastCache(250000,16) for i=1:len]
  unpairedcache = FastCache[FastCache(500000,4) for i=1:params.siteCats]
  unpairedcolumncache = ColumnCache(500000)
  pairedcolumncache = ColumnCache(2500000)
  numsamples = 250
  grammar = KH99()
  v = Float64[]
  for s=1:numsamples
    paired = zeros(Int, dataset.numcols)
    priorll = samplestructure(rng, inside, pairedlogprobs, unpairedlogprobs, 1, dataset.numcols, paired, grammar, 1.0)
    ll = 0.0
    if !calcpriorentropy
      ll = calculatecachedlikelihood(dataset, params, museparams, params.states, paired, cache, unpairedcache, true,unpairedcolumncache,pairedcolumncache)
    end
    log2pd = (ll +  priorll - Z)/log(2.0)
    push!(v, log2pd)
  end
  return v
end


function estimateinformationentropy(rng::AbstractRNG, dataset::Dataset, params::ModelParameters, maxbasepairdistance::Int=500, usecuda::Bool=true, calcpriorentropy::Bool=false)
  unpairedlogprobs = zeros(Float64, dataset.numcols)
  pairedlogprobs = zeros(Float64, dataset.numcols, dataset.numcols)
  for i=1:dataset.numcols
    pairedlogprobs[i,i] = -Inf
  end
  ret = nothing

  if !calcpriorentropy
    unpairedlogprobs = computeunpairedlikelihoods(dataset, params)
    pairedlogprobs, ret = coevolutionall(dataset,params,true,true,true,maxbasepairdistance,usecuda)
    maskgapped!(pairedlogprobs,dataset.gapfrequency,0.5,-Inf)
  end
  inside = computeinsideKH99(unpairedlogprobs, pairedlogprobs, 1.0,false,usecuda)
  Z = inside[1,1,dataset.numcols]
  println("inside Z =", Z)
  e = 0.0

  #tic()
  total = -Inf

  v = Float64[]
  refs = []
  for p=1:16
    #ref = @spawn estimateinformationentropyhelper(MersenneTwister(rand(1:typemax(Int))), dataset, params, inside, Z, pairedlogprobs, unpairedlogprobs, calcpriorentropy)
    ref = estimateinformationentropyhelper(MersenneTwister(rand(1:typemax(Int))), dataset, params, inside, Z, pairedlogprobs, unpairedlogprobs, calcpriorentropy)
    push!(refs, ref)
  end

  for ref in refs
    append!(v, ref)
    #append!(v, fetch(ref))
  end

  #elapsed = toq()
  #println("Est time= ", elapsed)
  Hmax = 0.142 - 1.5*log2(dataset.numcols) + 1.388*dataset.numcols
  H = -mean(v)
  Hstdev = std(v)
  Hstderr = std(v)/sqrt(length(v))
  return H, Hstdev, Hstderr, Hmax, H/Hmax
end

function estimateinformationentropy2(rng::AbstractRNG, dataset::Dataset, params::ModelParameters, maxbasepairdistance::Int=500, usecuda::Bool=true)
  unpairedlogprobs = computeunpairedlikelihoods(dataset, params)
  pairedlogprobs, ret = coevolutionall(dataset,params,true,true,true,maxbasepairdistance,usecuda)
  maskgapped!(pairedlogprobs,dataset.gapfrequency,0.5,-Inf)
  inside = computeinsideKH99(unpairedlogprobs, pairedlogprobs, 1.0,false,usecuda)
  Z = inside[1,1,dataset.numcols]
  e = 0.0

  #tic()
  total = -Inf
  museparams = getmusespecificparamsarray(params)
  len = params.siteCats*params.siteCats*length(params.lambdarates)
  cache = FastCache[FastCache(250000,16) for i=1:len]


  unpairedcache = FastCache[FastCache(500000,4) for i=1:params.siteCats]
  unpairedcolumncache = ColumnCache(500000)
  pairedcolumncache = ColumnCache(2500000)
  v = Float64[]
  numsamples = 5000
  grammar = KH99()
  for s=1:numsamples
    paired = zeros(Int, dataset.numcols)
    priorll = samplestructure(rng, inside, pairedlogprobs, unpairedlogprobs, 1, dataset.numcols, paired, grammar, 1.0)
    ll = calculatecachedlikelihood(dataset, params, museparams, params.states, paired, cache, unpairedcache, true,unpairedcolumncache,pairedcolumncache)
    log2pd = (ll +  priorll - Z)/log(2.0)
    push!(v, log2pd)
  end

  #=
  final = Float64[]
  @sync for p=1:8
     w = @async begin
      unpairedcache = FastCache[FastCache(100000,4) for i=1:params.siteCats]
      unpairedcolumncache = ColumnCache(200000)
      pairedcolumncache = ColumnCache(5000000)
      v = Float64[]
      numsamples = 500
      for s=1:numsamples
        paired = zeros(Int, dataset.numcols)
        priorll = samplestructure(rng, inside, pairedlogprobs, unpairedlogprobs, 1, dataset.numcols, paired, KH99(), 1.0)
        ll = calculatecachedlikelihood(dataset, params, museparams, params.states, paired, cache, unpairedcache, true,unpairedcolumncache,pairedcolumncache)
        log2pd = (ll +  priorll - Z)/log(2.0)
        push!(v, log2pd)
      end
      return v
    end
    append!(final, w)
  end=#

  #elapsed = toq()
  #println("Est time= ", elapsed)
  Hmax = 0.142 - 1.5*log2(dataset.numcols) + 1.388*dataset.numcols
  H = -mean(v)
  Hstdev = std(v)
  Hstderr = std(v)/sqrt(length(v))
  return H, Hstdev, Hstderr, Hmax, H/Hmax
end

function computeinformationentropy(dataset::Dataset, params::ModelParameters, maxbasepairdistance::Int=500, usecuda::Bool=true, calcpriorentropy::Bool=false)
  unpairedlogprobs = zeros(Float64, dataset.numcols)
  pairedlogprobs = zeros(Float64, dataset.numcols, dataset.numcols)
  for i=1:dataset.numcols
    pairedlogprobs[i,i] = -Inf
  end
  ret = nothing

  unpairedlogprobs = computeunpairedlikelihoods(dataset, params)
  pairedlogprobs, ret = coevolutionall(dataset,params,true,true,true,maxbasepairdistance,usecuda)

  maskgapped!(pairedlogprobs,dataset.gapfrequency,0.5,-Inf)
  inside = computeinsideKH99(unpairedlogprobs, pairedlogprobs, 1.0,false,usecuda)
  outside = computeoutsideKH99(inside,unpairedlogprobs, pairedlogprobs, usecuda)
  unpairedposteriorprobs,pairedposteriorprobs = computebasepairprobs2(inside, outside, unpairedlogprobs, pairedlogprobs, KH99())
  Z = inside[1,1,dataset.numcols]
  println("inside Z =", Z)
  #println(exp(unpairedposteriorprobs-Z))
  #println(pairedlogprobs[1:10,1:10])
  #println(exp(pairedposteriorprobs[1:10,1:10]-Z))

  pairedentropy = 0.0
  unpairedentropy = 0.0
  for i=1:dataset.numcols
    if !isnan(unpairedposteriorprobs[i]) && unpairedlogprobs[i] != -Inf
      unpairedentropy += ((unpairedlogprobs[i])/log(2.0))*exp(unpairedposteriorprobs[i]-Z)
      #pairedentropy += ((unpairedlogprobs[i])/log(2.0))*exp(log(unpairedposteriorprobs[i])/log(2.0))
      #unpairedentropy += unpairedposteriorprobs[i]
      #unpairedentropy += log(unpairedposteriorprobs[i])/log(2.0)
    end
    for j=i+1:dataset.numcols
      if i != j
        if !isnan(pairedposteriorprobs[i,j]) && pairedlogprobs[i,j] != -Inf
          #pairedentropy += (log(pairedposteriorprobs[i,j])/log(2.0))*exp(log(pairedposteriorprobs[i,j])/log(2.0))
          pairedentropy += ((pairedlogprobs[i,j])/log(2.0))*exp(pairedposteriorprobs[i,j]-Z)
          #pairedentropy += log(pairedposteriorprobs[i,j])/log(2.0)
        end
      end
    end
  end

  println(Z)
  println(Z/log(2.0))
  #h = -(pairedentropy + unpairedentropy) + (Z/log(2.0))
  h = Z/log(2.0) - (pairedentropy + unpairedentropy)
  hmax = 0.142 - 1.5*log2(dataset.numcols) + 1.388*dataset.numcols
  return h, hmax
end

function runchain(runiters::Int, chain::Chain, dataset::Dataset, grammar::KH99, outputprefix::AbstractString, outputprefixchain::AbstractString, mcmcoptions::MCMCoptions, burnins::Array{Int,1}, burnindata::Array{Float64,2}=nothing, thermodynamicsamples::Dict{Int, Array{Array{Int,1}}}=Dict{Int, Array{Array{Int,1}}}(), tuningvectors::Array{Array{Float64,1},1}=Array{Float64,1}[], branchtuningvectors::Array{Array{Float64,1},1}=Array{Float64,1}[], covmatin=nothing; siteCats::Int=siteCats, lambdaCats::Int=lambdaCats, fixGU::Bool=false, fixGCAU::Bool=false)
  usemodelswitch = false
  mode = mcmcoptions.mode
  M = mcmcoptions.M
  samplebranchlengths = mcmcoptions.samplebranchlengths
  integratesiterates = mcmcoptions.integratesiterates
  integratestructure = mcmcoptions.integratestructure
  maxbasepairdistance = mcmcoptions.maxbasepairdistance
  usecuda = mcmcoptions.usecuda

  covmat = covmatin
  mcmclogfile = string(outputprefixchain, ".mcmc.log")
  mcmclog = open(mcmclogfile, chain.iter == 0 ? "w" : "a")
  structurelog = open(string(outputprefixchain, ".mcmc.structures"), chain.iter == 0 ? "w" : "a")
  consensuslog = open(string(outputprefixchain, ".mcmc.consensus"), chain.iter == 0 ? "w" : "a")
  siterateslog = open(string(outputprefixchain, ".mcmc.siterates"), chain.iter == 0 ? "w" : "a")
  if chain.iter == 0
    write(mcmclog, string("iter\tcurrentll\tproposedll\tpropratio\tlambdaGC\tlambdaAT\tlambdaGT\tlambdazeroweight\tlambdaGammaShapeGC\tlambdaGammaScaleGC\tlambdaGammaShapeAT\tlambdaGammaScaleAT\tlambdaGammaShapeGT\tlambdaGammaScaleGT\tsiteGammaShape\tsiteGammaScale\tqAC\tqAG\tqAT\tqCG\tqCT\tqGT\t"))
    if samplebranchlengths
      for b=2:length(chain.currentparams.branchlengths)
        write(mcmclog,string("t",b,"\t"))
      end
    end
    write(mcmclog,"piA\tpiC\tpiG\tpiU\t")
    for r=1:length(chain.currentparams.lambdaweightsGC)
      write(mcmclog, string("lambdaweights",r,"\t"))
    end
    for r=1:length(chain.currentparams.lambdaratesGC)
      write(mcmclog, string("lambdarateGC",r,"\t"))
    end
    for r=1:length(chain.currentparams.lambdaratesAT)
      write(mcmclog, string("lambdarateAT",r,"\t"))
    end
    for r=1:length(chain.currentparams.lambdaratesGT)
      write(mcmclog, string("lambdarateGT",r,"\t"))
    end
    for r=1:length(chain.currentparams.siteRates)
      write(mcmclog, string("siterate",r,"\t"))
    end
    for r=1:length(chain.currentparams.siteWeights)
      write(mcmclog, string("siteweights",r,"\t"))
    end
    write(mcmclog, string("iterspersecond\t"))
    write(mcmclog, string("hashindex\n"))
    write(siterateslog,"iter\t")
    for col=1:dataset.numcols
      write(siterateslog,string("site",col,"\t"))
    end
    write(siterateslog, "\n")
  end
  acceptancelog = open(string(outputprefixchain, ".mcmc.acceptance.log"), chain.iter == 0 ? "w" : "a")
  tempacceptancelog = open(string(outputprefixchain, ".mcmc.tempacceptance.log"), chain.iter == 0 ? "w" : "a")
  if usemodelswitch
      modelswitchlog = open(string(outputprefixchain, ".mcmc.modelswitch.log"), chain.iter == 0 ? "w" : "a")
      if chain.B == 1.0 && chain.iter == 0
        write(modelswitchlog, string("iter\tM\tU\tq0\tq1\tdelta\n"))
      end
  end

  paramfile = string(outputprefix,".params")
  if chain.iter == 0
    if !isfile(paramfile)
      jsondict = JSON.parse(open(joinpath(@__DIR__, "default.params")))
      out = open(paramfile,"w")
      ret = replace(JSON.json(jsondict),",\"" => ",\n\"")
      ret = replace(ret, "],[" => "],\n[")
      ret = replace(ret, "{" => "{\n")
      ret = replace(ret, "}" => "\n}")
      write(out,ret)
      close(out)
    end
  end

  jsondict = JSON.parse(open(paramfile))
  thishash = hash(jsondict)
  if chain.lasthash != thishash
    println("HASH CHANGED ", chain.lasthash, " -> ", thishash)
    chain.hashindex += 1
    println("HASH INDEX, ", chain.hashindex)
  end
  chain.lasthash = thishash

  moveweights = onedarray(jsondict["moveweights"])
  if samplebranchlengths
    push!(moveweights, 2.0)
  end
  if lambdaCats == 1
    moveweights[5] = 0.0
    moveweights[6] = 0.0
  end

  sigmas = ones(Float64,6+length(chain.currentparams.branchlengths)-1)*0.25
  for b=2:length(chain.currentparams.branchlengths)
    sigmas[6+b-1] = chain.currentparams.branchlengths[b]*0.2
  end

  lambdaSigma = eye(3)*1.2
  lambdaSigma[3,3] = 0.4
  lambdaSigma = twodarray(jsondict["lambdaSigma"])
  lambdaScalar = jsondict["lambdaScalar"]
  lambdaMv = MvNormal(lambdaSigma*lambdaScalar)

  lambdaGCSigma = jsondict["lambdaGCsigma"]
  lambdaATSigma = jsondict["lambdaATsigma"]
  lambdaGTSigma = jsondict["lambdaGTsigma"]

  QSigma = eye(6)*0.08
  if samplebranchlengths
      QSigma = eye(5)*0.1
  end
  QSigma = twodarray(jsondict["QSigma"])
  QMv = MvNormal(QSigma)

  gibbsInterval = jsondict["gibbsInterval"]
  qScalar = jsondict["qScalar"]

  q1sigma = jsondict["q1sigma"]
  q2sigma = jsondict["q2sigma"]
  q3sigma = jsondict["q3sigma"]
  q4sigma = jsondict["q4sigma"]
  q5sigma = jsondict["q5sigma"]
  q6sigma = jsondict["q6sigma"]

  lambdaGammaShape = jsondict["lambdaGammaShape"]
  lambdaZeroSigma = jsondict["lambdaZeroSigma"]

  maxTuningIter =  jsondict["maxTuningIter"]

  structures = Array{Int,1}[]
  for iter=chain.iter:chain.iter+runiters
    propratio = 0.0
    cont = true
    gibbs = false
    if length(burnins) > 0 && iter == burnins[1] && burnindata != nothing && chain.B == 1.0
      burnin = popfirst!(burnins)
      out = open(string(outputprefixchain,".",iter,".params"),"w")
      jsonout = deepcopy(jsondict)
      jsonout["lambdaSigma"] = cov(burnindata[div(burnin,3):burnin,1:3])
      jsonout["QSigma"] = cov(burnindata[div(burnin,3):burnin,4:9])
      Qvec = Float64[]
      for z=1:6
        push!(Qvec, std(burnindata[div(burnin,3):burnin,3+z]))
      end
      jsonout["Qvec"] = Qvec
      ret = replace(JSON.json(jsonout),",\"" => ",\n\"")
      ret = replace(ret, "],[" => "],\n[")
      ret = replace(ret, "{" => "{\n")
      ret = replace(ret, "}" => "\n}")
      write(out, ret)
      close(out)
      #=
      lambdaSigma = cov(burnindata[div(burnin,3):burnin,1:3])
      lambdaMv = MvNormal(lambdaScalar*lambdaSigma)
      QSigma = cov(burnindata[div(burnin,3):burnin,4:9])
      QMv = MvNormal(qScalar*QSigma)
      if samplebranchlengths
          QSigma = cov(burnindata[div(burnin,3):burnin,4:8])
          QMv = MvNormal(qScalar*QSigma)
      end
      for i=1:length(sigmas)
        if i <= 6
          sigmas[i] = std(burnindata[div(burnin,3):burnin,3+i])*4.0
        else
          #sigmas[i] = std(burnindata[div(burnin,3):burnin,3+i])*2.0
        end
      end
      jsondict["lambdaSigma"] = lambdaSigma
      jsondict["QSigma"] = QSigma
      ret = replace(JSON.json(jsondict),",\"", ",\n\"")
      ret = replace(ret, "],[", "],\n[")
      ret = replace(ret, "{", "{\n")
      ret = replace(ret, "}", "\n}")
      #ret = replace(JSON.json(jsondict),",[", ",\n[")
      out = open(string(outputprefixchain,".params"),"w")
      write(out, ret)
      close(out)=#
    end

    if mode == MODE_IO_SAMPLE && chain.B == 1.0 && iter % (gibbsInterval*2) == 0
      #=
      pairedlogprobs2, ret = coevolutionall(dataset,chain.currentparams,true,true,true,maxbasepairdistance)
      meanlambdas = zeros(Float64, dataset.numcols, dataset.numcols)
      posteriorone = zeros(Float64, dataset.numcols, dataset.numcols)
      totals = ones(Float64, dataset.numcols, dataset.numcols)*-Inf
      for (probs, lambdarate, prior) in zip(ret,chain.currentparams.lambdarates,chain.currentparams.lambdaweights)
        logprior = log(prior)
        for p=1:dataset.numcols
          for q=p+1:dataset.numcols
            totals[p,q] = CommonUtils.logsumexp(totals[p,q], logprior+probs[p,q])
          end
        end
      end
      for (probs, lambdarate, prior) in zip(ret,chain.currentparams.lambdarates,chain.currentparams.lambdaweights)
        logprior = log(prior)
        for p=1:dataset.numcols
          for q=p+1:dataset.numcols
            weight = exp(logprior+probs[p,q]-totals[p,q])
            meanlambdas[p,q] += weight*lambdarate
            meanlambdas[q,p] += weight*lambdarate
          end
        end
      end
      #h5write(string(outputprefixchain,".h5"), string("meanlambdas/",iter), meanlambdas)
      lambdaslog = open(string(outputprefixchain,".lambdas"), iter == 0 ? "w" : "a")
      write(lambdaslog, string(">",iter,"\n"))
      write(lambdaslog, string(meanlambdas,"\n"))
      close(lambdaslog)



      for p=1:dataset.numcols
        for q=p+1:dataset.numcols
          posteriorone[p,q] = 1.0 - exp(log(chain.currentparams.lambdaweights[1])+ret[1][p,q]-totals[p,q])
          posteriorone[q,p] = posteriorone[p,q]
        end
      end
      #h5write(string(outputprefixchain,".h5"), string("posteriorone/",iter), posteriorone)
      posteriorfile = string(outputprefixchain,".posterior")
      posteriorlog = open(posteriorfile, iter == 0 ? "w" : "a")
      write(posteriorlog, string(">",iter,"\n"))
      write(posteriorlog, string(posteriorone,"\n"))
      close(posteriorlog)

      basepairfile = string(outputprefixchain,".baseprobs")
      basepairlog = open(basepairfile, iter == 0 ? "w" : "a")
      outside = computeoutside(chain.inside,chain.unpairedlogprobs, chain.pairedlogprobs, grammar)
      unpairedposteriorprobs,pairedposteriorprobs = computebasepairprobs(chain.inside,outside,chain.unpairedlogprobs, chain.pairedlogprobs, grammar)
      #h5write(string(outputprefixchain,".h5"), string("unpairedposteriorprobs/",iter), unpairedposteriorprobs)
      #h5write(string(outputprefixchain,".h5"), string("pairedposteriorprobs/",iter), pairedposteriorprobs)
      write(basepairlog, string(">",iter,"\n"))
      write(basepairlog, string(pairedposteriorprobs,"\n"))
      close(basepairlog)
      =#
    end

    if mode == MODE_FIXED_STRUCTURE || iter % 50 == 0
      #=
      museparams = getmusespecificparamsarray(chain.currentparams)
      sitelikelihoods,sampledstates,museconditionals = testpaired(dataset, chain.currentparams, museparams, chain.currentparams.states, chain.paired, chain.B, false, chain.rng)
      GUarray = zeros(Float64,dataset.numcols)
      for i=1:dataset.numcols
        if chain.paired[i] > i
          lambdaGUll = -Inf
          total = -Inf
          for musespecificparams in museparams
            if musespecificparams.lambdarate != 0.0
              if musespecificparams.lambdaGT == 1.0
                lambdaGUll = CommonUtils.logsumexp(lambdaGUll, museconditionals[i,musespecificparams.lambdacat])
              end
              total = CommonUtils.logsumexp(total, museconditionals[i,musespecificparams.lambdacat])
            end
          end
          GUarray[i] = exp(lambdaGUll-total)
        end
      end
      posteriorfile = string(outputprefixchain,".gu")
      posteriorlog = open(posteriorfile, iter == 0 ? "w" : "a")
      if iter == 0
        write(posteriorlog, "iter\t")
        for i=1:dataset.numcols
          write(posteriorlog, string("site",i,"\t"))
        end
        write(posteriorlog, "\n")
      end
      write(posteriorlog, string(iter, "\t", join(GUarray, "\t"),"\n"))
      close(posteriorlog)
      =#
    end

    if (mode == MODE_FIXED_STRUCTURE || mode == MODE_VIENNA_SAMPLE) && iter % 50 == 0
      museparams = getmusespecificparamsarray(chain.currentparams)
      unpairedloglikelihoods,pairedloglikelihoods,chain.currentparams.states,museconditionals = computeuncachedlikelihood(dataset, chain.currentparams, museparams, chain.currentparams.states, chain.paired, 1.0, false, chain.rng)

      posteriorone = zeros(Float64, size(museconditionals,1))
      posteriorlambda = zeros(Float64, size(museconditionals,1))
      for i=1:size(museconditionals,1)
        if chain.paired[i] > i
          posteriorone[i] = 1.0 - exp(museconditionals[i,1]-pairedloglikelihoods[i])
          for j=1:chain.currentparams.lambdaCats
            #posteriorlambda[i] += exp(museconditionals[i,j]-pairedloglikelihoods[i])*chain.currentparams.lambdarates[j]
          end
          #posteriorlambda[i]  = sum(exp(likelihoods[i,:]-total).*chain.currentparams.lambdarates)
        else
          posteriorone[i] = NaN
          posteriorlambda[i] = NaN
        end
      end
      posteriorfile = string(outputprefixchain,".posterior")
      posteriorlambdafile = string(outputprefixchain,".posteriorlambda")
      posteriorlog = open(posteriorfile, iter == 0 ? "w" : "a")
      posteriorlambdalog = open(posteriorlambdafile, iter == 0 ? "w" : "a")
      if iter == 0
        write(posteriorlog, "iter\t")
        write(posteriorlambdalog, "iter\t")
        for i=1:size(museconditionals,1)
          write(posteriorlog, string("site",i,"\t"))
          write(posteriorlambdalog, string("site",i,"\t"))
        end
        write(posteriorlog, "\n")
        write(posteriorlambdalog, "\n")
      end
      write(posteriorlog, string(iter, "\t", join(posteriorone, "\t"),"\n"))
      write(posteriorlambdalog, string(iter, "\t", join(posteriorlambda, "\t"),"\n"))
      close(posteriorlog)
      close(posteriorlambdalog)

    end

    if !integratesiterates && integratestructure && rand(chain.rng, 1:gibbsInterval) == gibbsInterval
      gibbs = true
      unpairedlogprobs = computeunpairedlikelihoods(dataset, chain.currentparams)
      pairedlogprobs, ret = coevolutionall(dataset,chain.currentparams,true,true,true,maxbasepairdistance,usecuda)
      maskgapped!(pairedlogprobs,dataset.gapfrequency,0.5,-Inf)
      #inside = computeinsideparallel(unpairedlogprobs, pairedlogprobs, KH99(), chain.B)
      inside = computeinsideKH99(unpairedlogprobs, pairedlogprobs, 1.0, false,usecuda)
      Z = inside[1,1,dataset.numcols]
      savemaximum2(Z, getparamsvector2(chain.currentparams), string(outputprefix,".maxmcmc"))
      println("inside Z =", Z)
      outside = computeoutside(inside,unpairedlogprobs, pairedlogprobs, grammar)
      unpairedposteriorprobs,pairedposteriorprobs = computebasepairprobs(inside, outside, unpairedlogprobs, pairedlogprobs, grammar)
      unpairedpostlogprobs = log(unpairedposteriorprobs)
      pairedpostlogprobs = log(pairedposteriorprobs)

      museparams = getmusespecificparamsarray(chain.currentparams)
      params = chain.currentparams
      siteconditionals = ones(Float64,dataset.numcols, params.siteCats)*-Inf
      for x=1:dataset.numcols
        yend = min(x+maxbasepairdistance,dataset.numcols)
        for y=1:dataset.numcols
          i = 1
          for j=1:params.siteCats
            logj = log(params.siteWeights[j])
            for k=1:params.siteCats
              logk = log(params.siteWeights[k])
              for musespecificparam in museparams
                loglambdaweight = musespecificparam.logprob + logj + logk
                ll = loglambdaweight + ret[i][x,y] + pairedpostlogprobs[x,y]
                siteconditionals[x, ((j-1) % params.siteCats) + 1] = CommonUtils.logsumexp(siteconditionals[x, ((j-1) % params.siteCats) + 1], ll)
                siteconditionals[y, ((k-1) % params.siteCats) + 1] = CommonUtils.logsumexp(siteconditionals[y, ((k-1) % params.siteCats) + 1], ll)
                i += 1
              end
            end
          end
        end
      end

      for siteCat=1:params.siteCats
        temp = computeunpairedlikelihoodscat(dataset, params, Int[siteCat for i=1:dataset.numcols], Int[i for i=1:dataset.numcols], siteCat)
        logj = log(params.siteWeights[siteCat])
        for x=1:dataset.numcols
          siteconditionals[x, siteCat] = CommonUtils.logsumexp(siteconditionals[x, siteCat], logj+unpairedpostlogprobs[x]+temp[x])
        end
      end

      sampledstates = zeros(Int,dataset.numcols)
      for x=1:dataset.numcols
        total = -Inf
        for siteCat=1:params.siteCats
          total = CommonUtils.logsumexp(total, siteconditionals[x,siteCat])
        end
        siteconditionals[x,:] = exp(siteconditionals[x,:]-total)
        sampledstates[x] = CommonUtils.sample(chain.rng, siteconditionals[x,:])
      end

      tempparams = deepcopy(chain.currentparams)
      tempparams.states = copy(sampledstates)


      unpairedlogprobs = computeunpairedlikelihoods(dataset, chain.currentparams, chain.currentparams.states)
      pairedlogprobs, ret = coevolutionall(dataset,chain.currentparams,true,false,false,maxbasepairdistance,usecuda)
      propold = computeinsideKH99(unpairedlogprobs, pairedlogprobs, 1.0,false,true)[1,1,dataset.numcols]

      unpairedlogprobs = computeunpairedlikelihoods(dataset, tempparams, tempparams.states)
      pairedlogprobs, ret = coevolutionall(dataset,tempparams,true,false,false,maxbasepairdistance,usecuda)
      propnew = computeinsideKH99(unpairedlogprobs, pairedlogprobs, 1.0,false,true)[1,1,dataset.numcols]


      oldZ = computetotallikelihood(chain.rng, dataset, chain.currentparams, chain.paired, samplebranchlengths,false,true,M,true,maxbasepairdistance, usecuda)*chain.B
      newZ = computetotallikelihood(chain.rng, dataset, tempparams, chain.paired, samplebranchlengths,false,true,M,true,maxbasepairdistance, usecuda)*chain.B

      propratio = propold - propnew
      delta = newZ - oldZ
      println(tempparams.states)
      println("ratio=", exp(delta + propratio), "\t", delta, "\t", propratio, "\t", chain.B)
      if exp(delta + propratio) > rand(chain.rng)
        accept = true
        chain.currentparams.states = copy(tempparams.states)
        chain.currentll = newZ
        logAccept!(chain.logger,"sitegibbs")
      else
        logReject!(chain.logger,"sitegibbs")
      end
    elseif integratesiterates && !integratestructure &&  mode == MODE_IO_SAMPLE && rand(chain.rng, 1:gibbsInterval) == gibbsInterval
      gibbs = true
      if M == 1.0
        tempparams = deepcopy(chain.currentparams)
        unpairedlogprobs = computeunpairedlikelihoods(dataset, tempparams)
        #tic()
        pairedlogprobs, ret = coevolutionall(dataset,tempparams,true,true,true,maxbasepairdistance,usecuda)
        #elapsed = toq()
        if usecuda
          #logtime("allpairs_cuda", elapsed, dataset.numcols)
        else
          #logtime("allpairs_cpu", elapsed, dataset.numcols)
        end

        maskgapped!(pairedlogprobs,dataset.gapfrequency,0.5,-Inf)
        #tic()
        inside = computeinsideKH99(unpairedlogprobs, pairedlogprobs, 1.0,false,usecuda)
        Z = inside[1,1,dataset.numcols]
        savemaximum2(Z, getparamsvector2(tempparams), string(outputprefix,".maxmcmc"))
        println("inside Z =", Z)
        #elapsed = toq()
        if usecuda
          #logtime("inside_cuda", elapsed, dataset.numcols)
        else
          #logtime("inside_cpu", elapsed, dataset.numcols)
        end
        chain.pairedlogprior = samplestructure(chain.rng, inside, pairedlogprobs, unpairedlogprobs, 1, dataset.numcols, chain.paired, grammar, chain.B)
      elseif M == 0.0
        tempparams = deepcopy(chain.currentparams)
        tempparams.lambdaGT = 1.0
        unpairedlogprobs = computeunpairedlikelihoods(dataset, tempparams)
        pairedlogprobs, ret = coevolutionall(dataset,tempparams,true,true,true,maxbasepairdistance,usecuda)
        maskgapped!(pairedlogprobs,dataset.gapfrequency,0.5,-Inf)
        inside = computeinsideKH99(unpairedlogprobs, pairedlogprobs, 1.0,false,usecuda)
        chain.pairedlogprior = samplestructure(chain.rng, inside, pairedlogprobs, unpairedlogprobs, 1, dataset.numcols, chain.paired, grammar, chain.B)
      else
        tempparams = deepcopy(chain.currentparams)
        unpairedlogprobs1 = computeunpairedlikelihoods(dataset, tempparams)
        pairedlogprobs1, ret = coevolutionall(dataset,tempparams,true,true,true,maxbasepairdistance,usecuda)
        maskgapped!(pairedlogprobs1,dataset.gapfrequency,0.5,-Inf)
        inside1 = computeinsideKH99(unpairedlogprobs1, pairedlogprobs1, 1.0,false,usecuda)

        tempparams.lambdaGT = 1.0
        unpairedlogprobs0 = computeunpairedlikelihoods(dataset, tempparams)
        pairedlogprobs0, ret = coevolutionall(dataset,tempparams,true,true,true,maxbasepairdistance,usecuda)
        maskgapped!(pairedlogprobs0,dataset.gapfrequency,0.5,-Inf)
        inside0 = computeinsideKH99(unpairedlogprobs0, pairedlogprobs0, 1.0,false,usecuda)

        v = [log(1.0 - M)+inside0[1,1,dataset.numcols], log(M)+inside1[1,1,dataset.numcols]]*chain.B
        total = CommonUtils.logsumexp(v[1],v[2])
        v = exp(v-total)
        r = CommonUtils.sample(chain.rng, v)
        println(v,"\t",r)
        if r == 1
            chain.pairedlogprior = samplestructure(chain.rng, inside0, pairedlogprobs0, unpairedlogprobs0, 1, dataset.numcols, chain.paired, grammar, chain.B)
        elseif r == 2
            chain.pairedlogprior = samplestructure(chain.rng, inside1, pairedlogprobs1, unpairedlogprobs1, 1, dataset.numcols, chain.paired, grammar, chain.B)
        end
      end
      chain.currentll = computetotallikelihood(chain.rng, dataset, chain.currentparams, chain.paired, samplebranchlengths,integratesiterates,integratestructure,M,true,maxbasepairdistance, usecuda)
      logAccept!(chain.logger,"structuregibbs")
      #=

      tempparams = deepcopy(chain.currentparams)
      if rand(chain.rng) > M
        tempparams.lambdaGT = 1.0
      end

      gibbs = true

      unpairedlogprobs = computeunpairedlikelihoods(dataset, tempparams)
      pairedlogprobs, ret = coevolutionall(dataset,tempparams,true,true,true,500)
      maskgapped!(pairedlogprobs,dataset.gapfrequency,0.5,-Inf)
      inside = computeinsidecuda(unpairedlogprobs, pairedlogprobs, 1.0)
      #outside = computeoutside(inside,unpairedlogprobs, pairedlogprobs, grammar)

      for z=1:5
        oldpairedprior = chain.pairedlogprior
        oldpaired = copy(chain.paired)
        propold = computetotallikelihood(chain.rng, dataset, tempparams, chain.paired, samplebranchlengths,integratesiterates,false,1.0)+oldpairedprior
        oldll = computetotallikelihood(chain.rng, dataset, chain.currentparams, chain.paired, samplebranchlengths,integratesiterates,false,M)+oldpairedprior

        newpairedprior = samplestructure(chain.rng, inside, pairedlogprobs, unpairedlogprobs, 1, dataset.numcols, chain.paired, grammar, chain.B)
        propnew = computetotallikelihood(chain.rng, dataset, tempparams, chain.paired, samplebranchlengths,integratesiterates,false,1.0)+newpairedprior
        newll = computetotallikelihood(chain.rng, dataset, chain.currentparams, chain.paired, samplebranchlengths,integratesiterates,false,M)+newpairedprior

        propratio = (propold - propnew)*chain.B
        delta = (newll - oldll)*chain.B
        println("ratio=", exp(delta + propratio), "\t", delta, "\t", propratio, "\t", chain.B, "\t", oldpairedprior, "\t", newpairedprior)
        if exp(delta + propratio) > rand(chain.rng)
          accept = true
          chain.currentll = propnew
          chain.pairedlogprior = newpairedprior
          logAccept!(chain.logger,"structuregibbs")
        else
          chain.paired = copy(oldpaired)
          chain.pairedlogprior = oldpairedprior
          logReject!(chain.logger,"structuregibbs")
        end
      end

      =#

      #=
      #if rand() < 1.0
        #useprev = rand(chain.rng) < 0.80
        useprev = false
        if !useprev
          tempparams = deepcopy(chain.currentparams)
          tempparams.lambdaGT = M*chain.currentparams.lambdaGT + (1.0-M)*1.0

          samplestates(chain.rng, dataset, tempparams, chain.paired, chain.B)
          chain.unpairedlogprobs = computeunpairedlikelihoods(dataset, tempparams, tempparams.states)
          #tic()
          chain.pairedlogprobs, ret = coevolutionall(dataset,tempparams,true,false,false,maxbasepairdistance)
          maskgapped!(chain.pairedlogprobs,dataset.gapfrequency,0.5,-Inf)
          elapsed = #toc()
          colspersecond = (dataset.numcols*(dataset.numcols-1.0)/2.0) / elapsed
          println("coevolution=",elapsed)
          println(colspersecond)
          #tic()
          chain.inside = computeinsidecuda(chain.unpairedlogprobs, chain.pairedlogprobs, chain.B)
          chain.maxparams = deepcopy(tempparams)
          #chain.paired = zeros(Int,dataset.numcols)
          elapsed = #toc()
          println("inside=",elapsed)
        end


        if integratesiterates

           #cache = FastCache[FastCache(125000,16) for i=1:len]
           #unpairedcache = FastCache[FastCache(100000,4) for i=1:chain.currentparams.siteCats]
           #propratio  = calculatecachedlikelihood(dataset, chain.currentparams, museparams,chain.currentparams.states, chain.paired, cache, unpairedcache, false)
           accept = false
           if !useprev
             for z=1:5
               oldpaired = copy(chain.paired)
               samplestructure(chain.rng, chain.inside, chain.pairedlogprobs, chain.unpairedlogprobs, 1, dataset.numcols, chain.paired, grammar, chain.B)
               propratio = computetotallikelihood(chain.rng, dataset, chain.maxparams, oldpaired, samplebranchlengths,false,integratestructure,1.0,1.0)*chain.B - computetotallikelihood(chain.rng, dataset, chain.maxparams, chain.paired, samplebranchlengths,false,integratestructure,1.0,1.0)*chain.B
               delta = computetotallikelihood(chain.rng, dataset, chain.currentparams, chain.paired, samplebranchlengths,true,integratestructure,1.0,M)*chain.B - computetotallikelihood(chain.rng, dataset, chain.currentparams, oldpaired, samplebranchlengths,true,integratestructure,1.0,M)*chain.B
               println("Structure proposal=", exp(delta + propratio))
               if exp(delta + propratio) > rand(chain.rng)
                 accept = true
               else
                 chain.paired = copy(oldpaired)
               end
             end
           else
             #tic()
             len = chain.currentparams.siteCats*chain.currentparams.siteCats*length(chain.currentparams.lambdarates)
             cache1 = FastCache[FastCache(250000,16) for i=1:len]
             unpairedcache1 = FastCache[FastCache(100000,4) for i=1:chain.currentparams.siteCats]
             unpairedcolumncache1 = ColumnCache(100000)
             pairedcolumncache1 = ColumnCache(1000000)
             museparams1 = getmusespecificparamsarray(chain.maxparams)

             cache2 = FastCache[FastCache(250000,16) for i=1:len]
             unpairedcache2 = FastCache[FastCache(100000,4) for i=1:chain.currentparams.siteCats]
             unpairedcolumncache2 = ColumnCache(100000)
             pairedcolumncache2 = ColumnCache(1000000)
             museparams2 = getmusespecificparamsarray(chain.currentparams)

             for z=1:5
               oldpaired = copy(chain.paired)
               samplestructure(chain.rng, chain.inside, chain.pairedlogprobs, chain.unpairedlogprobs, 1, dataset.numcols, chain.paired, grammar, chain.B)
               #propratio = computetotallikelihood(chain.rng, dataset, chain.maxparams, oldpaired, samplebranchlengths,true,1.0) - computetotallikelihood(chain.rng, dataset, chain.maxparams, chain.paired, samplebranchlengths,true,1.0)
               propratio  = calculatecachedlikelihood(dataset, chain.maxparams, museparams1, chain.maxparams.states, oldpaired, cache1, unpairedcache1, false, unpairedcolumncache1, pairedcolumncache1)*chain.B
               propratio -= calculatecachedlikelihood(dataset, chain.maxparams, museparams1, chain.maxparams.states, chain.paired, cache1, unpairedcache1, false, unpairedcolumncache1, pairedcolumncache1)*chain.B
               delta = calculatecachedlikelihood(dataset, chain.currentparams, museparams2, chain.currentparams.states, chain.paired, cache2, unpairedcache2, true, unpairedcolumncache2, pairedcolumncache2)*chain.B
               delta -= calculatecachedlikelihood(dataset, chain.currentparams, museparams2, chain.currentparams.states, oldpaired, cache2, unpairedcache2, true, unpairedcolumncache2, pairedcolumncache2)*chain.B

               #delta = computetotallikelihood(chain.rng, dataset, chain.currentparams, chain.paired, samplebranchlengths,true,M)*chain.B - computetotallikelihood(chain.rng, dataset, chain.currentparams, oldpaired, samplebranchlengths,true,M)*chain.B
               #println("Fast structure proposal=", exp(delta + propratio))
               if exp(delta + propratio) > rand(chain.rng)
                 accept = true
                 #logAccept!(chain.logger,"structuregibbs")
               else
                 chain.paired = copy(oldpaired)
                 #logReject!(chain.logger,"structuregibbs")
               end
             end
             #toc()
             println("Fast structure gibbs=", accept)

           end

           if accept
             if useprev
               logAccept!(chain.logger,"faststructuregibbs")
             else
               logAccept!(chain.logger,"structuregibbs")
             end
           else
             if useprev
               logReject!(chain.logger,"faststructuregibbs")
             else
               logReject!(chain.logger,"structuregibbs")
             end
           end

           chain.currentll = computetotallikelihood(chain.rng, dataset, chain.currentparams, chain.paired, samplebranchlengths,true,integratestructure,1.0,M)=
        else
          samplestructure(chain.rng, chain.inside, chain.pairedlogprobs, chain.unpairedlogprobs, 1, dataset.numcols, chain.paired, grammar, chain.B)
          chain.currentll = computetotallikelihood(chain.rng, dataset, chain.currentparams, chain.paired, samplebranchlengths,integratesiterates,integratestructure,1.0,M)
        end
        =#
    elseif mode == MODE_VIENNA_SAMPLE && rand(chain.rng, 1:gibbsInterval) == gibbsInterval
      gibbs = true

      println("Thermodynamic sampling")

      #tic()
      museparams = getmusespecificparamsarray(chain.currentparams)
      len = chain.currentparams.siteCats*chain.currentparams.siteCats*length(chain.currentparams.lambdarates)
      cache = FastCache[FastCache(125000,16) for i=1:len]
      unpairedcache = FastCache[FastCache(100000,4) for i=1:chain.currentparams.siteCats]
      unpairedcolumncache = ColumnCache(100000)
      pairedcolumncache = ColumnCache(1000000)
      llcurrent = calculatecachedlikelihood(dataset, chain.currentparams, museparams,chain.currentparams.states, chain.paired, cache, unpairedcache, integratesiterates,unpairedcolumncache,pairedcolumncache)
      accept = 0.0
      ziter = 100
      for z=1:ziter
        temp = copy(chain.paired)
        chain.paired = samplethermodynamic(thermodynamicsamples, chain.rng,dataset.sequences)
        #chain.paired = maskgapped(chain.paired, dataset.gapfrequency, 0.5, 0)
        llsample = calculatecachedlikelihood(dataset, chain.currentparams,museparams,chain.currentparams.states, chain.paired, cache, unpairedcache, integratesiterates,unpairedcolumncache,pairedcolumncache)
        if exp((llsample-llcurrent)*chain.B) > rand(chain.rng)
          llcurrent = llsample
          accept += 1.0
          #println("accept")
        else
          chain.paired = temp
          #println("reject")
        end
      end
      println("gibbsratio=", (accept/Float64(ziter)))
      #toc()
      chain.currentll =  computetotallikelihood(chain.rng, dataset, chain.currentparams, chain.paired, samplebranchlengths,integratesiterates,integratestructure,M,true,maxbasepairdistance, usecuda)
    else
      move = CommonUtils.sample(chain.rng,moveweights)
      movename = string(move)
      if move == 1
        if rand(chain.rng) < 0.75
          randsigma = rand(lambdaMv)
          chain.proposedparams.lambdaGC = chain.currentparams.lambdaGC + randsigma[1]
          chain.proposedparams.lambdaAT = chain.currentparams.lambdaAT + randsigma[2]
          chain.proposedparams.lambdaGT = chain.currentparams.lambdaGT + randsigma[3]
          movename = "lambdacorr"
        else
          sel = rand(chain.rng, 1:3)
          if sel == 1
            d1 = Truncated(Normal(chain.currentparams.lambdaGC, lambdaGCSigma), 0.0, Inf)
            chain.proposedparams.lambdaGC  = rand(d1)
            d2 = Truncated(Normal(chain.proposedparams.lambdaGC, lambdaGCSigma), 0.0, Inf)
            propratio += logpdf(d2, chain.currentparams.lambdaGC) - logpdf(d1, chain.proposedparams.lambdaGC)
            movename = "lambdaGC"
          elseif sel == 2
            d1 = Truncated(Normal(chain.currentparams.lambdaAT, lambdaATSigma), 0.0, Inf)
            chain.proposedparams.lambdaAT  = rand(d1)
            d2 = Truncated(Normal(chain.proposedparams.lambdaAT, lambdaATSigma), 0.0, Inf)
            propratio += logpdf(d2, chain.currentparams.lambdaAT) - logpdf(d1, chain.proposedparams.lambdaAT)
            movename = "lambdaAT"
          elseif sel == 3
            d1 = Truncated(Normal(chain.currentparams.lambdaGT, lambdaGTSigma), 0.0, Inf)
            chain.proposedparams.lambdaGT  = rand(d1)
            d2 = Truncated(Normal(chain.proposedparams.lambdaGT, lambdaGTSigma), 0.0, Inf)
            propratio += logpdf(d2, chain.currentparams.lambdaGT) - logpdf(d1, chain.proposedparams.lambdaGT)
            movename = "lambdaGT"
          end
        end
      elseif move == 2
        r = rand(chain.rng, 1:6)
        qmultiplier = 2.0
        qall = false
        if rand() < 0.5
          qmultiplier *= jsondict["qmultiplier"]
          qall = true
        end
        Qvec = jsondict["Qvec"]
        if qall || r == 1
          d1 = Truncated(Normal(chain.currentparams.q1, Qvec[1]*qmultiplier), 0.0, Inf)
          chain.proposedparams.q1  = rand(d1)
          d2 = Truncated(Normal(chain.proposedparams.q1, Qvec[1]*qmultiplier), 0.0, Inf)
          propratio += logpdf(d2, chain.currentparams.q1) - logpdf(d1, chain.proposedparams.q1)
        end
        if qall || r == 2
          d1 = Truncated(Normal(chain.currentparams.q2, Qvec[2]*qmultiplier), 0.0, Inf)
          chain.proposedparams.q2  = rand(d1)
          d2 = Truncated(Normal(chain.proposedparams.q2, Qvec[2]*qmultiplier), 0.0, Inf)
          propratio += logpdf(d2, chain.currentparams.q2) - logpdf(d1, chain.proposedparams.q2)
        end
        if qall || r == 3
          d1 = Truncated(Normal(chain.currentparams.q3, Qvec[3]*qmultiplier), 0.0, Inf)
          chain.proposedparams.q3  = rand(d1)
          d2 = Truncated(Normal(chain.proposedparams.q3, Qvec[3]*qmultiplier), 0.0, Inf)
          propratio += logpdf(d2, chain.currentparams.q3) - logpdf(d1, chain.proposedparams.q3)
        end
        if qall || r == 4
          d1 = Truncated(Normal(chain.currentparams.q4, Qvec[4]*qmultiplier), 0.0, Inf)
          chain.proposedparams.q4  = rand(d1)
          d2 = Truncated(Normal(chain.proposedparams.q4, Qvec[4]*qmultiplier), 0.0, Inf)
          propratio += logpdf(d2, chain.currentparams.q4) - logpdf(d1, chain.proposedparams.q4)
        end
        if qall || r == 5
          d1 = Truncated(Normal(chain.currentparams.q5, Qvec[5]*qmultiplier), 0.0, Inf)
          chain.proposedparams.q5  = rand(d1)
          d2 = Truncated(Normal(chain.proposedparams.q5, Qvec[5]*qmultiplier), 0.0, Inf)
          propratio += logpdf(d2, chain.currentparams.q5) - logpdf(d1, chain.proposedparams.q5)
        end
        if qall || r == 6
          d1 = Truncated(Normal(chain.currentparams.q6, Qvec[6]*qmultiplier), 0.0, Inf)
          chain.proposedparams.q6  = rand(d1)
          d2 = Truncated(Normal(chain.proposedparams.q6, Qvec[6]*qmultiplier), 0.0, Inf)
          propratio += logpdf(d2, chain.currentparams.q6) - logpdf(d1, chain.proposedparams.q6)
        end
        movename = string("gtr_", r)
        if qall
          movename = string("gtrall")
        end
      elseif move == 3
        len = length(chain.proposedparams.freqs)
        r1 = rand(chain.rng, 1:len)
        r2 = rand(chain.rng, 1:len)
        while r1 == r2
          r2 = rand(chain.rng, 1:len)
        end
        rdelta = rand(chain.rng)*jsondict["pisigma"]
        chain.proposedparams.freqs[r1] += rdelta
        chain.proposedparams.freqs[r2] -= rdelta
        if chain.proposedparams.freqs[r1] >= 1.0 || chain.proposedparams.freqs[r2] < 0.0
          cont = false
        end
        movename = string("pi",r1)
      elseif move == 4
        if !integratesiterates && rand() < 0.5
          #tic()
          samplestates(chain.rng, dataset, chain.proposedparams, chain.paired, chain.B)
          #e3 = #toc()
          gibbs = true
          movename = "states"
        else
          #scale = jsondict["siteGammaShape"]
          #=scale = jsondict["siteGammaScale"]

          d1 = Gamma(chain.currentparams.siteGammaShape / scale, scale)
          chain.proposedparams.siteGammaShape  = rand(d1)
          d2 = Gamma(chain.proposedparams.siteGammaShape / scale, scale)
          propratio += logpdf(d2, chain.currentparams.siteGammaShape) - logpdf(d1, chain.proposedparams.siteGammaShape)
          chain.proposedparams.siteGammaScale = 1.0 / chain.proposedparams.siteGammaShape
          chain.proposedparams.siteRates = discretizegamma(chain.proposedparams.siteGammaShape, chain.proposedparams.siteGammaScale, chain.proposedparams.siteCats)
          movename = "siteGammaShape"
          =#

          sigma = jsondict["siteGammaShape"]
          d1 = Truncated(Normal(chain.currentparams.siteGammaShape, sigma), 0.0, Inf)
          chain.proposedparams.siteGammaShape  = rand(d1)
          d2 = Truncated(Normal(chain.proposedparams.siteGammaShape, sigma), 0.0, Inf)
          propratio += logpdf(d2, chain.currentparams.siteGammaShape) - logpdf(d1, chain.proposedparams.siteGammaShape)
          chain.proposedparams.siteGammaScale = 1.0 / chain.proposedparams.siteGammaShape
          chain.proposedparams.siteRates = discretizegamma(chain.proposedparams.siteGammaShape, chain.proposedparams.siteGammaScale, chain.proposedparams.siteCats)
          movename = "siteGammaShape"
        end
      elseif move == 5
        #=
        scale = jsondict["lambdaGammaScale"]
        d1 = Gamma(chain.currentparams.lambdaGammaShape / scale, scale)
        chain.proposedparams.lambdaGammaShape  = rand(d1)
        d2 = Gamma(chain.proposedparams.lambdaGammaShape / scale, scale)
        propratio += logpdf(d2, chain.currentparams.lambdaGammaShape) - logpdf(d1, chain.proposedparams.lambdaGammaShape)
        chain.proposedparams.lambdaGammaScale = 1.0 / chain.proposedparams.lambdaGammaShape
        chain.proposedparams.lambdaweights, chain.proposedparams.lambdarates = discretizegamma2(chain.proposedparams.lambdaGammaShape, chain.proposedparams.lambdaGammaScale, chain.proposedparams.lambdaCats)
        movename = "lambdaGammaShape"
        =#

        sigma = jsondict["lambdaGammaShape"]
        d1 = Truncated(Normal(chain.currentparams.lambdaGammaShape, sigma), 0.0, Inf)
        chain.proposedparams.lambdaGammaShape  = rand(d1)
        d2 = Truncated(Normal(chain.proposedparams.lambdaGammaShape, sigma), 0.0, Inf)
        propratio += logpdf(d2, chain.currentparams.lambdaGammaShape) - logpdf(d1, chain.proposedparams.lambdaGammaShape)
        movename = "lambdaGammaShape"
        #chain.proposedparams.lambdaGammaScale = 1.0 / chain.proposedparams.lambdaGammaShape
        #chain.proposedparams.lambdaweights, chain.proposedparams.lambdarates = discretizegamma2(chain.proposedparams.lambdazeroweight, chain.proposedparams.lambdaGammaShape, chain.proposedparams.lambdaGammaScale, chain.proposedparams.lambdaCats)
        chain.proposedparams.lambdaweightsGC, chain.proposedparams.lambdaratesGC, chain.proposedparams.lambdaweightsAT, chain.proposedparams.lambdaratesAT, chain.proposedparams.lambdaweightsGT, chain.proposedparams.lambdaratesGT = discretizegamma3(chain.proposedparams.lambdazeroweight, chain.proposedparams.lambdaGammaShapeGC, chain.proposedparams.lambdaGammaShapeAT, chain.proposedparams.lambdaGammaShapeGT, chain.proposedparams.lambdaGammaScaleGC, chain.proposedparams.lambdaCats)


        #=
        len = length(chain.proposedparams.lambdarates)
        r = rand(chain.rng, 2:len)
        sigma = jsondict["lambdaratesigma"]
        d1 = Truncated(Normal(chain.currentparams.lambdarates[r], sigma), 0.0, Inf)
        chain.proposedparams.lambdarates[r]  = rand(d1)
        d2 = Truncated(Normal(chain.proposedparams.lambdarates[r], sigma), 0.0, Inf)
        propratio += logpdf(d2, chain.currentparams.lambdarates[r]) - logpdf(d1, chain.proposedparams.lambdarates[r])
        movename = string("lambdarate",r)=#

        #=
        len = length(proposedparams.lambdarates)
        for r=2:len
          sigma = 10.0 / len
          d1 = Truncated(Normal(currentparams.lambdarates[r], sigma), 0.0, Inf)
          proposedparams.lambdarates[r]  = rand(d1)
          d2 = Truncated(Normal(proposedparams.lambdarates[r], sigma), 0.0, Inf)
          propratio += logpdf(d2, currentparams.lambdarates[r]) - logpdf(d1, proposedparams.lambdarates[r])
          movename = string("lambdarate",r)
        end=#
        #=
        r1 = rand(1:len)
        r2 = rand(1:len)
        while r1 == r2
          r2 = rand(1:len)
        end
        rdelta = abs(randn()*0.01)
        proposedparams.lambdaweights[r1] += rdelta
        proposedparams.lambdaweights[r2] -= rdelta
        if proposedparams.lambdaweights[r1] >= 1.0 || proposedparams.lambdaweights[r2] < 0.0
          cont = false
        end=#

        #=
        len = length(proposedparams.lambdarates)
        r = rand(2:len)
        lb = 1.0
        ub = Inf
        if r > 1
          lb = currentparams.lambdarates[r-1]
        end
        if r < len
          ub = currentparams.lambdarates[r+1]
        end
        sigma = 4.0
        d1 = Truncated(Normal(currentparams.lambdarates[r], sigma), lb, ub)
        proposedparams.lambdarates[r]  = rand(d1)
        d2 = Truncated(Normal(proposedparams.lambdarates[r], sigma), lb, ub)
        propratio += logpdf(d2, currentparams.lambdarates[r]) - logpdf(d1, proposedparams.lambdarates[r])
        movename = string("lambdarate",r)=#
        #=
        proposedparams.lambdarates[r]  = currentparams.lambdarates[r] + randn(rng)
        if proposedparams.lambdarates[r] <= 1.0
          cont = false
        end
        for q=2:len
          if proposedparams.lambdarates[q-1] >= proposedparams.lambdarates[q]
            cont = false
            break
          end
        end=#
        #=
        d1 = Truncated(Normal(currentparams.q4, sigmas[4]), 0.0, Inf)
        proposedparams.q4  = rand(d1)
        d2 = Truncated(Normal(proposedparams.q4, sigmas[4]), 0.0, Inf)
        propratio += logpdf(d2, currentparams.q4) - logpdf(d1, proposedparams.q4)=#
      elseif move == 6
        d1 = Truncated(Normal(chain.currentparams.lambdazeroweight, lambdaZeroSigma), 0.0, 1.0)
        chain.proposedparams.lambdazeroweight  = rand(d1)
        d2 = Truncated(Normal(chain.proposedparams.lambdazeroweight, lambdaZeroSigma), 0.0, 1.0)
        propratio += logpdf(d2, chain.currentparams.lambdazeroweight) - logpdf(d1, chain.proposedparams.lambdazeroweight)

        #chain.proposedparams.lambdaweights, chain.proposedparams.lambdarates = discretizegamma2(chain.proposedparams.lambdazeroweight, chain.proposedparams.lambdaGammaShape, chain.proposedparams.lambdaGammaScale, chain.proposedparams.lambdaCats)
        chain.proposedparams.lambdaweightsGC, chain.proposedparams.lambdaratesGC, chain.proposedparams.lambdaweightsAT, chain.proposedparams.lambdaratesAT, chain.proposedparams.lambdaweightsGT, chain.proposedparams.lambdaratesGT = discretizegamma3(chain.proposedparams.lambdazeroweight, chain.proposedparams.lambdaGammaShapeGC, chain.proposedparams.lambdaGammaShapeAT, chain.proposedparams.lambdaGammaShapeGT, chain.proposedparams.lambdaGammaScaleGC, chain.proposedparams.lambdaCats)
        movename = "lambdazeroweight"
        #=
        d1 = Truncated(Normal(chain.currentparams.q5, sigmas[5]), 0.0, Inf)
        chain.proposedparams.q5  = rand(d1)
        d2 = Truncated(Normal(chain.proposedparams.q5, sigmas[5]), 0.0, Inf)
        propratio += logpdf(d2, chain.currentparams.q5) - logpdf(d1, chain.proposedparams.q5)=#
      elseif move == 7 # sample branch lengths
        #proposedparams.q6 = currentparams.q6 + randn()*sigmas[6]
        selection = shuffle([b for b=2:length(chain.currentparams.branchlengths)])[1:min(5,length(chain.currentparams.branchlengths))]

        for b in selection
          d1 = Truncated(Normal(chain.currentparams.branchlengths[b], max(0.001, sigmas[6+b-1])), 0.0, Inf)
          chain.proposedparams.branchlengths[b]  = rand(d1)
          d2 = Truncated(Normal(chain.proposedparams.branchlengths[b], max(0.001, sigmas[6+b-1])), 0.0, Inf)
          propratio += logpdf(d2, chain.currentparams.branchlengths[b]) - logpdf(d1, chain.proposedparams.branchlengths[b])
        end
      elseif move == 8
        r = rand(chain.rng)
        if r < 0.3
          mask = zeros(Int,14)
          #counting = Int[i for i=8:14]
          counting = Int[8,9,10,11,13]
          if lambdaCats == 1
            counting = Int[9,11,13]
          end
          shuffle!(counting)
          for i=1:2
            mask[counting[i]] = 1
          end

          newparams, covmat = gettunedparams(gettuningvector(chain.currentparams), tuningvectors, mask, jsondict["tuningScale"], maxTuningIter, covmat)
          settuningparams(newparams, chain.proposedparams)
          movename = "adaptivemcmc_lambdas2"
        elseif r < 0.6
          mask = ones(Int, 14)
          mask[12] = 0
          mask[14] = 0
          newparams, covmat = gettunedparams(gettuningvector(chain.currentparams), tuningvectors, ones(Int, 14), jsondict["tuningScale"], maxTuningIter, covmat)
          settuningparams(newparams, chain.proposedparams)
          movename = "adaptivemcmc_all"
        elseif r < 0.8
          mask = zeros(Int,14)
          mask[8] = 1
          mask[9] = 1
          mask[10] = 1
          mask[11] = 1
          #mask[12] = 1
          mask[13] = 1
          #mask[14] = 1

          if lambdaCats == 1
            mask[8] = 0
            mask[10] = 0
          end



          #mask[rand(chain.rng,9:14)] = 1

          newparams, covmat = gettunedparams(gettuningvector(chain.currentparams), tuningvectors, mask, jsondict["tuningScale"], maxTuningIter, covmat)
          settuningparams(newparams, chain.proposedparams)
          movename = "adaptivemcmc_lambdas"
        elseif r < 0.90
          mask = zeros(Int,14)
          mask[1] = 1
          mask[2] = 1
          mask[3] = 1
          mask[4] = 1
          mask[5] = 1
          mask[6] = 1
          mask[7] = 1

          newparams, covmat = gettunedparams(gettuningvector(chain.currentparams), tuningvectors, mask, jsondict["tuningScale"], maxTuningIter, covmat)
          settuningparams(newparams, chain.proposedparams)
          movename = "adaptivemcmc_Q"
        else
          mask = zeros(Int,14)
          counting = [i for i=1:6]
          shuffle!(counting)
          for i=1:3
            mask[counting[i]] = 1
          end

          newparams, covmat = gettunedparams(gettuningvector(chain.currentparams), tuningvectors, mask, jsondict["tuningScale"], maxTuningIter, covmat)
          settuningparams(newparams, chain.proposedparams)
          movename = "adaptivemcmc_Q3"
        end
      elseif move == 9
        mask = zeros(Int, length(chain.currentparams.branchlengths))
        counting = [i for i=1:length(mask)]
        shuffle!(counting)
        for i=1:25
          mask[counting[i]] = 1
        end

        chain.proposedparams.branchlengths = getbranchtunedparams(chain.currentparams.branchlengths, branchtuningvectors, mask, jsondict["tuningScale"], maxTuningIter)
        #println(chain.proposedparams.branchlengths)

        #=
        for i=1:length(mask)
          if mask[i] == 1
            println(chain.currentparams.branchlengths[i], "\t", chain.proposedparams.branchlengths[i])
          end
        end=#
        #println(chain.proposedparams.branchlengths)
        movename = "branchlengths"
      end



      



      cont = isvalid(chain.proposedparams)
      #=
      if length(chain.proposedparams.lambdaweightsGC) == 1
        chain.proposedparams.lambdaratesGC[1] = chain.proposedparams.lambdaGC-1.0
        chain.proposedparams.lambdaratesAT[1] = chain.proposedparams.lambdaAT-1.0
        chain.proposedparams.lambdaratesGT[1] = chain.proposedparams.lambdaGT-1.0
        chain.proposedparams.lambdaweightsGC[1] = 1.0
        chain.proposedparams.lambdaweightsAT[1] = 1.0
        chain.proposedparams.lambdaweightsGT[1] = 1.0
        #chain.proposedparams.lambdazeroweight = 1.0
      end=#
      paramsvector = getparamsvector(chain.proposedparams)
      chain.proposedparams = getparams(paramsvector, dataset, siteCats, lambdaCats, 1, fixGU, fixGCAU)
      if cont
        chain.proposedll = computetotallikelihood(chain.rng, dataset, chain.proposedparams, chain.paired, samplebranchlengths,integratesiterates,integratestructure,M,true,maxbasepairdistance, usecuda)
      end

      #=
      if move == 9
        println(cont, "\t", chain.proposedll, "\t", chain.currentll, "\t", (chain.proposedll-chain.currentll))
      end=#
      if samplebranchlengths
        chain.currentparams.q6 = 1.0
        chain.proposedparams.q6 = 1.0
      end
      if cont && (gibbs || exp((chain.proposedll-chain.currentll)*chain.B+propratio) > rand(chain.rng))
        chain.currentll = chain.proposedll
        chain.currentparams = deepcopy(chain.proposedparams)
        logAccept!(chain.logger,string(movename))
      else
        chain.proposedparams = deepcopy(chain.currentparams)
        chain.proposedll = chain.currentll
        logReject!(chain.logger,string(movename))
      end

    end

    if chain.B == 1.0 && iter < maxTuningIter
      push!(tuningvectors, gettuningvector(chain.currentparams))
      if samplebranchlengths
        push!(branchtuningvectors, chain.currentparams.branchlengths)
      end
    end

    if integratestructure || iter % 10 == 0 || (chain.B == 1.0 && iter % 5 == 0)
      if iter % 20 == 0
        if usemodelswitch
            if chain.B == 1.0
              U, q0, q1, delta = getmodelswitchlikelihoods(chain.rng, dataset, chain.proposedparams, chain.paired, samplebranchlengths,integratesiterates,integratestructure,M,true,maxbasepairdistance, usecuda)
              write(modelswitchlog, string(iter,"\t",M,"\t",U,"\t",q0,"\t",q1,"\t",delta,"\n"))
              flush(modelswitchlog)
            end
        end
      end

      prevtime = chain.timings[max(1,div(length(chain.timings),2))]
      iterspersecond = (iter-prevtime[1]+1.0) / (time() - prevtime[2])
      #sum((currentparams.lambdaratesGC+1.0).*currentparams.lambdaweightsGC),"\t",sum((currentparams.lambdaratesAT+1.0).*currentparams.lambdaweightsAT), "\t", sum((currentparams.lambdaratesGT+1.0).*currentparams.lambdaweightsGT)
      meanlambdaGC = sum((chain.currentparams.lambdaratesGC.+1.0).*chain.currentparams.lambdaweightsGC)
      meanlambdaAT = sum((chain.currentparams.lambdaratesAT.+1.0).*chain.currentparams.lambdaweightsAT)
      meanlambdaGT = sum((chain.currentparams.lambdaratesGT.+1.0).*chain.currentparams.lambdaweightsGT)
      #write(mcmclog, string(iter,"\t",chain.currentll, "\t", chain.proposedll,"\t", propratio,"\t", meanlambdaGC, "\t", meanlambdaAT, "\t", meanlambdaGT, "\t", chain.currentparams.lambdazeroweight, "\t", chain.currentparams.lambdaGammaShapeGC, "\t", chain.currentparams.lambdaGammaScaleGC, "\t", chain.currentparams.lambdaGammaShapeAT, "\t", chain.currentparams.lambdaGammaScaleAT, "\t", chain.currentparams.lambdaGammaShapeGT, "\t", chain.currentparams.lambdaGammaScaleGT, "\t", chain.currentparams.siteGammaShape, "\t", chain.currentparams.siteGammaScale, "\t", chain.currentparams.q1, "\t", chain.currentparams.q2, "\t", chain.currentparams.q3, "\t", chain.currentparams.q4, "\t", chain.currentparams.q5, "\t", chain.currentparams.q6,"\t"))
      write(mcmclog, string(iter,"\t",chain.currentll, "\t", chain.proposedll,"\t", propratio,"\t", chain.currentparams.lambdaGC, "\t", chain.currentparams.lambdaAT, "\t", chain.currentparams.lambdaGT, "\t", chain.currentparams.lambdazeroweight, "\t", chain.currentparams.lambdaGammaShapeGC, "\t", chain.currentparams.lambdaGammaScaleGC, "\t", chain.currentparams.lambdaGammaShapeAT, "\t", chain.currentparams.lambdaGammaScaleAT, "\t", chain.currentparams.lambdaGammaShapeGT, "\t", chain.currentparams.lambdaGammaScaleGT, "\t", chain.currentparams.siteGammaShape, "\t", chain.currentparams.siteGammaScale, "\t", chain.currentparams.q1, "\t", chain.currentparams.q2, "\t", chain.currentparams.q3, "\t", chain.currentparams.q4, "\t", chain.currentparams.q5, "\t", chain.currentparams.q6,"\t"))
      if samplebranchlengths
        for b=2:length(chain.currentparams.branchlengths)
          write(mcmclog, string(chain.currentparams.branchlengths[b],"\t"))
        end
      end
      for a=1:4
        write(mcmclog, string(chain.currentparams.freqs[a],"\t"))
      end

      for r=1:length(chain.currentparams.lambdaweightsGC)
        write(mcmclog, string(chain.currentparams.lambdaweightsGC[r],"\t"))
      end
      for r=1:length(chain.currentparams.lambdaratesGC)
        write(mcmclog, string(chain.currentparams.lambdaratesGC[r],"\t"))
      end
      for r=1:length(chain.currentparams.lambdaratesAT)
        write(mcmclog, string(chain.currentparams.lambdaratesAT[r],"\t"))
      end
      for r=1:length(chain.currentparams.lambdaratesGT)
        write(mcmclog, string(chain.currentparams.lambdaratesGT[r],"\t"))
      end
      for r=1:length(chain.currentparams.siteRates)
        write(mcmclog, string(chain.currentparams.siteRates[r],"\t"))
      end
      for r=1:length(chain.currentparams.siteWeights)
        write(mcmclog, string(chain.currentparams.siteWeights[r],"\t"))
      end
      write(mcmclog, string(iterspersecond, "\t"))
      write(mcmclog, string(chain.hashindex, "\n"))


      #try
        if M == 1.0 && chain.B == 1.0 && iter % 10 == 0
          if integratestructure
            unpairedlogprobs = computeunpairedlikelihoods(dataset, chain.currentparams)
            pairedlogprobs, ret = coevolutionall(dataset,chain.currentparams,true,true,true,maxbasepairdistance,usecuda)
            maskgapped!(pairedlogprobs,dataset.gapfrequency,0.5,-Inf)
            inside = computeinsideKH99(unpairedlogprobs, pairedlogprobs, 1.0,false,usecuda)
            Z = inside[1,1,dataset.numcols]
            savemaximum2(Z, getparamsvector2(chain.currentparams), string(outputprefix,".maxmcmc"))
            println("inside Z =", Z)
            chain.pairedlogprior = samplestructure(chain.rng, inside, pairedlogprobs, unpairedlogprobs, 1, dataset.numcols, chain.paired, grammar, chain.B)
          end

          write(structurelog, string(getdotbracketstring(chain.paired), "\n"))
          flush(structurelog)

          samplestates(chain.rng, dataset, chain.currentparams, chain.paired, chain.B)
          write(siterateslog,string(iter,"\t"))
          for col=1:dataset.numcols
            if chain.currentparams.states[col] > 0
              write(siterateslog,string(chain.currentparams.siteRates[chain.currentparams.states[col]],"\t"))
            end
          end
          write(siterateslog, "\n")
        end
      #catch
        #println("An error occured 1787")
      #end

      if M == 1.0 && chain.B == 1.0 && iter % 100 == 0
        #=
        unpairedlogprobs = computeunpairedlikelihoods(dataset, chain.currentparams)
        pairedlogprobs, ret = coevolutionall(dataset,chain.currentparams,true,true,true,maxbasepairdistance)
        maskgapped!(pairedlogprobs,dataset.gapfrequency,0.5,-Inf)
        inside = computeinsidecuda(unpairedlogprobs, pairedlogprobs, 1.0)
        outside = computeoutside(inside,unpairedlogprobs, pairedlogprobs, grammar)

        outside2 = computeoutsidecuda(inside,unpairedlogprobs, pairedlogprobs)

        unpairedposteriorprobs,pairedposteriorprobs = computebasepairprobs(inside, outside, unpairedlogprobs, pairedlogprobs, grammar)
        #unpairedpostlogprobs = log(unpairedposteriorprobs)
        pairedpostlogprobs = log(pairedposteriorprobs)

        basepairfile = string(outputprefixchain,".baseprobs")
        basepairlog = open(basepairfile, iter == 0 ? "w" : "a")
        write(basepairlog, string(">",iter,"\n"))
        write(basepairlog, string(pairedposteriorprobs,"\n"))
        close(basepairlog)


        meanlambdas = zeros(Float64, dataset.numcols, dataset.numcols)
        totals = ones(Float64, dataset.numcols, dataset.numcols)*-Inf
        for (probs, lambdarate, prior) in zip(ret,chain.currentparams.lambdarates,chain.currentparams.lambdaweights)
          logprior = log(prior)
          for p=1:dataset.numcols
            for q=p+1:dataset.numcols
              totals[p,q] = CommonUtils.logsumexp(totals[p,q], logprior+probs[p,q]+pairedpostlogprobs[p,q])
            end
          end
        end
        for (probs, lambdarate, prior) in zip(ret,chain.currentparams.lambdarates,chain.currentparams.lambdaweights)
          logprior = log(prior)
          for p=1:dataset.numcols
            for q=p+1:dataset.numcols
              weight = exp(logprior+probs[p,q]+pairedpostlogprobs[p,q]-totals[p,q])
              meanlambdas[p,q] += weight*lambdarate
              meanlambdas[q,p] += weight*lambdarate
            end
          end
        end
        lambdaslog = open(string(outputprefixchain,".lambdas"), iter == 0 ? "w" : "a")
        write(lambdaslog, string(">",iter,"\n"))
        write(lambdaslog, string(meanlambdas,"\n"))
        close(lambdaslog)

        posteriornotone = zeros(Float64, dataset.numcols, dataset.numcols)
        for p=1:dataset.numcols
          for q=p+1:dataset.numcols
            posteriornotone[p,q] = 1.0 - exp(log(chain.currentparams.lambdaweights[1])+ret[1][p,q]+pairedpostlogprobs[p,q]-totals[p,q])
            posteriornotone[q,p] = posteriornotone[p,q]
          end
        end
        posteriorfile = string(outputprefixchain,".posteriornotone")
        posteriorlog = open(posteriorfile, iter == 0 ? "w" : "a")
        write(posteriorlog, string(">",iter,"\n"))
        write(posteriorlog, string(posteriornotone,"\n"))
        close(posteriorlog)=#
      end

      if chain.B == 1.0
        write(tempacceptancelog, string(list(chain.logger),"\n"))
        flush(tempacceptancelog)
      end
    end
    if iter % 100 == 0
      flush(mcmclog)
      flush(siterateslog)
    end

    if iter % 500 == 0
      close(mcmclog)
      summarizemcmc(mcmclogfile)
      mcmclog = open(mcmclogfile,"a")
      write(acceptancelog, string(list(chain.logger),"\n"))
      flush(acceptancelog)
      MCMCLogging.clear!(chain.logger)
      #=
      if chain.B == 1.0 && iter % 500 == 0 && mode == MODE_IO_SAMPLE
        close(structurelog)
        structures = readdbnstructures(string(outputprefixchain, ".structures"))
        structurelog = open(string(outputprefixchain, ".structures"), chain.iter == 0 ? "w" : "a")
        if length(structures) > 3
          singleprobs, pairedprobs = getpairprobs(structures[max(1,div(length(structures),3)):end])
          consensus = getPosteriorDecodingConsensusStructure(pairedprobs, singleprobs)
          write(consensuslog, string(getdotbracketstring(consensus),"\n"))
          flush(consensuslog)
        end
      end=#
    end

  end
  close(mcmclog)
  close(siterateslog)
  close(consensuslog)
  close(structurelog)
  close(tempacceptancelog)
  if usemodelswitch
      close(modelswitchlog)
  end
  close(acceptancelog)
  chain.iter = chain.iter+runiters+1
  #return chain, burnins, burnindata, covmat
  return covmat
end

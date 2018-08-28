include("InsideCUDA.jl")

function iterativeposteriordecoding(singleprobs::Array{Float64,1}, pairprobs::Array{Float64,2}, alpha::Float64=2.0)
    datalen = length(singleprobs)
    ematrix = zeros(Float64, datalen, datalen)
    smatrix = zeros(Int, datalen, datalen)
    for i=1:datalen
        #ematrix[i,i] = singleprobs[i]
    end

    for diag=2:datalen
        for i=1:datalen-diag+1
            j = i + diag - 1
            #=
            if i+1 == j
              println("ERROR")
          end=#
            e1 = singleprobs[i] + ematrix[i+1,j]
            e2 = alpha*pairprobs[i,j] + ematrix[i+1,j-1]
            maxe3 = -Inf
            maxk = 1
            for k=i+1:j-1
                v = alpha*pairprobs[i,k] + ematrix[i+1,k-1] + ematrix[k+1,j]
                if v > maxe3
                    maxe3 = v
                    maxk = k
                end
            end
            maxval, maxindex = findmax([e1,e2,maxe3])
            ematrix[i,j] = maxval
            if maxindex == 1 || maxindex == 2
                smatrix[i,j] = -maxindex
            else
                smatrix[i,j] = maxk
            end
        end
    end
    return ematrix, smatrix
end

function posteriordecodingtraceback(smatrix::Array{Int,2})
  datalen = size(smatrix,1)
  paired = zeros(Int, datalen)
  stack = Tuple{Int,Int}[]
  push!(stack, (1, datalen))
  while length(stack) > 0
    i,j = pop!(stack)
    k = smatrix[i,j]
    if k == -1
        push!(stack, (i+1,j))
    elseif k == -2
      paired[i] = j
      paired[j] = i
      if i+1 < j-1
        push!(stack, (i+1,j-1))
      end
    elseif k > 0
      if i < k
        paired[i] = k
        paired[k] = i
      end
      if i+1 <= k-1
        push!(stack, (i+1,k-1))
      end
      if k+1 <= j
        push!(stack, (k+1,j))
      end
    end
  end

  return paired
end

function getPosteriorDecodingConsensusStructure(pairedposteriorprobs::Array{Float64,2}, unpairedposteriorprobs::Array{Float64,1}, usecuda::Bool=false, alpha::Float64=2.0)
  if usecuda
    ematrix, smatrix = computeposteriordecodingcuda(unpairedposteriorprobs, pairedposteriorprobs, alpha)
    return posteriordecodingtraceback(smatrix)
  else
    ematrix, smatrix = iterativeposteriordecoding(pairedposteriorprobs, unpairedposteriorprobs, alpha)
    return posteriordecodingtraceback(smatrix)
  end
end

function getpairprobs(structures::Array{Array{Int,1}})
  len = length(structures[1])
  singleprobs = zeros(Float64, len)
  pairedprobs = zeros(Float64, len, len)
  for structure in structures
    for i=1:length(structure)
      if structure[i] > i
        pairedprobs[i,structure[i]] += 1.0
        pairedprobs[structure[i],i] += 1.0
      elseif structure[i] == 0
        singleprobs[i] += 1.0
      end
    end
  end
  return singleprobs/length(structures), pairedprobs/length(structures)
end


function parsedbn(dbn::AbstractString)
  stack = Int[]
  stackChars = Int[]
  spl = strip(dbn)
  paired = zeros(Int, length(spl))
  left = ['(','<', 'A', 'B', 'C', 'D', 'E']
  right = [')','>', 'a', 'b', 'c', 'd', 'e']
  for i=1:length(spl)
    l = find(left .== spl[i])
    r = find(right .== spl[i])
    if length(l) == 1
      push!(stack,i)
      push!(stackChars, l[1])
    elseif length(r) == 1
      if length(stackChars) > 0 && stackChars[end] == r[1]
        j = pop!(stack)
        pop!(stackChars)
        paired[i] = j
        paired[j] = i
      end
    end
  end

  #=
  for i=1:length(spl)
    if spl[i] == '('
      push!(stack,i)
    elseif spl[i] == ')'
      j = pop!(stack)
      paired[i] = j
      paired[j] = i
    end
  end=#
  return paired
end

function removelonelybasepairs(paired::Array{Int,1})
  ret = copy(paired)
  for i=1:length(paired)
    if paired[i] > i && (i==1 || paired[i-1] == 0) && (i == length(paired) || paired[i+1] == 0)
      ret[paired[i]] = 0
      ret[i] = 0
    end
  end
  return ret
end

function readdbnstructures(infile)
  structures = Array{Int,1}[]
  fin = open(infile, "r")
  for line in readlines(fin)
    if length(strip(line)) > 0
      push!(structures, parsedbn(line))
    end
  end
  close(fin)
  return structures
end

function readctfile(filename::String)
  f = open(filename)
  lines = readlines(f)
  close(f)

  spl = split(lines[1])
  len = parse(Int, spl[1])
  paired = zeros(Int,len)
  sequence = ""
  for i=2:length(lines)
    spl = split(lines[i])
    if length(spl) >= 5
      x = parse(Int,spl[1])
      c = spl[2]
      sequence = string(sequence, c)
      y = parse(Int,spl[5])
      if y > x
        if paired[x] != 0 || paired[y] != 0
          #println("ZZZ",x,"\t",y,"\t",paired[x],"\t",paired[y])
        else
          paired[x] = y
          paired[y] = x
        end
      end
    end
  end

  return sequence, paired
end

function precision(real::Array{Int,1}, predicted::Array{Int,1})
  count = 0.0
  total = 0.0
  for i=1:length(predicted)
    if predicted[i] > i
      total += 1.0
      if real[i] == predicted[i]
        count += 1.0
      end
    end
  end
  return count/total
end

function isconflict(real::Array{Int,1}, x::Int, y::Int)
  conflict = false
  for i=1:length(real)
    if real[i] > i
      if i == x && real[i] == y

      elseif i < x < real[i] && i < y < real[i]

      elseif x < i && y > real[i]

      else
        conflict = true
      end
    end
  end

  return conflict
end

function mcc(real::Array{Int,1}, predicted::Array{Int,1})
  TP = 0.0
  FN = 0.0
  FP = 0.0
  TN = 0.0
  E = 0.0
  len = length(real)
  for i=1:length(predicted)
    if predicted[i] > i && real[i] == predicted[i]
      TP += 1.0
    end

    if predicted[i] > i && real[i] != predicted[i]
      if !isconflict(real, i, predicted[i])
        E += 1.0
      end
      FP += 1.0
    end

    if real[i] > i &&  real[i] != predicted[i]
      FN += 1.0
    end
  end
  TN = (len*len-len)/2.0  - TP - FN - FP

  FP = FP - E

  num = TP*TN - FP*FN
  den = sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  return num/den
end

function mcchelper(real::Array{Int,1}, predicted::Array{Int,1})
  TP = 0.0
  FN = 0.0
  FP = 0.0
  TN = 0.0
  E = 0.0
  len = length(real)
  for i=1:length(predicted)
    if predicted[i] > i && real[i] == predicted[i]
      TP += 1.0
    end

    if predicted[i] > i && real[i] != predicted[i]
      if !isconflict(real, i, predicted[i])
        E += 1.0
      end
      FP += 1.0
    end

    if real[i] > i &&  real[i] != predicted[i]
      FN += 1.0
    end
  end
  TN = (len*len-len)/2.0  - TP - FN - FP

  FP = FP - E

  num = TP*TN - FP*FN
  den = sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  return TP, FN, FP, TN, E
end

function recall(real::Array{Int,1}, predicted::Array{Int,1})
  return precision(predicted, real)
end

function f1score(real::Array{Int,1}, predicted::Array{Int,1})
  p = precision(real, predicted)
  r = recall(real, predicted)

  return 2.0*p*r/(p+r)
end

function weightedmountainvector(paired::Array{Int,1})
  f1 = zeros(Float64, length(paired))
  for i=1:length(paired)
    if i > 1
      f1[i] = f1[i-1]
    end
    if paired[i] != 0
      f1[i] += 1.0 / (paired[i] - i)
    end
  end
  return f1
end

function mountainvector(paired::Array{Int,1})
  f1 = zeros(Float64, length(paired))
  for i=1:length(paired)
    if i > 1
      f1[i] = f1[i-1]
    end
    if paired[i] != 0
      if paired[i] > i
        f1[i] += 1.0
      else
        f1[i] -= 1.0
      end
    end
  end
  return f1
end


function structurestar(len::Int)
  dbn = ""
  if len % 2 == 0
    dbn = string(repeat("(", div(len-2, 2)), "..", repeat(")", div(len-2, 2)))
  else
    dbn = string(repeat("(", div(len-1, 2)), ".", repeat(")", div(len-1, 2)))
  end
  return parsedbn(dbn)
end

function structurezero(len::Int)
  return zeros(Int,len)
end

function weightedmountaindistance(paired1::Array{Int,1}, paired2::Array{Int,1})
  f1 = weightedmountainvector(paired1)
  f2 = weightedmountainvector(paired2)
  d = 0.0
  for i=1:length(f1)
    d += abs(f1[i] - f2[i])
  end
  return d
end

function weightedmountaindiameter(len::Int)
  return weightedmountaindistance(structurestar(len), structurezero(len))
end

function mountaindistance(paired1::Array{Int,1}, paired2::Array{Int,1}, p::Float64)
  f1 = mountainvector(paired1)
  f2 = mountainvector(paired2)
  d = 0.0
  for i=1:length(f1)
    d += abs(f1[i] - f2[i])^p
  end
  return d
end

function mountaindiameter(len::Int, p::Float64)
  return mountaindistance(structurestar(len), structurezero(len), p)
end

function normalisedmountaindistance(paired1::Array{Int,1}, paired2::Array{Int,1}, p::Float64)
  return mountaindistance(paired1,paired2,p) / mountaindiameter(length(paired1), p)
end

function normalisedweightedmountaindistance(paired1::Array{Int,1}, paired2::Array{Int,1})
  return weightedmountaindistance(paired1,paired2) / weightedmountaindiameter(length(paired1))
end

function readdbnfile(filename::String)
  f = open(filename)
  lines = readlines(f)
  close(f)
  return parsedbn(lines[1])
end

function writectfile(paired::Array{Int,1}, seq::String, filename::String)
  f = open(filename, "w")
  write(f, string(length(paired),"\n"))
  for i=1:length(paired)
    write(f,string(i,"\t",seq[i], "\t", i-1,"\t", i+1,"\t",paired[i],"\t",i,"\n"))
  end
  close(f)
end

function runppfold(alignmentfile, treefile, outputdir)
  readstring(`java -jar ../binaries/PPfold3.1.1.jar $alignmentfile -t $treefile --onlyCT -o $outputdir`)
  m = match(r"^(.*/)([^/]*)$", alignmentfile)
  #outputdir = m[1]
  f = string(match(r"([^\\]*\.)(\w+)$",m[2])[1], "ct")
  return readctfile(string(outputdir, f))
end

function runrnaalifold(alignmentfile, treefile, outputdir)
  st = readstring(`RNAalifold $alignmentfile`)
  return parsedbn(split(st)[2])
  #println(st)
  #=
  m = match(r"^(.*/)([^/]*)$", alignmentfile)
  #outputdir = m[1]
  f = string(match(r"([^\\]*\.)(\w+)$",m[2])[1], "ct")
  return readctfile(string(outputdir, f))=#
end

function getpairprobs(structures::Array{Array{Int,1}})
  len = length(structures[1])
  singleprobs = zeros(Float64, len)
  pairedprobs = zeros(Float64, len, len)
  for structure in structures
    for i=1:length(structure)
      if structure[i] > i
        pairedprobs[i,structure[i]] += 1.0
        pairedprobs[structure[i],i] += 1.0
      elseif structure[i] == 0
        singleprobs[i] += 1.0
      end
    end
  end
  return singleprobs/length(structures), pairedprobs/length(structures)
end

function maskgapped!(probs::Array{Float64,2}, gapfrequency::Array{Float64,1}, cutoff::Float64, v::Float64=0.0)
  numcols = length(gapfrequency)
  for i=1:numcols
    for j=1:numcols
      if gapfrequency[i] >= cutoff || gapfrequency[j] >= cutoff
        probs[i,j] = v
      end
    end
  end
  return probs
end

function maskgapped(paired::Array{Int,1}, gapfrequency::Array{Float64,1}, cutoff::Float64, v::Int=0)
  structure = copy(paired)
  numcols = length(gapfrequency)
  for i=1:numcols
    if structure[i] > i
      j = structure[i]
      if gapfrequency[i] >= cutoff || gapfrequency[j] >= cutoff
        structure[i] = 0
        structure[j] = 0
      end
    end
  end
  return structure
end

function maskgapped!(structures::Array{Array{Int,1},1}, gapfrequency::Array{Float64,1}, cutoff::Float64, v::Int=0)
  numcols = length(gapfrequency)
  for structure in structures
    for i=1:numcols
      if structure[i] > i
        j = structure[i]
        if gapfrequency[i] >= cutoff || gapfrequency[j] >= cutoff
          structure[i] = 0
          structure[j] = 0
        end
      end
    end
  end
end


function samplethermodynamic(thermodynamicsamples, rng, sequences)
  index = rand(1:length(sequences))
  v = get(thermodynamicsamples, index, [])
  if length(v) == 0
    (so,si,pr) = readandwrite(`RNAsubopt -p 500`)
    write(si,string(sequences[index],"\n"))
    write(si, "@\n")
    spl = split(readstring(so))
    close(si)

    for i=2:length(spl)
      sample = parsedbn(strip(spl[i]))
      push!(v, sample)
    end
    thermodynamicsamples[index] = v
  end
  len = length(v)
  return splice!(v,rand(1:len))
end

#include("Combinatorics.jl")

#=
real = parsedbn(".(((.((...))..)))")
predicted = parsedbn("((...((...))...))")
TP, FN, FP, TN, E = mcchelper(real,predicted)
println(real)
println(predicted)
println((TP, FN, FP, TN, E))
println(precision(real,predicted))
println(recall(real,predicted))
println(f1score(real,predicted))
println(mcc(real,predicted))
println(mountaindistance(real,predicted,1.0))
println(mountaindiameter(length(real),1.0))
println(normalisedmountaindistance(real,predicted,1.0))=#

using FastaIO


function getdinuccode(c, d)
  nuc1 = get(nucmapping, c, 0)
  nuc2 = get(nucmapping, d, 0)
  return (nuc1-1)*4 + nuc2
end

canonicalGCAU = [getdinuccode('G','C'), getdinuccode('C','G'), getdinuccode('A','U'), getdinuccode('U','A')]
canonicalGCAUGU = [getdinuccode('G','C'), getdinuccode('C','G'), getdinuccode('A','U'), getdinuccode('U','A'), getdinuccode('G','U'), getdinuccode('U','G')]
#pairs = canonicalGCAUGU

function computeMutualInformation(sequences::Array{String,1}, dataset::Dataset, pos1::Int, pos2::Int, pairs)
  freqs = zeros(Float64, 4)
  numgaps = 0.0
  watsonCount = 0.0
  pairTotal = 0.0
  for s=1:dataset.numseqs
    nuc1 = get(nucmapping, sequences[s][pos1],0)
    if nuc1 > 0
      freqs[nuc1] += 1.0
    end
    nuc2 = get(nucmapping, sequences[s][pos2],0)
    if nuc2 > 0
      freqs[nuc2] += 1.0
    end

    dinuc = (nuc1-1)*4 + nuc2
    if nuc1 == 0 || nuc2 == 0
      numgaps += 1.0
    elseif dinuc in canonicalGCAUGU
      watsonCount += 1.0
    end
    pairTotal += 1.0
  end
  freqs /= sum(freqs)

  observedWatson = watsonCount / pairTotal
  expectedWatson = 0.0
  for i=1:4
    for j=1:4
      dinuc = (i-1)*4 + j
      if dinuc in pairs
        expectedWatson += freqs[i] * freqs[j]
      end
    end
  end

  return (observedWatson * ((observedWatson == 0.0 ? 0.0 : log(observedWatson / expectedWatson)) / log(2.0))) - (numgaps / pairTotal)
end

function computeRNAalifoldMutualInformation(sequences::Array{String,1}, dataset::Dataset, pos1::Int, pos2::Int, pairs)
  C = 0.0
  N2 = 0.0
  q = 0.0
  N = 0.0

  numseqs = length(sequences)
  dinucs = Int[]
  g = zeros(Float64, numseqs)
  for a=1:numseqs
    dinuc = getdinuccode(sequences[a][pos1], sequences[a][pos2])
    push!(dinucs, dinuc)
    if getdinuccode(sequences[a][pos1], sequences[a][pos2]) in pairs
     g[a] = 1.0
    end
  end

  for a=1:numseqs
    for b=a+1:numseqs
      hamming = 0.0
      if dinucs[a] == dinucs[b]
        hamming = 0.0
      elseif sequences[a][pos1] == sequences[b][pos1] || sequences[a][pos2] == sequences[b][pos2]
        hamming = 1.0
      else
        hamming = 2.0
      end

      C += hamming * g[a] * g[b]
      N2 += 1.0
    end
    q += 1.0 - g[a]
    N += 1.0
  end

  C /= N2
  q /= N

  return C-q
end

function computeObservedFreqsWithStacking(sequences::Array{String,1}, pos1::Int, pos2::Int, pos3::Int, pos4::Int, stacking::Bool, pairs)
    watsonCount = 0.0
    pairTotal = 0.0
    freqs = zeros(Float64, 4)
    numGaps = 0.0
    numseqs = length(sequences)
    for i=1:numseqs
        dinuc = getdinuccode(sequences[i][pos1], sequences[i][pos2])
        dinucStacking = getdinuccode(sequences[i][pos3], sequences[i][pos4])
        if dinuc in pairs && (stacking == (dinucStacking in pairs))
            watsonCount += 1
        end
        pairTotal += 1
    end

    return watsonCount / pairTotal
end

function computeMutualInformationWithStacking(sequences::Array{String,1}, pos1::Int, pos2::Int, pos3::Int, pos4::Int, pairs)
  freqs = zeros(Float64, 4)
  pairTotal = 0.0
  numGaps = 0.0
  numStackingGaps = 0.0
  numseqs = length(sequences)
  for i=1:numseqs
      dinuc = getdinuccode(sequences[i][pos1], sequences[i][pos2])
      dinucStacking = getdinuccode(sequences[i][pos3], sequences[i][pos4])
      if dinuc <= 0 || dinucStacking <= 0
          numStackingGaps += 1
      elseif dinuc <= 0
          numGaps += 1.0
      end
      pairTotal += 1.0

      nuc1 = get(nucmapping, sequences[i][pos1],0)
      if nuc1 > 0
        freqs[nuc1] += 1.0
      end
      nuc2 = get(nucmapping, sequences[i][pos2],0)
      if nuc2 > 0
        freqs[nuc2] += 1.0
      end
  end
  freqs /= sum(freqs)

  q = 0.0
  for i=1:4
    for j=1:4
      dinuc = (i-1)*4 + j
      if dinuc in pairs
        q += freqs[i] * freqs[j]
      end
    end
  end

  s = 0.0
  for y=1:2
      p = computeObservedFreqsWithStacking(sequences, pos1, pos2, pos3, pos4, y != 2, pairs)
      if (p != 0.0)
          s += p * log(p / q)
      end
  end

  return s - ((numGaps + numStackingGaps) / (pairTotal * 2.0))
end

function computeCovariation(paired, fastafile, newickfile)
  dataset = Dataset(fastafile, newickfile)
  names =  String[]
  sequences = String[]
  FastaIO.FastaReader(fastafile) do fr
   for (desc, seq) in fr
     len = length(seq)
     push!(names,desc)
     push!(sequences, seq)
   end
  end

  len = length(paired)
  mutualinformation = zeros(Float64, len)
  rnaalifold_mutualinformation = zeros(Float64, len)
  mutualinformation_stacking = zeros(Float64, len)
  for i=1:len
    if paired[i] > i
      mutualinformation[i] = computeMutualInformation(sequences,dataset,i,paired[i], canonicalGCAUGU)
      rnaalifold_mutualinformation[i] = computeRNAalifoldMutualInformation(sequences,dataset,i,paired[i], canonicalGCAUGU)
      mutualinformation_stacking[i] = computeMutualInformationWithStacking(sequences, i, paired[i], i+1, paired[i]-1, canonicalGCAUGU)
    end
  end

    return mutualinformation, rnaalifold_mutualinformation, mutualinformation_stacking
end

#=
fastafile = "/media/michael/Sandisk500GB/data/hiv1.fas.norm50.structure.align"
newickfile = "/media/michael/Sandisk500GB/data/hiv1_4/hiv.nwk"
sequence, pairedsites = readctfile("/media/michael/Sandisk500GB/data/hiv1.fas.norm50.structure.align.ct")
mutualinformation, rnaalifold_mutualinformation, mutualinformation_stacking = computeCovariation(pairedsites, fastafile, newickfile)
i = 1
for (a,b,c) in zip(mutualinformation, rnaalifold_mutualinformation, mutualinformation_stacking)
  if pairedsites[i] > i
    println(i, "\t", pairedsites[i], "\t", a,"\t", b, "\t", c)
  end
  i += 1
end=#

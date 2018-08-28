using FastaIO
using Random
using SharedArrays

RULETYPE1 = 1 # A->. (aka U->.)
RULETYPE2 = 2 # A -> (B) (aka U->(V))
RULETYPE3 = 3 # A->. A->BC (aka U->VW)

mutable struct Rule
  left::Char
  right::String
  prob::Float64
  logprob::Float64
  ruletype::Int
  leftindex::Int
  rightindices::Array{Int,1}

  function Rule(left::Char,right::String,prob::Float64,ruletype::Int)
    ruleindex = Dict('S' => 1, 'F' => 2, 'L' => 3)
    leftindex = ruleindex[left]
    rightindices = Int[]
    for rchar in right
      push!(rightindices, get(ruleindex,rchar,0))
    end
    new(left,right,prob, log(prob),ruletype,leftindex,rightindices)
  end
end

export KH99
mutable struct KH99
  rules::Array{Rule,1}
  type1rules::Array{Rule,1}
  type2rules::Array{Rule,1}
  type3rules::Array{Rule,1}
  ruleindex::Dict{Char,Int}

  export KH99
  function KH99()
    rules = Rule[]
    push!(rules, Rule('S', "LS",0.868534, 3))
    push!(rules, Rule('S', "s",0.117609877998, 1))
    push!(rules, Rule('S', "dFd",0.013856122002, 2))
    push!(rules, Rule('F', "dFd",0.787640, 2))
    push!(rules, Rule('F', "LS",0.21236, 3))
    push!(rules, Rule('L', "s",0.894603, 1))
    push!(rules, Rule('L', "dFd",0.105397, 2))
    type1rules = Rule[rules[2],rules[6]]
    type2rules = Rule[rules[3],rules[4],rules[7]]
    type3rules = Rule[rules[1],rules[5]]
    ruleindex = Dict('S' => 1, 'F' => 2, 'L' => 3)
    return new(rules, type1rules, type2rules, type3rules, ruleindex)
  end
end

nucmapping = Dict('A' => 1, 'C' => 2, 'G' => 3, 'T' => 4, 'U' => 4)
function safelog(x::Float64)
  if x < 0.0
    return -1e10
  else
    return log(x)
  end
end


function orderpair(a::Float64, b::Float64)
  if a < b
    return a,b
  end
  return b,a
end



function quickExp(v::Float64)
  if v < -57.0
    return 0.0
  elseif v == 0.0
    return 1.0
  end

  return exp(v)
end

function calculateKH99prior(pairedsites::Array{Int,1})
  rules = Rule[]
  push!(rules, Rule('S', "sS",0.868534*0.894603, 0))
  push!(rules, Rule('S', "dFdS",0.868534*0.105397, 0))
  push!(rules, Rule('S', "s",0.117609877998, 0))
  push!(rules, Rule('S', "dFd",0.013856122002, 0))
  push!(rules, Rule('F', "dFd",0.787640, 0))
  push!(rules, Rule('F', "sS",0.21236*0.894603, 0))
  push!(rules, Rule('F', "dFdS",0.21236*0.105397, 0))
  push!(rules, Rule('L', "s",0.894603, 0))
  push!(rules, Rule('L', "dFd",0.105397, 0))

  s = 1
  e = length(pairedsites)
  stack = [('S',s,e)]
  ll = 0.0
  #println(pairedsites)
  while length(stack) > 0
    elem = pop!(stack)
    symbol = elem[1]
    s = elem[2]
    e = elem[3]
    #println(">",elem)
    if pairedsites[s] == 0
      for r in rules
        if r.left == symbol && r.right[1] == 's'
          if s == e && length(r.right) == 1
            ll += r.logprob
          elseif s < e && length(r.right) == 2
            ll += r.logprob
            push!(stack, (r.right[2],s+1,e))
            #println("A",stack[end])
            break
          end
        end
      end
    elseif pairedsites[s] > s
      for r in rules
        if r.left == symbol && r.right[1] == 'd'
          if pairedsites[s] == e && length(r.right) == 3
            ll += r.logprob
            push!(stack, (r.right[2],s+1,pairedsites[s]-1))
            #println("B",stack[end])
            break
          elseif pairedsites[s] != e && length(r.right) > 3
            ll += r.logprob
            push!(stack, (r.right[2],s+1,pairedsites[s]-1))
            #println("C",stack[end])
            push!(stack, (r.right[4],pairedsites[s]+1,e))
            #println("D",stack[end])
            break
          end
        end
      end
    end
  end

  return ll
end

function readpairedlogprobs(infile)
  m = 0
  rin = open(infile,"r")
  for line in eachline(rin)
    spl = split(line,",")
    x = parse(Int,spl[1])
    y = parse(Int,spl[2])
    m = max(max(m,x),y)
    p = parse(Float64,spl[3])
  end
  close(rin)
  pairedlogprobs = ones(Float64,m+1,m+1)*-Inf
  m = 0
  rin = open(infile,"r")
  for line in eachline(rin)
    spl = split(line,",")
    x = parse(Int,spl[1])
    y = parse(Int,spl[2])
    m = max(max(m,x),y)
    p = parse(Float64,spl[3])
    pairedlogprobs[x+1,y+1] = p
    pairedlogprobs[y+1,x+1] = p
  end
  close(rin)

  return pairedlogprobs
end

function getrnaprobs(fastafile)
  unpairedprior = Float64[0.337 0.207 0.202 0.254]
  pairedprior = Float64[0.007 0.004 0.011 0.127; 0.006 0.005 0.283 0.005; 0.023 0.275 0.006 0.045; 0.132 0.004 0.060 0.001]

  sequences = AbstractString[]
  len = 0
  FastaReader(fastafile) do fr
    for (desc, seq) in fr
      len = length(seq)
      push!(sequences, seq)
    end
  end

  unpairedlogprobs = zeros(Float64, len)
  pairedlogprobs = zeros(Float64, len, len)
  for seq in sequences
    for i=1:len
      nuc1 = get(nucmapping, seq[i], 0)
      if nuc1 > 0
        unpairedlogprobs[i] += log(unpairedprior[nuc1])
        pairedlogprobs[i,i] = -Inf
      end
      for j=i+1:len
        nuc2 = get(nucmapping, seq[j], 0)
        if nuc1 > 0 && nuc2 > 0
          pairedlogprobs[i,j] += log(pairedprior[nuc1,nuc2])
        elseif nuc1 > 0
          pairedlogprobs[i,j] += log(unpairedprior[nuc1])
        elseif nuc2 > 0
          pairedlogprobs[i,j] += log(unpairedprior[nuc2])
        end
        pairedlogprobs[j,i] = pairedlogprobs[i,j]
      end
    end
  end
  for i=1:len
    for j=max(1,i-2):min(len,i+2)
      pairedlogprobs[i,j] = -Inf
    end
  end
  return unpairedlogprobs, pairedlogprobs
end

function calculatedinucfreqs(obsFreqs::Array{Float64,1}, lambdaGC::Float64, lambdaAT::Float64, lambdaGT::Float64)
  dinucfreqs = zeros(Float64,16)

  piGC = obsFreqs[2]*obsFreqs[3]
  piAT = obsFreqs[1]*obsFreqs[4]
  piGT = obsFreqs[3]*obsFreqs[4]

  kappa = 1.0/(1.0 + 2.0*(piAT*(lambdaAT/lambdaATinv-1.0)) + 2.0*(piGC*(lambdaGC/lambdaGCinv-1.0)) + 2.0*(piGT*(lambdaGT/lambdaGTinv-1.0)))

  basepairings = [0,0,0,2,0,0,1,0,0,1,0,3,2,0,3,0]

  for h=1:4
    for v=1:4
      idx = (h-1)*4+v
      if basepairings[idx] == 1
        dinucfreqs[idx] = kappa*lambdaGC*lambdaGC*obsFreqs[h]*obsFreqs[v]
      elseif basepairings[idx] == 2
        dinucfreqs[idx] = kappa*lambdaAT*lambdaAT*obsFreqs[h]*obsFreqs[v]
      elseif basepairings[idx] == 3
        dinucfreqs[idx] = kappa*lambdaGT*lambdaGT*obsFreqs[h]*obsFreqs[v]
      else
        dinucfreqs[idx] = kappa*obsFreqs[h]*obsFreqs[v]
      end
    end
  end
  return dinucfreqs
end

function calculatenonevolutionaryprobs(fastafile, unpairedlogprobs::Array{Float64,1}, pairedlogprobs::Array{Float64,2}, obsfreqs::Array{Float64,1}, lambdaGC::Float64, lambdaAT::Float64, lambdaGT::Float64)
  sequences = AbstractString[]
  len = 0
  FastaReader(fastafile) do fr
    for (desc, seq) in fr
      len = length(seq)
      push!(sequences, seq)
    end
  end
  data = zeros(Int,length(sequences),len)
  for s=1:length(sequences)
    seq = sequences[s]
    for j=1:len
      data[s,j] = get(nucmapping,seq[j],0)
    end
  end

end

function readpairedlogprobs(pairedlogprobs::Array{Float64,2},infile)
  m = 0
  rin = open(infile,"r")
  for line in eachline(rin)
    spl = split(line,",")
    x = parse(Int,spl[1])
    y = parse(Int,spl[2])
    m = max(max(m,x),y)
    p = parse(Float64,spl[3])
    pairedlogprobs[x+1,y+1] = p
    pairedlogprobs[y+1,x+1] = p
  end
  close(rin)
end

function readunpairedlogprobs(infile)
  rin = open(infile,"r")
  unpairedlogprobs = Float64[]
  for line in eachline(rin)
    spl = split(line,",")
    x = parse(Int,spl[1])
    r = parse(Int,spl[2])
    p = parse(Float64,spl[3])
    push!(unpairedlogprobs,p)
  end
  return unpairedlogprobs
end

function computeinside(unpairedlogprobs::Array{Float64,1}, pairedlogprobs::Array{Float64,2}, grammar::KH99, B::Float64=1.0)
  type1rules = grammar.type1rules
  type2rules = grammar.type2rules
  type3rules = grammar.type3rules
  numnonterminals = 3
  len = length(unpairedlogprobs)
  inside = ones(Float64, numnonterminals, len, len)*-Inf
  for type1rule in type1rules
    for i=1:len
      inside[type1rule.leftindex, i,i] = (unpairedlogprobs[i] + type1rule.logprob)*B
    end
  end

  tmpvec = zeros(Float64,len)
  #order = zeros(Int,len,len)
  #c = 1
  for b=2:len
    for j=1:len-b+1
      for type3rule in type3rules
        #
        #tmp = -Inf
        for h=j:j+b-2
          prob1 = inside[type3rule.rightindices[1],j,h]
          prob2 = inside[type3rule.rightindices[2],h+1,j+b-1]
          tmpvec[h] = prob1+prob2
          #tmp = logsumexp(tmp, prob1 + prob2)
        end
        tmp = logsumexp(tmpvec,j,j+b-2)
        #=
        tmp = @parallel (logsumexp) for h=j:j+b-2
        tmpvec[h]
      end=#
      inside[type3rule.leftindex,j,j+b-1] = logsumexp(inside[type3rule.leftindex,j,j+b-1], type3rule.logprob*B+tmp)
    end

    for type2rule in type2rules
      inside[type2rule.leftindex,j,j+b-1] = logsumexp(inside[type2rule.leftindex,j,j+b-1], (type2rule.logprob + pairedlogprobs[j,j+b-1])*B + inside[type2rule.rightindices[2],j+1,j+b-2])
    end
    #order[j,j+b-1] = c
    #c += 1
    #println("£\t",j, "\t", b, "\t",j+b-1)
  end
end

return inside
end

function computinsideparallel2(unpairedlogprobs::Array{Float64,1}, pairedlogprobs::Array{Float64,2}, grammar::KH99, insidetemp::Array{Float64,3}, B::Float64=1.0)
  numnonterminals = 3
  len = length(unpairedlogprobs)
  inside = ones(Float64, numnonterminals, len, len)*-Inf
  type1rules = grammar.type1rules
  for type1rule in type1rules
    for i=1:len
      inside[type1rule.leftindex, i,i] = (unpairedlogprobs[i] + type1rule.logprob)*B
    end
  end

  #=
  blen = 100
  blocks = Tuple{Int,Int,Int,Int}[]
  b = 2
  while b < len
  bstart = 1
  bend = b + blen - 1

  j = 1
  while j < len
  jstart = j
  jend = jstart + blen - 1
  block = (jstart,min(len,jend),max(2,bstart),min(len,bend-1))
  println(block)
  push!(blocks, block)
  j += blen
end
b += blen
end=#

blen = 100
blocks = Tuple{Int,Int,Int,Int}[]
b = 1
while b < len
  bstart = b
  bend = b + blen - 1

  j = 1
  while j < len
    jstart = j
    jend = jstart + blen - 1
    block = (jstart,min(len,jend),max(2,bstart),min(len,bend))
    #println(block)
    push!(blocks, block)
    j += blen
  end
  b += blen
end

c = 1
for (jstart,jend,bstart,bend) in blocks
  if c>=6
    bstart = 1
    #bend = blen
  end
  println("Q",jstart,"\t",jend,"\t",bstart,"\t",bend)
  computeinsideblock(inside, unpairedlogprobs, pairedlogprobs, grammar, jstart, jend, bstart, bend, insidetemp, B)
  c += 1
end

return inside
end


function computeinsideblock(inside::Array{Float64,3}, unpairedlogprobs::Array{Float64,1}, pairedlogprobs::Array{Float64,2}, grammar::KH99, jstart::Int,jend::Int,bstart::Int, bend::Int, insidetemp::Array{Float64,3}, B::Float64=1.0)
  type1rules = grammar.type1rules
  type2rules = grammar.type2rules
  type3rules = grammar.type3rules
  len = length(unpairedlogprobs)

  tmpvec = zeros(Float64,len)

  for b=bstart:bend
    for j=jstart:jend
      #inside[:,j,j+b-1] = -Inf
      for type3rule in type3rules
        for h=j:j+b-2
          prob1 = inside[type3rule.rightindices[1],j,h]
          prob2 = inside[type3rule.rightindices[2],h+1,j+b-1]
          tmpvec[h] = prob1+prob2
        end
        tmp = logsumexp(tmpvec,j,j+b-2)
        inside[type3rule.leftindex,j,j+b-1] = logsumexp(inside[type3rule.leftindex,j,j+b-1], type3rule.logprob*B+tmp)
      end

      for type2rule in type2rules
        inside[type2rule.leftindex,j,j+b-1] = logsumexp(inside[type2rule.leftindex,j,j+b-1], (type2rule.logprob + pairedlogprobs[j,j+b-1])*B + inside[type2rule.rightindices[2],j+1,j+b-2])
      end

      for symbol=1:3
        if abs(insidetemp[symbol,j,j+b-1] - inside[symbol,j,j+b-1]) > 0.0001
          println("%\t",jstart,"\t",jend,"\t",bstart,"\t",bend)
          println(insidetemp[symbol,j,j+b-1],"\t",inside[symbol,j,j+b-1])
        end
      end
    end
  end
  return inside
end

function computeoutside(inside::Array{Float64,3}, unpairedlogprobs::Array{Float64,1}, pairedlogprobs::Array{Float64,2}, grammar::KH99)
  type1rules = grammar.type1rules
  type2rules = grammar.type2rules
  type3rules = grammar.type3rules
  numnonterminals = 3
  len = length(unpairedlogprobs)
  outside = ones(Float64, numnonterminals, len, len)*-Inf
  outside[1, 1, len] = 0.0
  for c=2:len
    b = len - c
    for j=0:len-b-1
      for type3rule in type3rules
        tmp = -Inf
        for k=j+b+1:len-1
          oprob = outside[type3rule.leftindex,j+1,k+1]
          iprob = inside[type3rule.rightindices[2],j+b+1+1,k+1]
          tmp = logsumexp(tmp, oprob+iprob)
        end
        outside[type3rule.rightindices[1], j+1, j+b+1] = logsumexp(outside[type3rule.rightindices[1], j+1, j+b+1], tmp + type3rule.logprob)
        tmp = -Inf
        for k=0:j-1
          oprob = outside[type3rule.leftindex,k+1,j+b+1]
          iprob = inside[type3rule.rightindices[1], k+1,j-1+1]
          tmp = logsumexp(tmp, oprob+iprob)
        end
        outside[type3rule.rightindices[2], j+1, j+b+1] = logsumexp(outside[type3rule.rightindices[2], j+1, j+b+1], tmp + type3rule.logprob)
      end

      if j >= 1 && j+b+1 <len
        for type2rule in type2rules
          outside[type2rule.rightindices[2],j+1,j+b+1] = logsumexp(outside[type2rule.rightindices[2],j+1,j+b+1], type2rule.logprob + pairedlogprobs[j-1+1,j+b+1+1] + outside[type2rule.leftindex,j-1+1,j+b+1+1])
        end
      end
    end
  end
  return outside
end

function computebasepairprobs(inside::Array{Float64,3}, outside::Array{Float64,3}, unpairedlogprobs::Array{Float64,1}, pairedlogprobs::Array{Float64,2}, grammar::KH99)
  type1rules = grammar.type1rules
  type2rules = grammar.type2rules
  type3rules = grammar.type3rules

  numnonterminals = 3
  len = length(unpairedlogprobs)
  res = ones(Float64, len)*-Inf
  Z = inside[1,1,len]
  for j=1:len
    tmp = -Inf
    for type1rule in type1rules
      #inside[ruleindex[type1rule.left], i,i] = unpairedprobs[i]*type1rule.prob
      tmp = logsumexp(tmp, outside[type1rule.leftindex,j,j]+type1rule.logprob)
    end
    res[j] = tmp+unpairedlogprobs[j] - Z
  end

  pairres = ones(Float64, len, len)*-Inf
  for j=0:len-1
    for k=j+1:len-1
      tmp = -Inf
      for type2rule in type2rules
        inres = inside[type2rule.rightindices[2], j+1+1,k-1+1]
        outres = outside[type2rule.leftindex, j+1,k+1]
        tmp = logsumexp(tmp, outres+type2rule.logprob+inres)
      end
      pairres[j+1,k+1] = tmp+pairedlogprobs[j+1,k+1]-Z
      pairres[k+1,j+1] = pairres[j+1,k+1]
    end
  end

  m = max(maximum(res), maximum(pairres))

  res = exp.(res.-m)
  pairres = exp.(pairres.-m)

  for x=1:len
    normalise = res[x]+sum(pairres[x,:])
    res[x] /= normalise
    pairres[x,:] /= normalise
  end
  #=
  normalise = res[1]+sum(pairres[1,:])
  res /= normalise
  pairres /= normalise=#
  return res,pairres
end

function computebasepairprobs2(inside::Array{Float64,3}, outside::Array{Float64,3}, unpairedlogprobs::Array{Float64,1}, pairedlogprobs::Array{Float64,2}, grammar::KH99)
  type1rules = grammar.type1rules
  type2rules = grammar.type2rules
  type3rules = grammar.type3rules

  numnonterminals = 3
  len = length(unpairedlogprobs)
  res = ones(Float64, len)*-Inf
  Z = inside[1,1,len]
  for j=1:len
    tmp = -Inf
    for type1rule in type1rules
      #inside[ruleindex[type1rule.left], i,i] = unpairedprobs[i]*type1rule.prob
      tmp = logsumexp(tmp, outside[type1rule.leftindex,j,j]+type1rule.logprob)
    end
    res[j] = tmp+unpairedlogprobs[j]
  end

  pairres = ones(Float64, len, len)*-Inf
  for j=0:len-1
    for k=j+1:len-1
      tmp = -Inf
      for type2rule in type2rules
        inres = inside[type2rule.rightindices[2], j+1+1,k-1+1]
        outres = outside[type2rule.leftindex, j+1,k+1]
        tmp = logsumexp(tmp, outres+type2rule.logprob+inres)
      end
      pairres[j+1,k+1] = tmp+pairedlogprobs[j+1,k+1]
      pairres[k+1,j+1] = pairres[j+1,k+1]
    end
  end

  return res, pairres
end

function samplestructurehelper(rng::AbstractRNG, inside::Array{Float64,3}, pairedlogprobs::Array{Float64,2}, unpairedlogprobs::Array{Float64,1}, parentsymbol::Char, x::Int, y::Int, paired::Array{Int,1}, grammar::KH99, B::Float64)
  type1rules = grammar.type1rules
  type2rules = grammar.type2rules
  type3rules = grammar.type3rules
  stack = Tuple{Char,Int,Int}[]
  push!(stack, (parentsymbol, x, y))
  while length(stack) > 0
    parentsymbol, x, y = pop!(stack)
    rulearr = Rule[]
    if x == y  # unpaired
    elseif x < y
      j = x
      b = y-j+1

      logliks = Float64[]
      rule = Tuple[]

      for type3rule in type3rules  # bifurcations / branching
        if type3rule.left == parentsymbol
          for h=j:j+b-2
            prob1 = inside[type3rule.rightindices[1],j,h]
            prob2 = inside[type3rule.rightindices[2],h+1,j+b-1]
            push!(logliks, type3rule.logprob + prob1 + prob2)
            push!(rule, (3,j,h,h+1,j+b-1,type3rule.right[1],type3rule.right[2]))
            push!(rulearr, type3rule)
          end
        end
      end

      for type2rule in type2rules # base-pairings
        if type2rule.left == parentsymbol
          push!(logliks, (type2rule.logprob + pairedlogprobs[j,j+b-1]) + inside[type2rule.rightindices[2],j+1,j+b-2])
          push!(rule, (2,j,j+b-1,type2rule.right[2]))
          push!(rulearr, type2rule)
        end
      end

      logliks *= B
      s = -Inf
      for v in logliks
        s = logsumexp(s,v)
      end
      r = CommonUtils.sample(rng, exp.(logliks.-s))  # sample rule used to generate bifurcation or base-pairing

      if rule[r][1] == 2 # base-pairing, recursively sample the sub-alignment
        a = rule[r][2]
        b = rule[r][3]
        paired[a] = b # store pairing
        paired[b] = a # store pairing
        push!(stack, (rule[r][4], a+1, b-1))
      elseif rule[r][1] == 3  # bifurcation, recursively sample the two sub-alignments
        push!(stack, (rule[r][6], rule[r][2], rule[r][3]))
        push!(stack, (rule[r][7], rule[r][4], rule[r][5]))
      end
    end
  end

  return calculateKH99prior(paired)
end

function samplestructure(rng::AbstractRNG, inside::Array{Float64,3}, pairedlogprobs::Array{Float64,2}, unpairedlogprobs::Array{Float64,1}, x::Int, y::Int, paired::Array{Int,1}, grammar::KH99, B::Float64)
  fill!(paired,0)
  priorll = samplestructurehelper(rng, inside, pairedlogprobs, unpairedlogprobs, 'S', x, y, paired, grammar, B)
  return priorll
end

mutable struct Block
  blockid::Int
  gridid::Int
  blockxstart::Int
  blockystart::Int

  function Block(blockid::Int,gridid::Int, blockxstart::Int, blockystart::Int)
    new(blockid,gridid, blockxstart,blockystart)
  end
end

function computeblock(block::Block, blen::Int, inside::Array{Float64,3}, unpairedlogprobs::Array{Float64,1}, pairedlogprobs::Array{Float64,2}, type2rules::Array{Rule,1}, type3rules::Array{Rule,1}, B::Float64=1.0)
#function computeblock(block::Block, blen::Int, inside::SharedArray{Float64,3}, unpairedlogprobs::SharedArray{Float64,1}, pairedlogprobs::SharedArray{Float64,2}, type2rules::Array{Rule,1}, type3rules::Array{Rule,1}, B::Float64=1.0)
  len = length(unpairedlogprobs)
  ##tmpvec = zeros(Float64,len)
  for q=1:blen
    for r=1:blen
      x = block.blockxstart + (blen - q + 1)
      y = block.blockystart + r
      if x != y &&  x <= len && y <= len && x < y
        j = x
        b = y - j + 1

        for type3rule in type3rules
          tmp = -Inf
          for h=j:j+b-2
            prob1 = inside[type3rule.rightindices[1],j,h]
            prob2 = inside[type3rule.rightindices[2],h+1,y]
            #tmpvec[h] = prob1+prob2
            tmp = logsumexp(tmp, prob1+prob2)
          end
          #tmp = logsumexp(tmpvec,j,y-1)
          inside[type3rule.leftindex,j,y] = logsumexp(inside[type3rule.leftindex,j,y], type3rule.logprob*B+tmp)
        end

        for type2rule in type2rules
          inside[type2rule.leftindex,j,y] = logsumexp(inside[type2rule.leftindex,j,y], (type2rule.logprob + pairedlogprobs[j,y])*B + inside[type2rule.rightindices[2],j+1,y-1])
        end
      end
    end
  end

  return inside
end

function computeinsideparallel(unpairedlogprobsin::Array{Float64,1}, pairedlogprobsin::Array{Float64,2}, grammar::KH99, B::Float64=1.0)
  type1rules = grammar.type1rules
  type2rules = grammar.type2rules
  type3rules = grammar.type3rules
  numnonterminals = 3
  len = length(unpairedlogprobsin)
  #inside = SharedArray(Float64, numnonterminals, len, len)
  inside = zeros(Float64, numnonterminals, len, len)
  for symbol=1:size(inside,1)
    for i=1:len
      for j=1:len
        inside[symbol,i,j] = -Inf
      end
    end
  end

  #unpairedlogprobs = convert(SharedArray,unpairedlogprobsin)
  #pairedlogprobs = convert(SharedArray,pairedlogprobsin)
  unpairedlogprobs = unpairedlogprobsin
  pairedlogprobs = pairedlogprobsin


  for type1rule in type1rules
    for i=1:len
      inside[type1rule.leftindex, i,i] = (unpairedlogprobs[i] + type1rule.logprob)*B
    end
  end

  blen = max(50, div(len,32))
  startb = 1
  startj = 1
  blocks = Block[]
  blockx = 0
  blocky = 0
  blocky2 = 0
  blockid = 1
  gridid = 1
  grids = []
  while blocky2 <= len
    grid = Block[]
    while blocky <= len
      block = Block(0,0, blockx,blocky)
      push!(grid, block)
      blockx += blen
      blocky += blen
    end
    for block in grid
      block.gridid = gridid
      block.blockid = blockid
      blockid += 1
    end
    push!(grids, grid)
    gridid += 1
    blocky2 += blen
    blockx = 0
    blocky = blocky2
  end

  for grid in grids
    #@sync @parallel for block in grid
    for block in grid
      computeblock(block, blen, inside, unpairedlogprobs, pairedlogprobs, type2rules, type3rules, B)
    end
  end

  return convert(Array,inside)
end

#include("InsideOutsideAlgorithm.jl")

type Block
  blockid::Int
  gridid::Int
  blockxstart::Int
  blockystart::Int

  function Block(blockid::Int,gridid::Int, blockxstart::Int, blockystart::Int)
    new(blockid,gridid, blockxstart,blockystart)
  end
end

function computeblock(block::Block, blen::Int, inside::SharedArray{Float64,3}, unpairedlogprobs::SharedArray{Float64,1}, pairedlogprobs::SharedArray{Float64,2}, type2rules::Array{Rule,1}, type3rules::Array{Rule,1}, B::Float64=1.0)
  len = length(unpairedlogprobs)
  tmpvec = zeros(Float64,len)
  for q=1:blen
    for r=1:blen
      x = block.blockxstart + (blen - q + 1)
      y = block.blockystart + r
      if x != y &&  x <= len && y <= len && x < y
        #matrix[x,y] = c
        #c += 1
        j = x
        b = y - j + 1

        for type3rule in type3rules
          for h=j:j+b-2
            prob1 = inside[type3rule.rightindices[1],j,h]
            prob2 = inside[type3rule.rightindices[2],h+1,y]
            tmpvec[h] = prob1+prob2
          end
          tmp = logsumexp(tmpvec,j,y-1)
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
  inside = SharedArray(Float64, numnonterminals, len, len)
  for symbol=1:size(inside,1)
    for i=1:len
      for j=1:len
        inside[symbol,i,j] = -Inf
      end
    end
  end

  unpairedlogprobs = convert(SharedArray,unpairedlogprobsin)
  pairedlogprobs = convert(SharedArray,pairedlogprobsin)


  for type1rule in type1rules
    for i=1:len
      inside[type1rule.leftindex, i,i] = (unpairedlogprobs[i] + type1rule.logprob)*B
    end
  end

  blen = div(len,16)
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


  #=
  matrix = zeros(Int, len, len)

  c = 1
  for i=1:len
    matrix[i,i] = -1
    c += 1
  end==#

  for grid in grids
    @sync @parallel for block in grid
      computeblock(block, blen, inside, unpairedlogprobs, pairedlogprobs, type2rules, type3rules, B)
    end
  end

  #println(inside[1,:,:])
  #=
  v = inside[1,:,:]
  for i=1:size(v,1)
    ret = ""
    for j=1:size(v,2)
      ret = string(ret,v[i,j],"\t")
    end
    println(ret)
  end=#

  return convert(Array,inside)
end

function computeblock2(block::Block, blen::Int, inside::Array{Float64,3}, unpairedlogprobs::Array{Float64,1}, pairedlogprobs::Array{Float64,2}, type2rules::Array{Rule,1}, type3rules::Array{Rule,1}, B::Float64=1.0)
  len = length(unpairedlogprobs)
  tmpvec = zeros(Float64,len)
  for q=1:blen
    for r=1:blen
      x = block.blockxstart + (blen - q + 1)
      y = block.blockystart + r
      if x != y &&  x <= len && y <= len && x < y
        #matrix[x,y] = c
        #c += 1
        j = x
        b = y - j + 1

        for type3rule in type3rules
          for h=j:j+b-2
            prob1 = inside[type3rule.rightindices[1],j,h]
            prob2 = inside[type3rule.rightindices[2],h+1,y]
            tmpvec[h] = prob1+prob2
          end
          tmp = logsumexp(tmpvec,j,y-1)
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

function computeinsideparallel2(unpairedlogprobs::Array{Float64,1}, pairedlogprobs::Array{Float64,2}, grammar::KH99, B::Float64=1.0)
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

  blen = div(len,16)
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


  #=
  matrix = zeros(Int, len, len)

  c = 1
  for i=1:len
    matrix[i,i] = -1
    c += 1
  end==#

  for grid in grids
    temps = Array{Float64,3}[]
    refs = []
    for block in grid
      println(block)
      ref = @spawn computeblock(block, blen, copy(inside), unpairedlogprobs, pairedlogprobs, type2rules, type3rules, B)
      push!(refs,ref)
    end

    index = 1
    for ref in refs
      block = grid[index]
      tempinside = fetch(ref)
      for symbol=1:size(inside,1)
        for q=1:blen
          for r=1:blen
            x = block.blockxstart + (blen - q + 1)
            y = block.blockystart + r
            if x != y &&  x <= len && y <= len && x < y
                j = x
                b = y - j + 1
                inside[symbol,j,y] = tempinside[symbol,j,y]
            end
          end
        end
      end
      index += 1
    end

  end

  #println(inside[1,:,:])
  #=
  v = inside[1,:,:]
  for i=1:size(v,1)
    ret = ""
    for j=1:size(v,2)
      ret = string(ret,v[i,j],"\t")
    end
    println(ret)
  end=#

  return inside
end
#=
len = 14
unpairedlogprobs = zeros(Float64,len)
pairedlogprobs = zeros(Float64,len,len)

inside = computeinsideparallel(unpairedlogprobs, pairedlogprobs, KH99(), 1.0)
println(inside[1,:,:])=#

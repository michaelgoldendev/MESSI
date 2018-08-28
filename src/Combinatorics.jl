function log2approxcount(n::Int)
  return 0.142 - 1.5*log2(n) + 1.388*n
end

function countstructures(n::Int, mindist::Int=2, cache::Dict{Int,BigInt}=Dict{Int,BigInt}())
  for z=0:mindist
    cache[z] = BigInt(1)
  end
  for j=1:n
    count = S(j,mindist,cache)
    println(j, "\t", count, "\t", 2.0^log2approxcount(j))
  end
  return cache
end


function S(n::Int, mindist::Int=2, cache::Dict{Int,BigInt}=Dict{Int,BigInt}())
  if haskey(cache,n)
    return cache[n]
  else
    v = BigInt(1)
    if n >= mindist + 1
      v = S(n-1,mindist,cache)
      for k=1:n-(mindist+1)
          v += S(k-1,mindist, cache)*S(n-k-1,mindist, cache)
      end
    end
    cache[n] = v
    return v
  end
end



function logsumexp(a::Float64, b::Float64)
    d = a-b
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

function logcountstructures(n::Int, mindist::Int=2)
  cache = Dict{Int,Float64}()
  for j=1:n
    count = logS(j,mindist,cache)
    exact = count/log(2.0)
    approx = log2approxcount(j)
    println("logapprox\t", j, "\t", exact, "\t", approx, "\t", exact/approx, "\t", count/log(10.0))
  end
  return cache
end

function logS(n::Int, mindist::Int=2, cache::Dict{Int,Float64}=Dict{Int,Float64}())
  if haskey(cache,n)
    return cache[n]
  else
    v = 0.0
    if n >= mindist + 1
      v = logS(n-1,mindist,cache)
      for k=1:n-(mindist+1)
          v = logsumexp(v, logS(k-1,mindist, cache) + logS(n-k-1,mindist, cache))
      end
    end
    cache[n] = v
    return v
  end
end

#=
function S(i::Int, j::Int, mindist::Int=2, cache::Dict{Tuple{Int,Int},BigInt}=Dict{Tuple{Int,Int},BigInt}())
  key = (i,j)
  if haskey(cache,key)
    return cache[key]
  elseif abs(i-j) > mindist+1
    v = S(i+1,j,mindist,cache)

  end
end=#
using Formatting
upton = 200

logcountstructures(200, 3)

cache0 = Dict{Int,BigInt}()
countstructures(upton, 0, cache0)

cache1 = Dict{Int,BigInt}()
countstructures(upton, 1, cache1)

cache2  = Dict{Int,BigInt}()
countstructures(upton, 2, cache2)

cache3  = Dict{Int,BigInt}()
countstructures(upton, 3, cache3)
for i=1:20
  #println(i," & ", sprintf1("%.1e", BigFloat(cache3[i])), "\\tabularnewline")
  println(i," & ", cache3[i], "\\tabularnewline")
end
exit()

ret = ""
for n=1:upton
  if n % 2 == 1
    #ret = string(ret, "\\rowcolor{black!20} ")
  end
  #ret = string(ret, n, " & ", cache0[n], " & ", cache1[n], " & ", cache2[n], " & ", cache3[n], "\\tabularnewline\n")
  ret = string(ret, n, "\t", cache0[n], "\t", cache1[n], "\t", cache2[n], "\t", cache3[n], "\n")
end
println(ret)

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

using Printf
push!(LOAD_PATH,@__DIR__)
using Mapping

function enumeratesubstructures(paired::Array{Int,1}, minlength::Int, maxlength::Int)
    structures = []
    #lastStructureAdded = false
    i = 0
    while i+1 < length(paired)
        i += 1
        if paired[i] > i
            slen = paired[i] - i + 1
            if slen <= maxlength
                push!(structures, (i,paired[i]))
                i = paired[i]
            end
        end
    end
    return structures
end


function generatesubstructures(paired::Array{Int,1}, minlength::Int, maxsublength::Int, maxlength::Int)
    endpoints = Tuple{Int,Int}[]
    len = length(paired)
    i = 1
    startindex = 1
    endindex = 1
    while i <= len
        if paired[i] > i
            sublen = paired[i]-i+1
            if sublen <= maxsublength
                startindex = i
                endindex = paired[i]
                push!(endpoints, (startindex,endindex))
                i = endindex
            end
        end
        i += 1
    end
    endindex = endpoints[end][end]+1
    println("Z",paired[endindex:end], "\t", endpoints)

    finalendpoints = Tuple{Int,Int}[]
    startr = 1
    endr = 1
    r = 1
    while r <= length(endpoints)
        sublen = endpoints[r][2] - endpoints[startr][1] + 1
        println(endpoints)
        println(r,"\t", startr, "\t", endr, "\t", sublen, "\t", maxlength)
        if sublen <= maxlength

        else
            println("E:", endpoints[startr][1])
            if r-1 > 0

                push!(finalendpoints, (endpoints[startr][1], endpoints[r-1][2]))
            end
            startr = r
        end
        r += 1
    end
    if startr <= length(endpoints)
        push!(finalendpoints, (endpoints[startr][1], endpoints[end][2]))
    end
    println("Y", finalendpoints, "\t", startr, "\t", length(endpoints))

    substructures = Array{Int,1}[]
    starta = 1
    enda = 1
    finalendpoints2 = Tuple{Int,Int}[]
    for (startindex,endindex) in finalendpoints
        enda = startindex-1
        zlen = enda-starta+1
        if zlen > 0
            push!(finalendpoints2, (starta, enda))
            push!(substructures, zeros(Int,zlen))
            #println(starta,"\t",enda,"\t", getdotbracketstring(substructures[end]))
        end
        substructure = copy(paired[startindex:endindex])
        for j=1:length(substructure)
            if substructure[j] != 0
                substructure[j] = substructure[j] - startindex + 1
            elseif !(startindex <= substructure[j] <= endindex)
                substructure[j]  = 0
            end
        end
        starta = endindex+1

        push!(finalendpoints2, (startindex, endindex))
        #println(startindex,"\t",endindex,"\t", getdotbracketstring(substructure))
        push!(substructures, substructure)
    end
    starta = finalendpoints[end][2]+1
    enda = length(paired)
    zlen = enda-starta+1
    if enda-starta > 0
        push!(finalendpoints2, (starta, enda))
        push!(substructures, zeros(Int,zlen))
        #println(starta,"\t",enda,"\t", getdotbracketstring(substructures[end]))
    end
    println("A", length(finalendpoints2), "\t", length(substructures))

    for substructure in substructures
        for i=1:length(substructure)
            substructure[i] = max(0, substructure[i])
        end
    end

    return finalendpoints2, substructures
end

using HypothesisTests
function mannwhitneyu(x,y)
    nx = length(x)
    ny = length(y)
    bigN = nx+ny
    n = ((nx * ny) / (bigN * (bigN - 1)))
    tieCorrectionFactor = 0.0
    d = ((bigN^3.0 - bigN) / 12.0 - tieCorrectionFactor)
    variance = sqrt(n*d)
    if length(x) > 2 && length(y) > 2
        mann = MannWhitneyUTest(x, y)
        zscore =  (mann.U - (nx * ny / 2.0)) / variance
        return zscore, pvalue(mann)
    else
        return 0.0, 1.0
    end
end

function substructureranking(paired, data, references,mapping,revmapping, ascending=true, bothnucleotides::Bool=false, unpairedsites::Bool=false)
    for ref in references
        refpos = revmapping[ref[1]]
        println(ref[2])
        println("mapping: ", ref[1], " -> ", refpos)
    end

    substructures = enumeratesubstructures(paired,20,250)
    ls = []
    finalls = []
    id = 1
    for substructure in substructures
        startx = substructure[1]
        endx = substructure[2]
        subpaired = paired[startx:endx]
        outer = Float64[]
        inner = Float64[]
        for i=1:length(paired)
            x = i
            if x > 0 && (paired[i] > i || (bothnucleotides && paired[i] > 0) || unpairedsites) && !isnan(data[x])
                if (startx <= i <= endx) && (startx <= paired[i] <= endx)
                    push!(inner,data[x])
                else
                    push!(outer,data[x])
                end
            end
        end

        if length(inner) > 0
            zscore,pvalue = mannwhitneyu(inner,outer)
            slen = endx-startx+1
            reference = ""
            repeatcount = 0
            #=
            for ref in references
            refpos = revmapping[ref[1]]
            if startx <= refpos <= endx
            reference = ref[2]
            repeatcount += 1
        end
    end=#

    for ref in references
        refstart = ref[1][1]
        refend = ref[1][2]
        for r=refstart:refend
            refpos = revmapping[r]
            if startx <= refpos <= endx
                reference = ref[2]
            end
        end
    end

    if repeatcount > 1
        println("REGION REPEATS ",repeatcount, "\t", reference)
    end
    med = "-"
    if !isnan(median(inner))
        med = @sprintf("%.2f", median(inner))
    end
    mappedstring = string(get(mapping, startx, 0), " - ", get(mapping, endx,0))
    if get(mapping, startx, 0) == 0 || get(mapping, endx,0) == 0
        mappedstring = "NA"
    end
    ret = string(startx, " - ", endx, " & ",  mappedstring, " & ", slen, " & ", reference, " & ", med, " & ", @sprintf("%.2f", zscore), "\\tabularnewline\n")
    sign = 1.0
    if !ascending
        sign = -1.0
    end
    push!(ls,(sign*zscore,ret,id))
    push!(finalls, (sign*zscore, median(inner), startx, endx))
    id += 1
end
end
#ls = reverse(sort(ls))
finalls = sort(finalls)
ls = sort(ls)
texls = []
csv = ""
for i=1:length(ls)
    csv = string(csv, i, ",", replace(replace(ls[i][2], " & " => ","), "\\tabularnewline" => ""))
    if i % 2 == 1
        push!(texls, string("\t\t\\rowcolor{black!20} ",i," & ", ls[i][2]))
    else
        push!(texls, string("\t\t", i," & ", ls[i][2]))
    end
end
tex = join(texls,"")
println(tex)
rankings = Float64[Float64(ls[i][3]) for i=1:length(ls)]
return rankings, finalls, tex, csv
end

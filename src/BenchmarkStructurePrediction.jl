include("RNATools.jl")

push!(LOAD_PATH,joinpath(@__DIR__))
#include("/media/michael/Sandisk500GB/juliamolev/src/MolecularEvolution.jl")
using MolecularEvolution
using Formatting

v = []

resultsdir = "results/"
rfamdir = "../datasets/RFAM/3D/"

for f in filter(x -> match(r".+\.stockholm\.txt", x) != nothing, readdir(rfamdir))
  rfamid = match(r"(.+)\.stockholm\.txt", f)[1]
  alignmentfile = joinpath(resultsdir, rfamid, string(rfamid,".fas.norm"))
  treefile = joinpath(resultsdir, rfamid, string(rfamid,".nwk"))
  predctfile = joinpath(resultsdir, rfamid, string(rfamid,".maxstructuretrunc"))
  truectfile = joinpath(rfamdir, string(rfamid,".ct"))

  if isfile(alignmentfile) && isfile(treefile) && isfile(truectfile) # && isfile(predctfile)
    println(rfamid,"\t",isfile(alignmentfile),"\t",isfile(treefile),"\t",isfile(predctfile),"\t",isfile(truectfile))
    push!(v, (alignmentfile,treefile,predctfile,truectfile))
  end
end
#exit()
#push!(v, ("/media/michael/Sandisk500GB/data/RF01854/RF01854.fas.norm", "/media/michael/Sandisk500GB/data/RF01854/RF01854.nwk", "/media/michael/Sandisk500GB/data/RF01854/RF01854.maxstructuretrunc", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF01854.pfam.txt.ct"))


function getbenchmarkvalue(method::Int, metric::Int, ourmethod, ppfoldmethod, rnaalifoldmethod)
  str = ""
  maxval, maxindex = findmax([ourmethod[metric][end],ppfoldmethod[metric][end],rnaalifoldmethod[metric][end]])
  if method == 1
    str = Formatting.fmt(".3f",ourmethod[metric][end])
  elseif method == 2
    str = Formatting.fmt(".3f",ppfoldmethod[metric][end])
  elseif method == 3
    str = Formatting.fmt(".3f",rnaalifoldmethod[metric][end])
  end
  if maxindex == method
    str = string("\\textbf{",str, "}")
  end
  return str
end



dir = "temp/"
if !isdir(dir)
  mkdir(dir)
end

nmetrics = 6
ourmethod = [Float64[] for m=1:nmetrics]
ppfoldmethod = [Float64[] for m=1:nmetrics]
rnaalifoldmethod = [Float64[] for m=1:nmetrics]

ret = ""
retmax = ""
for x in v
  global ret
  global retmax
  alignmentfile = x[1]
  m = match(r"^(.*/)([^/]*)$", alignmentfile)
  filename = m[2]
  println(filename)
  treefile = x[2]
  infile = x[3]
  structurefile = x[4]
  dataset = Dataset(getalignmentfromfastaandnewick(alignmentfile, treefile))

  ismaxstructure = false
  consensus = nothing
  if endswith(infile, ".structures")
    structures = readdbnstructures(infile)
    maskgapped!(structures, dataset.gapfrequency,0.5)
    #=for s=1:length(structures)
      structures[s] = removelonelybasepairs(structures[s])
    end=#
    singleprobs, pairedprobs = getpairprobs(structures[max(1,div(length(structures),3)):end])
    #pairedprobs = maskgapped!(pairedprobs,dataset.gapfrequency,0.5)
    consensus = getPosteriorDecodingConsensusStructure(pairedprobs, singleprobs)
    consensus = removelonelybasepairs(consensus)
  elseif endswith(infile,".max.consensus.ct") || endswith(infile, ".maxstructure") ||  endswith(infile, ".maxstructuretrunc")
    if isfile(infile)
      consensusseq, consensus = readctfile(infile)
      ismaxstructure = true
    end
  end

  if endswith(structurefile, ".dbn")
    real = readdbnfile(structurefile)
  elseif endswith(structurefile, ".ct")
    realseq, real = readctfile(structurefile)
  end

  ppfoldfile = string(alignmentfile, ".ppfold.ct")
  ppfoldseq = nothing
  ppfold = nothing
  if !isfile(ppfoldfile)
    ppfoldseq,ppfold = runppfold(alignmentfile, treefile, dir)
    writectfile(ppfold, ppfoldseq, ppfoldfile)
  else
    ppfoldseq,ppfold = readctfile(ppfoldfile)
  end

  rnaalifoldfile = string(alignmentfile, ".rnaalifoldfile.ct")
  rnalifoldseq = nothing
  rnaalifold = nothing
  if !isfile(rnaalifoldfile)
    rnaalifold = runrnaalifold(alignmentfile, treefile, dir)
    writectfile(rnaalifold, ppfoldseq, rnaalifoldfile)
  else
    rnalifoldseq,rnaalifold = readctfile(rnaalifoldfile)
  end

  if consensus != nothing
    push!(ourmethod[1], precision(real, consensus))
    push!(ourmethod[2], recall(real, consensus))
    push!(ourmethod[3], f1score(real, consensus))
    push!(ourmethod[5], 1.0 - normalisedweightedmountaindistance(real, consensus))
    push!(ourmethod[6], 1.0 - normalisedmountaindistance(real,consensus,1.0))
    push!(ourmethod[4], mcc(real,consensus))
    push!(ppfoldmethod[1], precision(real, ppfold))
    push!(ppfoldmethod[2], recall(real, ppfold))
    push!(ppfoldmethod[3], f1score(real, ppfold))
    push!(ppfoldmethod[5], 1.0 - normalisedweightedmountaindistance(real, ppfold))
    push!(ppfoldmethod[6], 1.0 - normalisedmountaindistance(real,ppfold,1.0))
    push!(ppfoldmethod[4], mcc(real,ppfold))
    push!(rnaalifoldmethod[1], precision(real, rnaalifold))
    push!(rnaalifoldmethod[2], recall(real, rnaalifold))
    push!(rnaalifoldmethod[3], f1score(real, rnaalifold))
    push!(rnaalifoldmethod[5], 1.0 - normalisedweightedmountaindistance(real, rnaalifold))
    push!(rnaalifoldmethod[6], 1.0 - normalisedmountaindistance(real,rnaalifold,1.0))
    push!(rnaalifoldmethod[4], mcc(real,rnaalifold))


    println(ourmethod)
    println(ppfoldmethod)
    println(rnaalifoldmethod)

    println(alignmentfile)
    println("FARCE\t", Formatting.fmt(".3f",precision(real, consensus)), "\t", Formatting.fmt(".3f",recall(real, consensus)), "\t", Formatting.fmt(".3f",f1score(real, consensus)), "\t", Formatting.fmt(".3f",1.0 - normalisedweightedmountaindistance(real,consensus)), "\t", Formatting.fmt(".3f",1.0 - normalisedmountaindistance(real,consensus,1.0)))
    println("PPfold\t", Formatting.fmt(".3f",precision(real, ppfold)), "\t", Formatting.fmt(".3f",recall(real, ppfold)), "\t", Formatting.fmt(".3f",f1score(real, ppfold)), "\t", Formatting.fmt(".3f",1.0 - normalisedweightedmountaindistance(real,ppfold)), "\t", Formatting.fmt(".3f",1.0 - normalisedmountaindistance(real,ppfold,1.0)))
    println("RNAalifold\t", Formatting.fmt(".3f",precision(real, rnaalifold)), "\t", Formatting.fmt(".3f",recall(real, rnaalifold)), "\t", Formatting.fmt(".3f",f1score(real, rnaalifold)), "\t", Formatting.fmt(".3f",1.0 - normalisedweightedmountaindistance(real,rnaalifold)), "\t", Formatting.fmt(".3f",1.0 - normalisedmountaindistance(real,rnaalifold,1.0)))



    str = ""
    str = string(str, " & ", "Our method & ", getbenchmarkvalue(1, 1, ourmethod, ppfoldmethod, rnaalifoldmethod), " & ", getbenchmarkvalue(1, 2, ourmethod, ppfoldmethod, rnaalifoldmethod), " & ", getbenchmarkvalue(1, 3, ourmethod, ppfoldmethod, rnaalifoldmethod), " & ", getbenchmarkvalue(1, 4, ourmethod, ppfoldmethod, rnaalifoldmethod), " & ", getbenchmarkvalue(1, 5, ourmethod, ppfoldmethod, rnaalifoldmethod), " & ", getbenchmarkvalue(1, 6, ourmethod, ppfoldmethod, rnaalifoldmethod), "\\\\*\n")
    str = string(str, filename, " & ", "PPfold & ", getbenchmarkvalue(2, 1, ourmethod, ppfoldmethod, rnaalifoldmethod), " & ", getbenchmarkvalue(2, 2, ourmethod, ppfoldmethod, rnaalifoldmethod), " & ", getbenchmarkvalue(2, 3, ourmethod, ppfoldmethod, rnaalifoldmethod), " & ", getbenchmarkvalue(2, 4, ourmethod, ppfoldmethod, rnaalifoldmethod), " & ", getbenchmarkvalue(2, 5, ourmethod, ppfoldmethod, rnaalifoldmethod), " & ", getbenchmarkvalue(2, 6, ourmethod, ppfoldmethod, rnaalifoldmethod), " \\\\*\n")
    str = string(str, " & ", "RNAalifold & ", getbenchmarkvalue(3, 1, ourmethod, ppfoldmethod, rnaalifoldmethod), " & ",getbenchmarkvalue(3, 2, ourmethod, ppfoldmethod, rnaalifoldmethod), " & ", getbenchmarkvalue(3, 3, ourmethod, ppfoldmethod, rnaalifoldmethod), " & ", getbenchmarkvalue(3, 4, ourmethod, ppfoldmethod, rnaalifoldmethod), " & ", getbenchmarkvalue(3, 5, ourmethod, ppfoldmethod, rnaalifoldmethod), " & ", getbenchmarkvalue(3, 6, ourmethod, ppfoldmethod, rnaalifoldmethod), " \\tabularnewline\n")
    str = string(str, "\\midrule\n")

   
    if ismaxstructure
      retmax = string(retmax, str)
    else
      ret = string(ret, str)
    end
  end
end

println(ret)
println("----------------------------------------------------------")
println(retmax)

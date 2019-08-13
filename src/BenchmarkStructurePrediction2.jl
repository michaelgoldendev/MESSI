include("RNATools.jl")

include("/media/michael/Sandisk500GB/juliamolev/src/MolecularEvolution.jl")
using MolecularEvolution
using Formatting

v = []

push!(v, ("/media/michael/Sandisk500GB/data/RF01854/RF01854.fas.norm", "/media/michael/Sandisk500GB/data/RF01854/RF01854.nwk", "/media/michael/Sandisk500GB/data/RF01854/RF01854.maxstructuretrunc", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF01854.pfam.txt.ct"))
#push!(v, ("/media/michael/Sandisk500GB/data/RF01854/RF01854.fas.norm", "/media/michael/Sandisk500GB/data/RF01854/RF01854.nwk", "/media/michael/Sandisk500GB/data/RF01854/RF01854.maxstructure", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF01854.pfam.txt.ct"))

push!(v, ("/media/michael/Sandisk500GB/data/RF01846/RF01846.fas.norm", "/media/michael/Sandisk500GB/data/RF01846/RF01846.nwk", "/media/michael/Sandisk500GB/data/RF01846/RF01846.maxstructuretrunc", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF01846.pfam.txt.ct"))
#push!(v, ("/media/michael/Sandisk500GB/data/RF01846/RF01846.fas.norm", "/media/michael/Sandisk500GB/data/RF01846/RF01846.nwk", "/media/michael/Sandisk500GB/data/RF01846/RF01846.maxstructure", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF01846.pfam.txt.ct"))

push!(v, ("/media/michael/Sandisk500GB/data/RF00002/RF00002.fas.norm", "/media/michael/Sandisk500GB/data/RF00002/RF00002.nwk", "/media/michael/Sandisk500GB/data/RF00002/RF00002.maxstructuretrunc", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00002.pfam.txt.ct"))
#push!(v, ("/media/michael/Sandisk500GB/data/RF00002/RF00002.fas.norm", "/media/michael/Sandisk500GB/data/RF00002/RF00002.nwk", "/media/michael/Sandisk500GB/data/RF00002/RF00002.maxstructure", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00002.pfam.txt.ct"))

push!(v, ("/media/michael/Sandisk500GB/data/RF02001/RF02001.fas.norm", "/media/michael/Sandisk500GB/data/RF02001/RF02001.nwk", "/media/michael/Sandisk500GB/data/RF02001/RF02001.maxstructuretrunc", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF02001.pfam.txt.ct"))
#push!(v, ("/media/michael/Sandisk500GB/data/RF02001/RF02001.fas.norm", "/media/michael/Sandisk500GB/data/RF02001/RF02001.nwk", "/media/michael/Sandisk500GB/data/RF02001/RF02001.maxstructure", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF02001.pfam.txt.ct"))


push!(v, ("/media/michael/Sandisk500GB/data/RF00379/RF00379.fas.norm", "/media/michael/Sandisk500GB/data/RF00379/RF00379.nwk", "/media/michael/Sandisk500GB/data/RF00379/RF00379.maxstructuretrunc", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00379.pfam.txt.ct"))
#push!(v, ("/media/michael/Sandisk500GB/data/RF00379/RF00379.fas.norm", "/media/michael/Sandisk500GB/data/RF00379/RF00379.nwk", "/media/michael/Sandisk500GB/data/RF00379/RF00379.maxstructure", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00379.pfam.txt.ct"))

push!(v, ("/media/michael/Sandisk500GB/data/RF00380/RF00380.fas.norm", "/media/michael/Sandisk500GB/data/RF00380/RF00380.nwk", "/media/michael/Sandisk500GB/data/RF00380/RF00380.maxstructuretrunc", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00380.pfam.txt.ct"))
#push!(v, ("/media/michael/Sandisk500GB/data/RF00380/RF00380.fas.norm", "/media/michael/Sandisk500GB/data/RF00380/RF00380.nwk", "/media/michael/Sandisk500GB/data/RF00380/RF00380.maxstructure", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00380.pfam.txt.ct"))


push!(v, ("/media/michael/Sandisk500GB/data/RF00012/RF00012.fas.norm", "/media/michael/Sandisk500GB/data/RF00012/RF00012.nwk", "/media/michael/Sandisk500GB/data/RF00012/RF00012.maxstructuretrunc", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00012.pfam.txt.ct"))
#push!(v, ("/media/michael/Sandisk500GB/data/RF00012/RF00012.fas.norm", "/media/michael/Sandisk500GB/data/RF00012/RF00012.nwk", "/media/michael/Sandisk500GB/data/RF00012/RF00012.maxstructure", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00012.pfam.txt.ct"))

push!(v, ("/media/michael/Sandisk500GB/data/RF00100/RF00100.fas.norm", "/media/michael/Sandisk500GB/data/RF00100/RF00100.nwk", "/media/michael/Sandisk500GB/data/RF00100/RF00100.maxstructuretrunc", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00100.pfam.txt.ct"))
#push!(v, ("/media/michael/Sandisk500GB/data/RF00100/RF00100.fas.norm", "/media/michael/Sandisk500GB/data/RF00100/RF00100.nwk", "/media/michael/Sandisk500GB/data/RF00100/RF00100.maxstructure", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00100.pfam.txt.ct"))


push!(v, ("/media/michael/Sandisk500GB/data/RF00020/RF00020.fas.norm", "/media/michael/Sandisk500GB/data/RF00020/RF00020.nwk", "/media/michael/Sandisk500GB/data/RF00020/RF00020.maxstructuretrunc", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00020.pfam.txt.ct"))
#push!(v, ("/media/michael/Sandisk500GB/data/RF00020/RF00020.fas.norm", "/media/michael/Sandisk500GB/data/RF00020/RF00020.nwk", "/media/michael/Sandisk500GB/data/RF00020/RF00020.maxstructure", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00020.pfam.txt.ct"))


push!(v, ("/media/michael/Sandisk500GB/data/RF00174/RF00174.fas.norm", "/media/michael/Sandisk500GB/data/RF00174/RF00174.nwk", "/media/michael/Sandisk500GB/data/RF00174/RF00174.maxstructuretrunc", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00174.pfam.txt.ct"))
#push!(v, ("/media/michael/Sandisk500GB/data/RF00174/RF00174.fas.norm", "/media/michael/Sandisk500GB/data/RF00174/RF00174.nwk", "/media/michael/Sandisk500GB/data/RF00174/RF00174.maxstructure", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00174.pfam.txt.ct"))


push!(v, ("/media/michael/Sandisk500GB/data/RF00001/RF00001.fas.norm", "/media/michael/Sandisk500GB/data/RF00001/RF00001.nwk", "/media/michael/Sandisk500GB/data/RF00001/RF00001.maxstructuretrunc", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00001.pfam.txt.ct"))
#push!(v, ("/media/michael/Sandisk500GB/data/RF00001/RF00001.fas.norm", "/media/michael/Sandisk500GB/data/RF00001/RF00001.nwk", "/media/michael/Sandisk500GB/data/RF00001/RF00001.maxstructure", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00001.pfam.txt.ct"))


push!(v, ("/media/michael/Sandisk500GB/data/RF00004/RF00004.fas.norm", "/media/michael/Sandisk500GB/data/RF00004/RF00004.nwk", "/media/michael/Sandisk500GB/data/RF00004/RF00004.maxstructuretrunc", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00004.pfam.txt.ct"))
#push!(v, ("/media/michael/Sandisk500GB/data/RF00004/RF00004.fas.norm", "/media/michael/Sandisk500GB/data/RF00004/RF00004.nwk", "/media/michael/Sandisk500GB/data/RF00004/RF00004.maxstructure", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00004.pfam.txt.ct"))


push!(v, ("/media/michael/Sandisk500GB/data/RNaseMRP_mcmc/RNaseMRP.fas.norm", "/media/michael/Sandisk500GB/data/RNaseMRP_mcmc/RNaseMRP.nwk", "/media/michael/Sandisk500GB/data/RNaseMRP_mcmc/RNaseMRP.maxstructuretrunc", "/media/michael/Sandisk500GB/dataclean/final/experimental_noncoding/RNase_MRP.dat.ct"))
#push!(v, ("/media/michael/Sandisk500GB/data/RNaseMRP_mcmc/RNaseMRP.fas.norm", "/media/michael/Sandisk500GB/data/RNaseMRP_mcmc/RNaseMRP.nwk", "/media/michael/Sandisk500GB/data/RNaseMRP_mcmc/RNaseMRP.maxstructure", "/media/michael/Sandisk500GB/dataclean/final/experimental_noncoding/RNase_MRP.dat.ct"))

push!(v, ("/media/michael/Sandisk500GB/data/RF02542/RF02542.fas.norm", "/media/michael/Sandisk500GB/data/RF02542/RF02542.nwk", "/media/michael/Sandisk500GB/data/RF02542/RF02542.maxstructuretrunc", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF02542.pfam.txt.ct"))
#push!(v, ("/media/michael/Sandisk500GB/data/RF02542/RF02542.fas.norm", "/media/michael/Sandisk500GB/data/RF02542/RF02542.nwk", "/media/michael/Sandisk500GB/data/RF02542/RF02542.maxstructure", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF02542.pfam.txt.ct"))

push!(v, ("/media/michael/Sandisk500GB/data/RF00026/RF00026.fas.norm", "/media/michael/Sandisk500GB/data/RF00026/RF00026.nwk", "/media/michael/Sandisk500GB/data/RF00026/RF00026.maxstructuretrunc", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00026.pfam.txt.ct"))
#push!(v, ("/media/michael/Sandisk500GB/data/RF00026/RF00026.fas.norm", "/media/michael/Sandisk500GB/data/RF00026/RF00026.nwk", "/media/michael/Sandisk500GB/data/RF00026/RF00026.maxstructure", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00026.pfam.txt.ct"))

push!(v, ("/media/michael/Sandisk500GB/data/RF00017/RF00017.fas.norm", "/media/michael/Sandisk500GB/data/RF00017/RF00017.nwk", "/media/michael/Sandisk500GB/data/RF00017/RF00017.maxstructuretrunc", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00017.pfam.txt.ct"))
#push!(v, ("/media/michael/Sandisk500GB/data/RF00017/RF00017.fas.norm", "/media/michael/Sandisk500GB/data/RF00017/RF00017.nwk", "/media/michael/Sandisk500GB/data/RF00017/RF00017.maxstructure", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00017.pfam.txt.ct"))


push!(v, ("/media/michael/Sandisk500GB/data/RF00011/RF00011.fas.norm", "/media/michael/Sandisk500GB/data/RF00011/RF00011.nwk", "/media/michael/Sandisk500GB/data/RF00011/RF00011.maxstructuretrunc", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00011.pfam.txt.ct"))
#push!(v, ("/media/michael/Sandisk500GB/data/RF00011/RF00011.fas.norm", "/media/michael/Sandisk500GB/data/RF00011/RF00011.nwk", "/media/michael/Sandisk500GB/data/RF00011/RF00011.maxstructure", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00011.pfam.txt.ct"))

push!(v, ("/media/michael/Sandisk500GB/data/RF00010/RF00010.fas.norm", "/media/michael/Sandisk500GB/data/RF00010/RF00010.nwk", "/media/michael/Sandisk500GB/data/RF00010/RF00010.maxstructuretrunc", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00010.pfam.txt.ct"))
#push!(v, ("/media/michael/Sandisk500GB/data/RF00010/RF00010.fas.norm", "/media/michael/Sandisk500GB/data/RF00010/RF00010.nwk", "/media/michael/Sandisk500GB/data/RF00010/RF00010.maxstructure", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00010.pfam.txt.ct"))


push!(v, ("/media/michael/Sandisk500GB/data/ires.fas", "/media/michael/Sandisk500GB/data/ires_max5/ires.nwk", "/media/michael/Sandisk500GB/data/ires_mcmc6/ires.maxstructuretrunc", "/media/michael/Sandisk500GB/data/ires.dbn"))
#push!(v, ("/media/michael/Sandisk500GB/data/ires.fas", "/media/michael/Sandisk500GB/data/ires_max5/ires.nwk", "/media/michael/Sandisk500GB/data/ires_mcmc6/ires.maxstructure", "/media/michael/Sandisk500GB/data/ires.dbn"))

push!(v, ("/media/michael/Sandisk500GB/data/RF02540/RF02540.fas.norm", "/media/michael/Sandisk500GB/data/RF02540/RF02540.nwk", "/media/michael/Sandisk500GB/data/RF02540/RF02540.maxstructuretrunc", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF02540.pfam.txt.ct"))
#push!(v, ("/media/michael/Sandisk500GB/data/RF02540/RF02540.fas.norm", "/media/michael/Sandisk500GB/data/RF02540/RF02540.nwk", "/media/michael/Sandisk500GB/data/RF02540/RF02540.maxstructure", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF02540.pfam.txt.ct"))

push!(v, ("/media/michael/Sandisk500GB/data/RF02541/RF02541.fas.norm", "/media/michael/Sandisk500GB/data/RF02541/RF02541.nwk", "/media/michael/Sandisk500GB/data/RF02541/RF02541.maxstructuretrunc", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF02541.pfam.txt.ct"))
#push!(v, ("/media/michael/Sandisk500GB/data/RF02541/RF02541.fas.norm", "/media/michael/Sandisk500GB/data/RF02541/RF02541.nwk", "/media/michael/Sandisk500GB/data/RF02541/RF02541.maxstructure", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF02541.pfam.txt.ct"))

push!(v, ("/media/michael/Sandisk500GB/data/RF00003/RF00003.fas.norm", "/media/michael/Sandisk500GB/data/RF00003/RF00003.nwk", "/media/michael/Sandisk500GB/data/RF00003/RF00003.maxstructuretrunc", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00003.pfam.txt.ct"))
#push!(v, ("/media/michael/Sandisk500GB/data/RF00003/RF00003.fas.norm", "/media/michael/Sandisk500GB/data/RF00003/RF00003.nwk", "/media/michael/Sandisk500GB/data/RF00003/RF00003.maxstructure", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00003.pfam.txt.ct"))

#=
push!(v, ("/media/michael/Sandisk500GB/data/RF00010/RF00010.fas.norm", "/media/michael/Sandisk500GB/data/RF00010/RF00010.nwk", "/media/michael/Sandisk500GB/data/RF00010/RF00010.B1.0.M1.0.structures", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF00010.pfam.txt.ct"))
#push!(v, ("/media/michael/Sandisk500GB/data/RF02540/RF02540.fas.norm", "/media/michael/Sandisk500GB/data/RF02540/RF02540.nwk", "/media/michael/Sandisk500GB/data/RF02540/RF02540.B1.0.M1.0.structures", "/media/michael/Sandisk500GB/dataclean/final/RFAM/RF02540.pfam.txt.ct"))
push!(v, ("/media/michael/Sandisk500GB/data/ires.fas", "/media/michael/Sandisk500GB/data/ires_max5/ires.nwk", "/media/michael/Sandisk500GB/data/ires_mcmc6/ires.B1.0.M1.0.structures", "/media/michael/Sandisk500GB/data/ires.dbn"))
#push!(v, ("/media/michael/Sandisk500GB/data/ires.fas", "/media/michael/Sandisk500GB/data/ires_max5/ires.nwk", "/media/michael/Sandisk500GB/data/ires_mcmc_tree/ires.B1.0.M1.0.structures", "/media/michael/Sandisk500GB/data/ires.dbn"))
#push!(v, ("/media/michael/Sandisk500GB/data/RNaseMRP_mcmc/RNaseMRP.fas.norm", "/media/michael/Sandisk500GB/data/RNaseMRP_mcmc/RNaseMRP.nwk", "/media/michael/Sandisk500GB/data/RNaseMRP_mcmc_lambda/RNaseMRP.B1.0.M1.0.structures", "/media/michael/Sandisk500GB/dataclean/final/experimental_noncoding/RNase_MRP.dat.ct"))
push!(v, ("/media/michael/Sandisk500GB/data/RNaseMRP_mcmc/RNaseMRP.fas.norm", "/media/michael/Sandisk500GB/data/RNaseMRP_mcmc/RNaseMRP.nwk", "/media/michael/Sandisk500GB/data/RNaseMRP_mcmc/RNaseMRP.B1.0.M1.0.structures", "/media/michael/Sandisk500GB/dataclean/final/experimental_noncoding/RNase_MRP.dat.ct"))
#push!(v, ("/media/michael/Sandisk500GB/data/CsrB.dat.fas", "/media/michael/Sandisk500GB/data/CsrB/CsrB.nwk", "/media/michael/Sandisk500GB/data/CsrB_mcmc/CsrB.B1.0.structures", "/media/michael/Sandisk500GB/data/CsrB.dat.ct"))
push!(v, ("/media/michael/Sandisk500GB/data/U2.dat.fas", "/media/michael/Sandisk500GB/data/U2.dat.fas.nwk", "/media/michael/Sandisk500GB/data/U2_mcmc/U2.B1.0.structures", "/media/michael/Sandisk500GB/data/U2.dbn"))


#push!(v, ("/media/michael/Sandisk500GB/data/ires.fas", "/media/michael/Sandisk500GB/data/ires_thermo/ires.nwk", "/media/michael/Sandisk500GB/data/ires_thermo3/ires.B1.0.structures2", "/media/michael/Sandisk500GB/data/ires.dbn"))
#push!(v, ("/media/michael/Sandisk500GB/data/RNaseP_nuc.dat.fas", "/media/michael/Sandisk500GB/data/RNase/RNase.nwk", "/media/michael/Sandisk500GB/data/RNAseP_mcmc5/RNAseP.B1.0.M1.0.structures", "/media/michael/Sandisk500GB/Dropbox/Oxfold II/Cofold Supplementary/RNaseP_nuc.dat.ct"))
=#


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
  alignmentfile = x[1]
  m = match(r"^(.*/)([^/]*)$", alignmentfile)
  filename = m[2]
  println(filename)
  treefile = x[2]
  infile = x[3]
  structurefile = x[4]
  dataset = Dataset(alignmentfile, treefile)

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
  elseif endswith(infile, ".maxstructure") ||  endswith(infile, ".maxstructuretrunc")
    consensusseq, consensus = readctfile(infile)
    ismaxstructure = true
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

println(ret)
println("----------------------------------------------------------")
println(retmax)

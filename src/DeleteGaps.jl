using FastaIO
#fastafile = "D:\\Dropbox\\dev\\farce-julia\\src\\mcmc\\hivgroupm\\hiv.fas.norm"
fastafile = "/media/michael/Sandisk500GB/Dropbox/dev/farce-julia/src/mcmc/hiv1b/hiv1b.fas.norm"
names = []
sequences = []
seqnametoindex  = Dict{AbstractString,Int}()
FastaIO.FastaReader(fastafile) do fr
 seqindex = 1
 for (desc, seq) in fr
   len = length(seq)
   if len >= 0
     push!(names,desc)
     push!(sequences, seq)
     seqnametoindex[desc] = seqindex
     seqindex += 1
   end
 end
end

nucmapping = Dict('A' => 1, 'C' => 2, 'G' => 3, 'T' => 4, 'U' => 4)

numseqs = length(sequences)
numcols = length(sequences[1])
gapfrequency = zeros(Float64,numcols)
for s=1:numseqs
  seq = sequences[s]
  for j=1:numcols
    nuc = get(nucmapping,seq[j],0)
    if nuc > 0

    else
      gapfrequency[j] += 1.0
    end
  end
end
gapfrequency /= numseqs
println(gapfrequency)

fout = open(string(fastafile,".deletegaps.fas"), "w")
for s=1:numseqs
  seq = sequences[s]
  ret = ""
  for j=1:numcols
    if j <= 830 || j >= 10273 || gapfrequency[j] <= 0.5
      ret = string(ret, sequences[s][j])
    end
  end
  println(fout, ">", names[s])
  println(fout, ret)
end
close(fout)

include("MolecularEvolution.jl")
using MolecularEvolution

newickstring = "((AY366010:0.08102,((DQ878451:0.12030,AY771412:0.04363):0.01928,(AF347286:0.06751,DQ878674:0.12642):0.00453):0.01639):0.02989,DQ672623:0.11113,(AY787528:0.07062,(X52154:0.20926,EF517442:0.04835):0.02905):0.00445);"

root = gettreefromnewick(newickstring) # build tree object from newick string

println(MolecularEvolution.prettyprintstring(root)) # pretty print view of tree

println(length(root)) # number of children at `root`

roottree(root)
println(length(root)) # number of children at root after rooting.

# iterate every child of the root node, testing if a child is a leaf node
for child in root
  println(isleafnode(child))
end

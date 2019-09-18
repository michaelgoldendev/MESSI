dir = "../datasets/RFAM/3D/"

# start from RF00009

fileinfo = Tuple[]
for f in filter(x -> match(r".fas", x) != nothing, readdir(dir))
	filepath = joinpath(dir,f)
	push!(fileinfo, (filesize(filepath), filepath, f))
end
sort!(fileinfo)

for (index,info) in enumerate(fileinfo)
	println(index,"\t", info)
end
#exit()

for (sizeinbytes, filepath, f) in fileinfo
	#println("julia MESSI.jl --alignment $(filepath) --unpaired")
	#println("julia MESSI.jl --alignment $(filepath) --optimizecoevolutiononly --maxoptiter 400")
	println("julia MESSI.jl --alignment $(filepath)")
	println("julia MESSI.jl --alignment $(filepath) --maxstructure --gapcutoff 0.5")
	#println("julia MESSI.jl --alignment $(filepath) --numlambdacats 1 --optimizecoevolutiononly --maxoptiter 400")
	#println("julia MESSI.jl --alignment $(filepath) --numlambdacats 1")
	#println("julia MESSI.jl --alignment $(filepath) --numlambdacats 1 --maxstructure --gapcutoff 0.5")
end

#=
for f in filter(x -> match(r".fas", x) != nothing, readdir(dir))
	filepath = joinpath(dir,f)
	#println("julia MESSI.jl --alignment $(filepath) --unpaired")
	println("julia MESSI.jl --alignment $(filepath)")
	println("julia MESSI.jl --alignment $(filepath) --maxstructure")
end
println("")
=#
dir = "../datasets/RFAM/3D/"

fileinfo = Tuple[]
for f in filter(x -> match(r".fas", x) != nothing, readdir(dir))
	filepath = joinpath(dir,f)
	push!(fileinfo, (filesize(filepath), filepath, f))
end
sort!(fileinfo)

for (sizeinbytes, filepath, f) in fileinfo
	println("julia MESSI.jl --alignment $(filepath) --unpaired")
	println("julia MESSI.jl --alignment $(filepath)")
	println("julia MESSI.jl --alignment $(filepath) --maxstructure")
end

#=
for f in filter(x -> match(r".fas", x) != nothing, readdir(dir))
	filepath = joinpath(dir,f)
	#println("julia MESSI.jl --alignment $(filepath) --unpaired")
	#println("julia MESSI.jl --alignment $(filepath)")
	println("julia MESSI.jl --alignment $(filepath) --maxstructure")
end
println("")
=#
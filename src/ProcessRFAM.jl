include("RNATools.jl")

dir = "../datasets/RFAM/3D/"

for f in filter(x -> match(r".+\.stockholm\.txt", x) != nothing, readdir(dir))
	rfamid = match(r"(.+)\.stockholm\.txt", f)[1]
	filepath = joinpath(dir,f)
	
	fin = open(filepath)
	alignout = open(joinpath(dir,string(rfamid,".fas")),"w")
	lastseq = ""
	for ln in readlines(fin)
		line = strip(ln)

		if match(r"#=GC\s+SS_cons\s+(.+)", line) != nothing
			m = match(r"#=GC\s+SS_cons\s+(.+)", line)
			dotbracket = m[1]
			pairedsites = parsedbn(dotbracket)
			writectfile(pairedsites, lastseq, joinpath(dir,string(rfamid,".ct")))
		elseif startswith(line,"#")

		elseif match(r"(.+)\s+(.+)",line) != nothing
			m = match(r"(.+)\s+(.+)", line)
			println(alignout, ">", m[1])
			println(alignout, m[2])
			lastseq = m[2]
		end
	end	
	close(alignout)
	close(fin)
end

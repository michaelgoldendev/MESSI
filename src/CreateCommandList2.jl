filepaths = AbstractString[]

#push!(filepaths, "..\\datasets\\ssDNA\\beet_curly_mcmc\\bctv.fas.norm")
#push!(filepaths, "..\\datasets\\ssDNA\\beet_curly_mcmc\\bctv_seperated.fas")

#push!(filepaths, "..\\datasets\\ssDNA\\bocavirus_mcmc\\bocavirus.fas.norm")

#=
push!(filepaths, "..\\datasets\\ssDNA\\wdf_mcmc\\wdf.fas.norm")=#
#push!(filepaths, "..\\datasets\\ssDNA\\wdf_mcmc\\wdf_seperated.fas")

#push!(filepaths, "..\\datasets\\ssRNA\\rhinovirus_a_recombination\\rhinovirus_a_seperated_alignment.fas")



#push!(filepaths, "..\\datasets\\ncRNA\\RF00001_mcmc2\\RF00001_seperated.fas")
#push!(filepaths, "..\\datasets\\ncRNA\\RF00003_mcmc\\RF00003.fas.norm")
#push!(filepaths, "..\\datasets\\ncRNA\\RF00003_mcmc\\RF00003_seperated.fas")
#push!(filepaths, "..\\datasets\\ncRNA\\RF00010_mcmc\\RF00010.fas.norm")
#push!(filepaths, "..\\datasets\\ncRNA\\RF00010_mcmc\\RF00010_seperated.fas")
#push!(filepaths, "..\\datasets\\ncRNA\\RF00379_mcmc\\RF00379.fas.norm")
#push!(filepaths, "..\\datasets\\ncRNA\\RF00379_mcmc\\RF00379_seperated.fas")
#push!(filepaths, "..\\datasets\\ncRNA\\RF01846_mcmc\\RF01846.fas.norm")
#push!(filepaths, "..\\datasets\\ncRNA\\RF01846_mcmc\\RF01846_seperated.fas")
#push!(filepaths, "..\\datasets\\ssDNA\\bocavirus_mcmc\\bocavirus_seperated.fas")

#push!(filepaths, "..\\datasets\\ssDNA\\msv\\msv.fas.norm")
#push!(filepaths, "..\\datasets\\ssDNA\\msv\\msv_seperated.fas")

#push!(filepaths, "..\\datasets\\ssDNA\\tylcv\\tylcv.fas.norm")
#push!(filepaths, "..\\datasets\\ssDNA\\tylcv\\tylcv_seperated.fas")
push!(filepaths, "..\\datasets\\ssRNA\\rhinovirus_a_recombination\\rhinovirus_a_seperated_alignment.fas")
push!(filepaths, "..\\datasets\\ssRNA\\fmdv\\fmdv_seperated.fas")
push!(filepaths, "..\\datasets\\ssRNA\\hepa\\hepa_seperated.fas")
push!(filepaths, "..\\datasets\\ssRNA\\hpv1\\hpv1_seperated.fas")
push!(filepaths, "..\\datasets\\ssRNA\\tobamovirus\\tobamovirus_seperated.fas")


for filepath in filepaths
	if isfile(filepath)
		println("julia MESSI.jl --alignment $(filepath) --numlambdacats 1 --maxbasepairdistance 500 --mcmc")
		#=
		println("julia MESSI.jl --alignment $(filepath) --numlambdacats 1 --maxbasepairdistance 500 --unpaired")
		println("julia MESSI.jl --alignment $(filepath) --numlambdacats 1 --maxbasepairdistance 500 --optimizecoevolutiononly --maxoptiter 400")
		println("julia MESSI.jl --alignment $(filepath) --numlambdacats 1 --maxbasepairdistance 500")
		println("julia MESSI.jl --alignment $(filepath) --numlambdacats 1 --maxbasepairdistance 500 --fixgu --optimizecoevolutiononly --maxoptiter 300")
		println("julia MESSI.jl --alignment $(filepath) --numlambdacats 1 --maxbasepairdistance 500 --fixgu")	
		=#
	end
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
dir = "../datasets/RFAM/3D/"

alignmentpaths = []
push!(alignmentpaths, "..\\datasets\\ncRNA\\RF00003_mcmc\\RF00003_seperated.fas")
push!(alignmentpaths, "..\\datasets\\ncRNA\\RF00001_mcmc2\\RF00001_seperated.fas")
#push!(filepaths, "..\\datasets\\ncRNA\\RF00003_mcmc\\RF00003.fas.norm")
#push!(filepaths, "..\\datasets\\ncRNA\\RF00010_mcmc\\RF00010.fas.norm")
#push!(filepaths, "..\\datasets\\ncRNA\\RF00379_mcmc\\RF00379.fas.norm")
push!(alignmentpaths, "..\\datasets\\ncRNA\\RF00379_mcmc\\RF00379_seperated.fas")
#push!(filepaths, "..\\datasets\\ncRNA\\RF01846_mcmc\\RF01846.fas.norm")
push!(alignmentpaths, "..\\datasets\\ncRNA\\RF01846_mcmc\\RF01846_seperated.fas")
push!(alignmentpaths, "..\\datasets\\ncRNA\\RF00010_mcmc\\RF00010_seperated.fas")

#push!(alignmentpaths, "..\\datasets\\ssDNA\\wdf_mcmc\\wdf_seperated.fas")
#push!(alignmentpaths, "..\\datasets\\ssDNA\\beet_curly_mcmc\\bctv_seperated.fas")
#push!(alignmentpaths, "..\\datasets\\ssDNA\\bocavirus_mcmc\\bocavirus_seperated.fas")
#push!(alignmentpaths, "..\\datasets\\ssDNA\\msv\\msv_seperated.fas")
#push!(alignmentpaths, "..\\datasets\\ssDNA\\tylcv\\tylcv_seperated.fas")

#push!(filepaths, "..\\datasets\\ssDNA\\beet_curly_mcmc\\bctv.fas.norm")
#push!(filepaths, "..\\datasets\\ssDNA\\beet_curly_mcmc\\bctv_seperated.fas")

#push!(filepaths, "..\\datasets\\ssDNA\\bocavirus_mcmc\\bocavirus.fas.norm")


#push!(filepaths, "..\\datasets\\ssRNA\\rhinovirus_a_recombination\\rhinovirus_a_seperated_alignment.fas")


#alignmentpaths = ["../datasets/ssRNA/hepa/hepa.fas.norm"]
#alignmentpaths = ["../datasets/ssDNA/msv/msv.fas.norm"]

#=

numsims = 20
for alignmentpath in alignmentpaths
	m = match(r"^([^\.]+)\..+$", basename(alignmentpath))
	
	#println("julia MESSI.jl --alignment $(alignmentpath) --unpaired --maxbasepairdistance 500 --numlambdacats 1")
	#println("julia MESSI.jl --alignment $(alignmentpath) --fixgu  --maxbasepairdistance 500 --numlambdacats 1")
	#println("julia MESSI.jl --alignment $(alignmentpath)  --maxbasepairdistance 500 --numlambdacats 1")
	println("julia MESSI.jl --alignment $(alignmentpath)  --simulate $(numsims) --fixgu  --maxbasepairdistance 500 --numlambdacats 1")
	for sim=1:numsims
		outputprefix = "results/$(m[1])/simulations/sim$(sim)"
		println("julia MESSI.jl --alignment $(outputprefix).fas --tree $(outputprefix).nwk --outputprefix $(outputprefix) --fixgu  --maxbasepairdistance 500 --numlambdacats 1")
		println("julia MESSI.jl --alignment $(outputprefix).fas --tree $(outputprefix).nwk --outputprefix $(outputprefix) --maxbasepairdistance 500 --numlambdacats 1")
	end
end
exit()=#

using JSON
for alignmentpath in alignmentpaths
	m = match(r"^([^\.]+)\..+$", basename(alignmentpath))
	inputprefix = "results/$(m[1])/$(m[1])"
	
	inputmax = string(inputprefix,".max.sitecats3.lambdacats1")
	inputmaxfigu = string(inputprefix,".fixgu.max.sitecats3.lambdacats1")
	println(inputmax)
	if isfile(inputmax) && isfile(inputmaxfigu)
		inputmaxZ = JSON.parsefile(inputmax)["Z"]
		inputmaxguZ = JSON.parsefile(inputmaxfigu)["Z"]
		reallr = 2.0*(inputmaxguZ-inputmaxZ)
		println("real\t",reallr)
	end		 
	ls = Float64[]
	for sim=1:100000
		simulationsprefix = "results/$(m[1])/simulations/sim$(sim)"
		simmax = string("results/$(m[1])/simulations/sim$(sim)",".max.sitecats3.lambdacats1")
		simmaxfixgu = string("results/$(m[1])/simulations/sim$(sim)",".fixgu.max.sitecats3.lambdacats1")
		if isfile(simmax) && isfile(simmaxfixgu)			
			simmaxZ = JSON.parsefile(simmax)["Z"]
			simmaxguZ = JSON.parsefile(simmaxfixgu)["Z"]
			simlr = 2.0*(simmaxguZ-simmaxZ)
			println(sim,"\t",simlr)
			push!(ls, simlr)
		else
			break
		end
	end
	sort!(ls)
	println(ls)
	lower = 0
	upper = 1
	for z=1:length(ls)-1
		if ls[z] <= reallr < ls[z+1]
			lower = z
			upper = z + 1
		end
	end
	plower = lower / length(ls)
	pupper = upper / length(ls)
	println(plower," <= p <= ",pupper)
	println((pupper+plower)/2.0)
end
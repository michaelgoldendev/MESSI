function felsensteinstack(nodelist::Array{TreeNode,1}, likelihoods::Array{Float64,2}, logm::Array{Float64,1}, transprobs::Array{Array{Float64,2},1}, alphabet::Int)
	for node in nodelist
	    if !isleafnode(node)
	        for a=1:alphabet
	            likelihoods[node.nodeindex,a] = -Inf
	        end
	    end
	end
  stack = Int[1]
  while length(stack) > 0
    nodeindex = stack[end]
    node = nodelist[nodeindex]
    if isleafnode(node)
        pop!(stack)
    else
        leftchildindex = node.children[1].nodeindex
        rightchildindex = node.children[2].nodeindex

        cont = true
        if likelihoods[leftchildindex, 1] == -Inf
          push!(stack, leftchildindex)
          cont = false
        end
        if likelihoods[rightchildindex, 1] == -Inf
          push!(stack, rightchildindex)
          cont = false
        end

        if cont
          v = (transprobs[leftchildindex]*likelihoods[leftchildindex,:]).*(transprobs[rightchildindex]*likelihoods[rightchildindex,:])
          likelihoods[nodeindex,:] = (transprobs[leftchildindex]*likelihoods[leftchildindex,:]).*(transprobs[rightchildindex]*likelihoods[rightchildindex,:])
          m = maximum(likelihoods[nodeindex,:])
          if m < 1e-20
            likelihoods[nodeindex,:] /= m
            logm[nodeindex] = log(m) + logm[leftchildindex] + logm[rightchildindex]
          else
            logm[nodeindex] = logm[leftchildindex] + logm[rightchildindex]
          end
          pop!(stack)
        end
    end
  end

  return likelihoods
end

function backwardssampling(rng::AbstractRNG, nodelist::Array{TreeNode,1}, likelihoods::Array{Float64,2}, logm::Array{Float64,1}, transprobs::Array{Array{Float64,2},1}, freqs::Array{Float64,1}, nsamples::Int=1)
    rootliks = likelihoods[1,:].*freqs
    rootliks /= sum(rootliks)
    samples = zeros(Int, nsamples, length(nodelist))
    for s=1:nsamples
        samples[s, 1] = CommonUtils.sample(rng, rootliks)
        stack = Int[1]
        while length(stack) > 0
          nodeindex = stack[end]
          node = nodelist[nodeindex]
          parentc = samples[s, nodeindex]
          pop!(stack)
          for child in node.children
              if samples[s, child.nodeindex] == 0
                  samples[s, child.nodeindex] = CommonUtils.sample(rng, transprobs[child.nodeindex][parentc,:].*likelihoods[child.nodeindex,:])
                  if !isleafnode(child)
                      push!(stack, child.nodeindex)
                  end
              end
          end
      end
    end
    return samples
end

function simulatealignments(rng::AbstractRNG, dataset::Dataset, inparams::ModelParameters, outputdir::AbstractString, siteCats::Int=3, lambdacats::Int=5, parameterisation::Int=1; samples::Int=1,usecuda::Bool=true, maxbasepairdistance::Int=1000000, fixGU::Bool=false,fixGCAU::Bool=false,unpairedmodel::Bool=false)
	paramvector = getparamsvector(inparams)
	if fixGU
		paramvector[3] = 1.0
	end
	if fixGCAU
		paramvector[2] = paramvector[1]
	end
	if unpairedmodel
		paramvector[1] = 1.0
		paramvector[2] = 1.0
		paramvector[3] = 1.0
	end
	params = getparams(paramvector, dataset, siteCats, lambdacats, parameterisation,fixGU,fixGCAU)
	println(inparams)
	println(params)
	if !isdir(outputdir)
		mkpath(outputdir)
	end
	grammar = KH99()
    unpairedlogprobs = zeros(Float64, dataset.numcols)
    pairedlogprobs = zeros(Float64, dataset.numcols, dataset.numcols)

    maskgapped!(pairedlogprobs,dataset.gapfrequency,0.5,-Inf)
    inside = computeinsideKH99(unpairedlogprobs, pairedlogprobs, 1.0, false, usecuda)
    Z = inside[1,1,dataset.numcols]

	unpairedfreqlists = Array{Float64,1}[]
    unpairedtransprobslists = Array{Array{Float64,2},1}[]
    for rate in params.siteRates
    	push!(unpairedfreqlists, params.freqs)
	  	Q = gtr(params.q1, params.q2, params.q3, params.q4, params.q5, params.q6, params.freqs)
  		transprobs = gettransitionmatriceslist(params.branchlengths, Q*rate)
  		push!(unpairedtransprobslists, transprobs)
  	end
  	unpairedcatprobs = params.siteWeights


    museparams = getmusespecificparamsarray(params)
    musecatprobs = Float64[]
    musefreqlists = Array{Float64,1}[]
    musetransprobslists = Array{Array{Float64,2},1}[]
    for siteCat1=1:params.siteCats
	    for siteCat2=1:params.siteCats
	        for musespecificparams in museparams
	            musemodel = MuseModel(params.freqs, getGC(params,musespecificparams), getAT(params,musespecificparams), getGT(params,musespecificparams), params.q1, params.q2, params.q3, params.q4, params.q5, params.q6, params.siteRates[siteCat1], params.siteRates[siteCat2])	            
	            logw = log(params.siteWeights[siteCat1]) + log(params.siteWeights[siteCat2]) + musespecificparams.logprob
	            push!(musecatprobs, logw)

	            push!(musefreqlists, musemodel.freqs)
	            transprobslist = gettransitionmatriceslist(params.branchlengths, musemodel.Q)
	            push!(musetransprobslists, transprobslist)
	        end
	    end
	end
	musecatprobs = exp.(musecatprobs .- maximum(musecatprobs))
	musecatprobs /= sum(musecatprobs)
	
	println(unpairedcatprobs)
	println(musecatprobs)
	nodelist = getnodelist(dataset.root)

	#pairedlikelihoods = ones(Float64, length(nodelist), 16) 
    #pairedlogm = zeros(Float64, length(nodelist))
    #felsensteinstack(nodelist, pairedlikelihoods, pairedlogm, musetransprobslists[1], 16)
    nucleotides = "ACGT"
    for s=1:samples
        paired = zeros(Int,dataset.numcols)
        samplestructure(rng, inside, pairedlogprobs, unpairedlogprobs, 1, dataset.numcols, paired, grammar, 1.0)

        alignment = zeros(Int, dataset.numseqs, dataset.numcols)

        for i=1:length(paired)
	        if paired[i] == 0
	        	z = CommonUtils.sample(rng, unpairedcatprobs)
	        	samples = backwardssampling(rng, nodelist, ones(Float64, length(nodelist), 4), zeros(Float64, length(nodelist)), unpairedtransprobslists[z], unpairedfreqlists[z])
	        	for node in nodelist
	        		if node.seqindex > 0
	        			alignment[node.seqindex,i] = samples[1,node.nodeindex]
	        		end
	        	end
	        elseif paired[i] > i
	        	z = CommonUtils.sample(rng, musecatprobs)
	        	samples = backwardssampling(rng, nodelist, ones(Float64, length(nodelist), 16), zeros(Float64, length(nodelist)), musetransprobslists[z], musefreqlists[z])
	        	for node in nodelist
	        		if node.seqindex > 0
		        		alignment[node.seqindex,i] = ((samples[1,node.nodeindex]-1) % 4) + 1
		        		alignment[node.seqindex,paired[i]] = div(samples[1,node.nodeindex]-1, 4) + 1
		        	end
	        	end
	        end
	    end

        sequences = AbstractString[]
        for seq=1:dataset.numseqs
        	sequence = ""
        	for col=1:dataset.numcols
        		if dataset.sequences[seq][col] == '-'
        			sequence = string(sequence, '-')        		
        		else
        			sequence = string(sequence, nucleotides[alignment[seq,col]])        		
        		end
        	end
        	push!(sequences,sequence)
        end
        fout = open(joinpath(outputdir, "sim$(s).fas"), "w")
        for (seqindex,seq) in enumerate(sequences)
        	println(fout, ">seq$(seqindex)")
        	println(fout, seq)
        end
        close(fout)
    end
end
using Pkg
#=
Pkg.add(PackageSpec(name="ArgParse", uuid="c7e460c6-2fb9-53a9-8c5b-16f535851c63"))
Pkg.add(PackageSpec(name="CSV", uuid="336ed68f-0bac-5ca0-87d4-7b16caf5d00b"))
Pkg.add(PackageSpec(name="CUDAdrv", uuid="c5f51814-7f29-56b8-a69c-e4d8f6be1fde"))
Pkg.add(PackageSpec(name="DataFrames", uuid="a93c6f00-e57d-5684-b7b6-d8193f3e46c0"))
Pkg.add(PackageSpec(name="DataStructures", uuid="864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"))
Pkg.add(PackageSpec(name="Distributions", uuid="31c24e10-a181-5473-b8eb-7969acd0382f"))
Pkg.add(PackageSpec(name="FastaIO", uuid="a0c94c4b-ebed-5953-b5fc-82fe598ac79f"))
Pkg.add(PackageSpec(name="Formatting", uuid="59287772-0a20-5a39-b81b-1366585eb4c0"))
Pkg.add(PackageSpec(name="HypothesisTests", uuid="09f84164-cd44-5f33-b23f-e6b0d136a0d5"))
Pkg.add(PackageSpec(name="JSON", uuid="682c06a0-de6a-54ab-a142-c8b1cf79cde6"))
Pkg.add(PackageSpec(name="NLopt", uuid="76087f3c-5699-56af-9a33-bf431cd00edd"))
Pkg.add(PackageSpec(name="Nullables", uuid="4d1e1d77-625e-5b40-9113-a560ec7a8ecd"))
Pkg.add(PackageSpec(name="PyPlot", uuid="d330b81b-6aea-500a-939a-2ce795aea3ee"))
Pkg.add(PackageSpec(name="StatsBase", uuid="2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"))
Pkg.add(PackageSpec(name="Blosc", uuid="a74b3585-a348-5f62-a45c-50e91977d574"))
Pkg.add(PackageSpec(name="HDF5", uuid="f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"))
Pkg.add(PackageSpec(name="SHA", uuid="ea8e919c-243c-51af-8825-aaa63cd721ce"))
Pkg.add("CUDAnative")
Pkg.add("CuArrays")
=#

push!(LOAD_PATH,joinpath(@__DIR__))
using Binaries

include("FARCE.jl")
include("MCMC.jl")
include("Simulations.jl")
using DataFrames
using NLopt
include("Covariation.jl")
using Mapping
include("RankingAndVisualisation.jl")
#using ProfileView
#using GaussianProcesses
#include("BayesianOptimization.jl")
#using BayesianOptimizationCustom
#@everywhere include("InsideParallel.jl")
using Distributions
using Printf
using CommonUtils
using RNABinaries
using SparseArrays
using BayesianOptimization
using GaussianProcesses

include("StructurePlotsAndTables.jl")

using ArgParse
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--alignment"
        help = "alignment file to be used"
        arg_type = String
        required = true
        "--tree"
        help = "tree file to be used"
        arg_type = String
        required = false
        "--structure"
        help = "RNA structure to be used"
        arg_type = String
        required = false
        "--numsitecats"
        help = "number of rate categories to be used for site-to-site rate variation"
        arg_type = Int
        default = 3
        "--numlambdacats"
        help = "number of rate categories to be used for base-pair to base-pair rate variation"
        arg_type = Int
        default = 5
        "--outputprefix"
        help = "prefix to be used for outputting results files"
        arg_type = String
        "--maxbasepairdistance"
        help = ""
        arg_type = Int
        default = 100000000
        "--gapcutoff"
        help = ""
        arg_type = Float64
        default = 1.0
        "--M"
        help = ""
        arg_type = Float64
        default = 1.0
        "--draw"
        help = ""
        action = :store_true
        "--full"
        help = "Perform all Maximum Likelihood inference under all models and output Likelihood Ratio Tests"
        action = :store_true
        "--mcmc"
        help = "Perform Bayesian inference using MCMC instead of the default Maximum Likelihood inference"
        action = :store_true
        "--simulate"
        help = ""
        arg_type = Int
        default = 0
        #=
        "--optimize"
        help = ""
        action = :store_true
        =#
        "--samplesperiter"
        help = ""
        arg_type = Int
        default = 0
        "--maxoptiter"
        help = ""
        arg_type = Int
        default = 2000
        "--maxmcemiter"
        help = ""
        arg_type = Int
        default = 150
        "--shuffle"
        help = ""
        action = :store_true
        "--unpaired"
        help = ""
        action = :store_true
        "--seed"
        help = ""
        arg_type = Int
        default = 15039814119191
        "--cpuonly"
        help = ""
        action = :store_true
        "--fixgu"
        help = ""
        action = :store_true
        "--fixgcau"
        help = ""
        action = :store_true
        "--optimizecoevolutiononly"
        help = ""
        action = :store_true
        "--calcentropy"
        help = ""
        action = :store_true
        "--calcpriorentropy"
        help = ""
        action = :store_true
        "--maxstructure"
        help = ""
        action = :store_true
        "--processmax"
        help = ""
        action = :store_true
        "--lambdazeroweightmin"
        help = ""
        arg_type = Float64
        default = 0.05
        "--bfactor"
        help = ""
        arg_type = Float64
        default = 0.98
        "--numchains"
        help = ""
        arg_type = Int
        default = 1
        "--expcanonicalconstraints"
        help = "EXPERIMENTAL. Only consider pairs of sites with at least one GC, AU, or GU pair."
        action = :store_true
        "--expcanonicalrnafold"
        help = "EXPERIMENTAL. Fold sequences at specified temperature in degrees celsius and use base-pair probabilities to constrain which sites may pair."
        arg_type = Float64
        default = -1000.0
    end

    return parse_args(s)
end

lambdaweightprior = Beta(2.0,2.0)
function printmatrix(mat::Array{Float64,2})
    ret = ""
    for i=1:size(mat,1)
        ret = string(ret,join(AbstractString[string(@sprintf("%.4f", v)) for v in mat[i,:]], "\t"), "\n")
    end
    return ret
end

function mapstructure(dataset, fastafile, ctfile)
    sequence, pairedsites = readctfile(ctfile)
    mapping, revmapping = createmapping(fastafile, sequence)
    len = length(pairedsites)
    mappedsequence = ""
    mappedpairedsites = zeros(Int,dataset.numcols)
    for i=1:dataset.numcols
        j = get(mapping, i, 0)
        println("A",i,"\t",j)
        if j > 0 && pairedsites[j] > j
            k = get(revmapping, pairedsites[j], 0)
            println("B",pairedsites[j],"\t",k)
            if k > 0
                mappedpairedsites[i] = k
                mappedpairedsites[k] = i
            end
        end

        if j  > 0
            mappedsequence = string(mappedsequence, sequence[j])
        else
            mappedsequence = string(mappedsequence, "-")
        end
    end

    dbnstring = ""
    for i=1:length(mappedpairedsites)
        if mappedpairedsites[i] > i
            dbnstring = string(dbnstring,"(")
        elseif mappedpairedsites[i] == 0
            dbnstring = string(dbnstring,".")
        else
            dbnstring = string(dbnstring,")")
        end
    end
    println(mappedsequence)
    println(dbnstring)

    return mappedsequence,mappedpairedsites, mapping, revmapping
end

function getZ(params::Array{Float64,1}, dataset::Dataset, siteCats::Int=3, lambdaCats::Int=5, usecuda::Bool=true,unpaired::Bool=true;gapcutoff::Float64=1.0)
    #tic()
    currentparams = getparams(params,dataset,siteCats,lambdaCats,0)
    if unpaired
        currentparams.lambdaratesGC = zeros(Float64, length(currentparams.lambdaratesGT))
        currentparams.lambdaratesAT = zeros(Float64, length(currentparams.lambdaratesGT))
        currentparams.lambdaratesGT = zeros(Float64, length(currentparams.lambdaratesGT))
    end

    grammar = KH99()
    unpairedlogprobs = computeunpairedlikelihoods(dataset, currentparams)
    maxbasepairdistance = 1000000
    println("C. maxbasepairdistance: ", maxbasepairdistance)
    pairedlogprobs, ret = coevolutionall(dataset,currentparams,true,false,true,maxbasepairdistance,usecuda)
    #elapsed = toc()
    #println("P1\t", elapsed)
    maskgapped!(pairedlogprobs,dataset.gapfrequency,gapcutoff,-Inf)
    #tic()
    ll = computeinsideKH99(unpairedlogprobs, pairedlogprobs, 1.0,true,usecuda)
    println(ll)
    #elapsed = toc()
    #println("P2\t", elapsed)
    return ll
end


maxll = -1e20
maxparams = []
optiter = 0
function modellikelihood(inparams::Array{Float64,1}, dataset::Dataset, paired::Array{Int,1}, siteCats::Int=3, lambdaCats::Int=5, fixGU::Bool=false, fixGCAU::Bool=false, integratestructure::Bool=true, unpairedmodel::Bool=false, maxfile=nothing, usecuda::Bool=true; maxbasepairdistance::Int=1000000,gapcutoff=1.0)
    global maxll
    global maxparams
    global optiter
    optiter += 1

    ll = -1e20

    params = copy(inparams)
    #println(optiter,"\t",params)
    #=
    if unpairedmodel      
        params[1] = 1.0
        params[2] = 1.0
        params[3] = 1.0
    end
    if fixGU
        params[3] = 1.0
    end=#
    if fixGCAU
        params[2] = params[1]
    end
  
    #=
    lower,upper,dummy = getbounds(fixGU, false, unpairedmodel)
    for z=1:length(params) # correct any out of bounds values
      params[z] = max(lower[z], params[z])
      params[z] = min(upper[z], params[z])
    end
    =#
    if sum(params[12:14]) >= 0.95
        return -Inf
    end
    currentparams = getparams(params,dataset,siteCats,lambdaCats,0,fixGU,fixGCAU)

    if unpairedmodel
        currentparams.lambdaratesGC = zeros(Float64, length(currentparams.lambdaratesGT))
        currentparams.lambdaratesAT = zeros(Float64, length(currentparams.lambdaratesGT))
        currentparams.lambdaratesGT = zeros(Float64, length(currentparams.lambdaratesGT))
    end

    if fixGU
        currentparams.lambdaratesGT = zeros(Float64, length(currentparams.lambdaratesGT))
    end

    try
        if integratestructure
            grammar = KH99()
            unpairedlogprobs = computeunpairedlikelihoods(dataset, currentparams)
            #tic()
            #maxbasepairdistance = 1000000
            #println("E. maxbasepairdistance: ", maxbasepairdistance)
            pairedlogprobs, ret = coevolutionall(dataset,currentparams,true,false,true,maxbasepairdistance,usecuda)
            maskgapped!(pairedlogprobs,dataset.gapfrequency,gapcutoff,-Inf)
            ll = computeinsideKH99(unpairedlogprobs, pairedlogprobs, 1.0, true, usecuda)
            #elapsed = toc();
            #println("Iteration ", optiter, ", elapsed: ", elapsed)
            #println("Iteration ", optiter)
            #println("A",ll)
        else
            #maxbasepairdistance = 1000000
            #println("F. maxbasepairdistance: ", maxbasepairdistance)
            ll = computetotallikelihood(MersenneTwister(5043820111),dataset, currentparams, paired,false,true,false,1.0,false, maxbasepairdistance, usecuda)
        end

        if !unpairedmodel && lambdaCats > 1
            ll += logpdf(lambdaweightprior, params[15])
        else
            ll += logpdf(lambdaweightprior, 0.5)
        end

        savemaximum(ll, params, maxfile)
        m = match(r"^([^\.]+)\..+$", basename(maxfile))
        fout = open(joinpath(dirname(maxfile),string("$(m[1]).cache.list",".sitecats$(siteCats).lambdaCats$(lambdaCats)")),"a")
        println(fout, ll,"\t",params)
        close(fout)
    catch e
        println("Exception")
        println(e)
        println(stacktrace(catch_backtrace()))
        if occursin("CuError", string(e))
            exit()
        end
        return -1e20
        #return ll
    end

    if isnan(ll) || isinf(ll)
        return -1e20
    end

    if ll > maxll
        maxll = ll
        maxparams = params
        ret = "Maximum\t"
        ret = string(ret, maxll, "\t")
        ret = string(ret, @sprintf("%d", optiter), "\t")
        for x in params
            ret = string(ret, @sprintf("%0.4f", x), "\t")
        end
        println(ret)
    end

    return ll
end


function freqconstraint(x::Array{Float64,1})
    return (x[12] + x[13] + x[14] - 1.0)
    #return 1.0 - (x[12] + x[13] + x[14])
end


function computesamplelikelihoods(params::ModelParameters, dataset::Dataset, structures::Array{Array{Int,1},1})
    total = -Inf
    museparams = getmusespecificparamsarray(params)
    len = params.siteCats*params.siteCats*length(params.lambdaratesGC)
    cache = FastCache[FastCache(250000,16) for i=1:len]
    unpairedcache = FastCache[FastCache(100000,4) for i=1:params.siteCats]
    unpairedcolumncache = ColumnCache(100000)
    pairedcolumncache = ColumnCache(1000000)
    likelihoods = Float64[]
    for paired in structures
        ll = calculatecachedlikelihood(dataset, params, museparams, params.states, paired, cache, unpairedcache, true,unpairedcolumncache,pairedcolumncache)
        push!(likelihoods, ll)
    end
    return likelihoods
end

#=
function computesamplelikelihoods(params::ModelParameters, dataset::Dataset, structures::Array{Array{Int,1},1})
likelihoods = SharedArray(Float64,length(structures))

numthreads = 8
blocks = Array{Tuple{Int,Array{Int,1}},1}[Tuple{Int,Array{Int,1}}[] for p=1:numthreads]
for i=1:length(structures)
push!(blocks[(i-1) % numthreads+1], (i,structures[i]))
end

museparams = getmusespecificparamsarray(params)
len = params.siteCats*params.siteCats*length(params.lambdarates)

@sync @parallel for block in blocks
cache = FastCache[FastCache(250000,16) for i=1:len]
unpairedcache = FastCache[FastCache(100000,4) for i=1:params.siteCats]
unpairedcolumncache = ColumnCache(100000)
pairedcolumncache = ColumnCache(1000000)
for structurepair in block
i = structurepair[1]
paired = structurepair[2]
likelihoods[i] = calculatecachedlikelihood(dataset, params, museparams, params.states, paired, cache, unpairedcache, true, unpairedcolumncache, pairedcolumncache)
end
end
return likelihoods
end=#

function computeimportanceratio(currentlikelihoods, params::Array{Float64,1}, dataset::Dataset, structures::Array{Array{Int,1},1}, siteCats::Int=3, lambdaCats::Int=5, fixGU::Bool=false,fixGCAU::Bool=false)
    println("X", params)
    global maxll
    global maxparams
    global optiter
    optiter += 1

    #tic()

    llsum = -Inf
    if sum(params[12:14]) >= 0.99 || params[10] < 0.0 || params[11] < 0.0
        return -1e20
    end
    freqs = Float64[params[12],params[13],params[14],1.0 - sum(params[12:14])]
    freqs /= sum(freqs)

    proposedparams = getparams(params,dataset,siteCats,lambdaCats,0,fixGU,fixGCAU)
    if fixGU
        proposedparams.lambdaratesGT = zeros(Float64, length(proposedparams.lambdaratesGT))
    end
    if fixGCAU
        proposedparams.lambdaratesGC = proposedparams.lambdaratesAT
    end
    proposedlikelihoods = computesamplelikelihoods(proposedparams, dataset, structures)

    count = 0
    for (currentll, proposedll) in zip(currentlikelihoods, proposedlikelihoods)
        llsum = CommonUtils.logsumexp(llsum, proposedll-currentll)
        count += 1
    end
    #println(llsum, "\t", log(count))
    llsum = llsum - log(count)
    #println("Importance ratio=", llsum)

    if isnan(llsum) || isinf(llsum)
        return -1e20
    end

    if llsum > maxll
        maxll = llsum
        maxparams = params
        #println("Maximum\t",maxll,"\t",optiter,"\t", maxparams)
        ret = "Maximum\t"
        ret = string(ret, maxll, "\t")
        #ret = string(ret, @sprintf("%0.2f", maxll), "\t")
        ret = string(ret, @sprintf("%d", optiter), "\t")
        for x in maxparams
            ret = string(ret, @sprintf("%0.2f", x), "\t")
        end
        println(ret)
        currentparams = getparams(params,dataset,siteCats,lambdaCats,0,fixGU,fixGCAU)
        if fixGU
            currentparams.lambdaratesGT = zeros(Float64, length(currentparams.lambdaratesGT))
        end
        if fixGCAU
            currentparams.lambdaratesGC = currentparams.lambdaratesAT
        end
        #=
        meanlambdaGC = sum((currentparams.lambdaratesGC+1.0).*currentparams.lambdaweightsGC)
        meanlambdaAT = sum((currentparams.lambdaratesAT+1.0).*currentparams.lambdaweightsAT)
        meanlambdaGT = sum((currentparams.lambdaratesGT+1.0).*currentparams.lambdaweightsGT)
        meanlambdaGC2 = sum((currentparams.lambdaratesGC[2:end]).*currentparams.lambdaweightsGC[2:end])/sum(currentparams.lambdaweightsGC[2:end])
        meanlambdaAT2 = sum((currentparams.lambdaratesAT[2:end]).*currentparams.lambdaweightsAT[2:end])/sum(currentparams.lambdaweightsAT[2:end])
        meanlambdaGT2 = sum((currentparams.lambdaratesGT[2:end]).*currentparams.lambdaweightsGT[2:end])/sum(currentparams.lambdaweightsGT[2:end])
        println("BmeanlambdaGC", meanlambdaGC,"\t", meanlambdaGC2,"\t", currentparams.lambdaGammaShapeGC, "\t", currentparams.lambdaGammaScaleGC, "\t", currentparams.lambdaratesGC)
        println("BmeanlambdaAT", meanlambdaAT,"\t", meanlambdaAT2,"\t", currentparams.lambdaGammaShapeAT, "\t", currentparams.lambdaGammaScaleAT, "\t", currentparams.lambdaratesAT)
        println("BmeanlambdaGT", meanlambdaGT,"\t", meanlambdaGT2,"\t", currentparams.lambdaGammaShapeGT, "\t", currentparams.lambdaGammaScaleGT, "\t", currentparams.lambdaratesGT)=#
    end
    #toc()
    return llsum
end

function parsejsonfile(infile)
    jsondict = Dict()
    open(infile,"r") do file
        jsondict = JSON.parse(file)
    end
    return jsondict
end

function savemaximum(Z::Float64, maxparams::Array{Float64,1}, maxfile, override::Bool=false, statuslabel::String="")
    maximumZ  = -Inf
    maximumparams = maxparams

    if isfile(maxfile)
        jsondict = parsejsonfile(maxfile)
        if Z > jsondict["Z"] || override
            if Z > maximumZ  || override
                maximumZ = Z
                jsondict["Z"] = maximumZ
                jsondict["maxparams"] = maxparams
                jsondict["status"] = statuslabel
                maximumparams = jsondict["maxparams"]
                out = open(maxfile,"w")
                ret = replace(JSON.json(jsondict),",\"" => ",\n\"")
                ret = replace(ret, "],[" => "],\n[")
                ret = replace(ret, "{" => "{\n")
                ret = replace(ret, "}" => "\n}")
                write(out,ret)
                close(out)
            end
        end
    elseif Z > maximumZ || override
        maximumZ = Z
        jsondict = Dict()
        jsondict["Z"] = maximumZ
        jsondict["maxparams"] = maxparams
        jsondict["status"] = statuslabel
        maximumparams = jsondict["maxparams"]
        out = open(maxfile,"w")
        ret = replace(JSON.json(jsondict),",\"" => ",\n\"")
        ret = replace(ret, "],[" => "],\n[")
        ret = replace(ret, "{" => "{\n")
        ret = replace(ret, "}" => "\n}")
        write(out,ret)
        close(out)
    end

    return maximumZ, maximumparams
end


function optimizesamplelikelihood(rng::AbstractRNG, initialparams::Array{Float64,1}, dataset::Dataset, siteCats::Int=3, lambdaCats::Int=5, fixGU::Bool=false, fixGCAU::Bool=false, fixLambdaWeight::Bool=true, unpairedmodel::Bool=false,maxoptiter::Int=1000,samplesperiter::Int=50,maxbasepairdistance::Int=500, maxfile=nothing, usecuda::Bool=true)
    global maxll
    maxll = -1e20
    if unpairedmodel
        initialparams[1] = 1.0
        initialparams[2] = 1.0
        initialparams[3] = 1.0
        #initialparams[10] = 1.0
    end
    if fixGU
        initialparams[3] = 1.0
    end
    if fixLambdaWeight
        initialparams[15] = 0.5
    end
    initialparams[15] = max(0.0001, initialparams[15])

    currentparams = getparams(initialparams,dataset,siteCats,lambdaCats,0,fixGU,fixGCAU)
    if fixGU
        currentparams.lambdaratesGT = zeros(Float64, length(currentparams.lambdaratesGT))
    end
    if fixGCAU
        currentparams.lambdaratesGC = currentparams.lambdaratesAT
    end

    grammar = KH99()
    unpairedlogprobs = computeunpairedlikelihoods(dataset, currentparams)
    #tic()
    println("F. maxbasepairdistance: ", maxbasepairdistance)
    pairedlogprobs, ret = coevolutionall(dataset,currentparams,true,false,true,maxbasepairdistance, usecuda)
    #println("unpairedlogprobs",unpairedlogprobs)
    #println("pairedlogprobs",pairedlogprobs)
    #elapsed = toc()
    #println("P1\t", elapsed)
    maskgapped!(pairedlogprobs,dataset.gapfrequency,gapcutoff,-Inf)

    #tic()
    inside = computeinsideKH99(unpairedlogprobs, pairedlogprobs, 1.0, false, usecuda)
    #toc()
    Z = inside[1,1,dataset.numcols]
    if maxfile != nothing
        savemaximum(Z,initialparams,maxfile, false, "Optimisation incomplete")
    end


    structures = Array{Int,1}[]
    for i=1:samplesperiter
        paired = zeros(Int,dataset.numcols)
        samplestructure(rng, inside, pairedlogprobs, unpairedlogprobs, 1, dataset.numcols, paired, grammar, 1.0)
        push!(structures,paired)
    end
    currentlikelihoods = computesamplelikelihoods(currentparams, dataset, structures)

    #opt = Opt(:LN_NELDERMEAD, 20)
    opt = Opt(:LN_COBYLA, 20)
    localObjectiveFunction = ((param, grad) -> computeimportanceratio(currentlikelihoods, param, dataset, structures, siteCats, lambdaCats, fixGU, fixGCAU))
    lower = ones(Float64, 20)*0.0001
    lower[1] = 1.0
    lower[2] = 1.0
    lower[3] = 1.0
    #lower[10] = 0.0001
    lower[10] = 1.0 #GC
    lower[11] = 0.0001
    lower[12] = 0.01
    lower[13] = 0.01
    lower[14] = 0.01
    lower[15] = 0.0001
    if unpairedmodel
        lower[10] = 1.0
        lower[15] = 0.5
        #lower[11] = 1.0
    end
    lower[16] = 0.0001
    #lower[17] = 0.0001
    lower[17] = 1.0 # AT
    lower[18] = 0.0001
    lower[19] = 1.0 # GT
    #lower[19] = 0.0001
    lower[20] = 0.0001

    if fixLambdaWeight
        lower[15] = 0.5
    end
    lower_bounds!(opt, lower)

    upper = ones(Float64, 20)*50.0
    upper[1] = 1.0
    upper[2] = 1.0
    upper[3] = 1.0
    upper[12] = 2.0
    upper[13] = 2.0
    upper[14] = 2.0
    upper[15] = 0.9999
    upper[16] = 250.0
    if unpairedmodel
        upper[1] = 1.0
        upper[2] = 1.0
        upper[3] = 1.0
        upper[10] = 1.0
        upper[15] = 0.5
        #upper[11] = 1.0
    end

    if fixLambdaWeight
        upper[15] = 0.5
    end
    if fixGU
        upper[3] = lower[3]
        upper[19] = lower[19]
    end
    if fixGCAU
        upper[2] = lower[2]
        upper[10] = lower[10] # GC
    end
    #ZZZZZZZZZZZZZZZZZZ
    upper[18] = 0.0001
    upper[20] = 0.0001
    upper_bounds!(opt, upper)

    ftol_abs!(opt,1e-3)
    maxeval!(opt, maxoptiter)
    max_objective!(opt, localObjectiveFunction)


    (minf,minx,ret) = optimize(opt, initialparams)
    return minx, initialparams, Z + logpdf(lambdaweightprior, minx[15])
end

function mcem(dataset::Dataset, siteCats::Int=3, lambdaCats::Int=5, fixGU::Bool=false, fixGCAU::Bool=false, unpairedmodel::Bool=false,maxmcemiter::Int=150,samplesperiter::Int=50,maxbasepairdistance::Int=500, initparams::ModelParameters=nothing, maxfile=nothing, usecuda::Bool=true)
    global maxll
    maxll = -1e20
    fixLambdaWeight = false
    initialparams = 2.0*ones(Float64,20)
    if initialparams == nothing
        initialparams[1] = 1.5
        initialparams[2] = 1.5
        initialparams[3] = 1.5
        initialparams[12] = dataset.obsfreqs[1]
        initialparams[13] = dataset.obsfreqs[2]
        initialparams[14] = dataset.obsfreqs[3]
        initialparams[15] = 0.5
        initialparams[1] = 1.0
        initialparams[2] = 1.0
        initialparams[3] = 1.0
    else
        initialparams = getparamsvector(initparams)
    end

    #=
    if fixGU
    initialparams[3] = 1.0
end
if fixGCAU
initialparams[2] = initialparams[1]
end=#
if fixLambdaWeight
    initialparams[15] = 0.5
end

initialparams[15] = max(0.0001, initialparams[15])

rng = MersenneTwister(757494371317)
maxZ = -1e20
maxparams = nothing
noimprovement = 0
maxoptiter = 100
for i=1:maxoptiter
    initialparams, params, Z = optimizesamplelikelihood(rng, initialparams, dataset, siteCats, lambdaCats, fixGU, fixGCAU, fixLambdaWeight, unpairedmodel, maxmcemiter,samplesperiter,maxbasepairdistance, maxfile, usecuda)
    if Z > maxZ
        maxZ = Z
        maxparams = copy(params)
        noimprovement = 0
    else
        noimprovement += 1
    end
    println(i,"\t",maxZ,"\t", noimprovement,"\t", maxparams)

    if noimprovement >= 3
        break
    end
end

return getparams(maxparams, dataset, siteCats, lambdaCats,0,fixGU,fixGCAU)
end

function shufflealignment(rng::AbstractRNG, fastafile, outputfile)
    sequences = AbstractString[]
    names = AbstractString[]

    FastaIO.FastaReader(fastafile) do fr
        for (desc, seq) in fr
            len = length(seq)
            push!(names,desc)
            push!(sequences, seq)
        end
    end

    columnindices = Int[i for i=1:length(sequences[1])]
    shuffle!(rng, columnindices)

    newsequences = AbstractString[]
    for seq in sequences
        newseq = ""
        for col in columnindices
            newseq = string(newseq, seq[col])
        end
        push!(newsequences, newseq)
    end

    fout = open(outputfile, "w")
    seqindex = 1
    for seq in newsequences
        write(fout, string(">seq", seqindex, "\n"))
        write(fout, string(seq, "\n"))
        seqindex += 1
    end
    close(fout)
    return outputfile
end

function maskarraysel(arr::Array{Float64,1}, sel::Array{Int,1})
    outarr = Float64[]
    for s in sel
        push!(outarr, arr[s])
    end
    return outarr
end

function maskarraysel(arr::Array{StepRangeLen{Float64},1}, sel::Array{Int,1})
    outarr = StepRangeLen{Float64}[]
    for s in sel
        push!(outarr, arr[s])
    end
    return outarr
end

function unmaskarraysel(arr::Array{Float64,1}, sel::Array{Int,1}, len::Int)
    outarr = zeros(Float64,len)
    i = 1
    for s in sel
        outarr[s] =  arr[i]
        i += 1
    end
    return outarr
end

function sumlessthan(x::Vector, maxsum::Float64=0.95)
    return x[12] + x[13] + x[14] - maxsum
end

function optimizemodel(dataset::Dataset, pairedsites::Array{Int,1}, siteCats::Int=3, lambdaCats::Int=5, fixGU::Bool=false, fixGCAU::Bool=false, integratestructure::Bool=true, unpairedmodel::Bool=false, initparams::ModelParameters=nothing, maxfile=nothing, usecuda::Bool=true; maxbasepairdistance::Int=1000000, maxoptiter::Int=2000, gapcutoff::Float64=1.0, coevolutionparamsonly::Bool=false)
    global optiter
    optiter = 0
    global maxll
    maxll = -1e20
    fixLambdaWeight = false
    paired = copy(pairedsites)
    if unpairedmodel
        fill!(paired, 0)
    end
    lower,upper,initialparams = getbounds(fixGU, fixLambdaWeight, unpairedmodel,siteCats,lambdaCats, getparamsvector(initparams), coevolutionparamsonly=coevolutionparamsonly)    
    for i=1:length(initialparams)
        initialparams[i] = min(upper[i], initialparams[i])
        initialparams[i] = max(lower[i], initialparams[i])
    end

    println("lower: ", lower)
    println("upper: ", upper)
    println("init: ", initialparams)
    println("unpairedmodel: ", unpairedmodel)

    usebayesopt = false
    minf = -Inf
    minx = nothing
    ret = ""
    if usebayesopt
        println("usebayesopt")
        f = (param -> modellikelihood(param, dataset, paired, siteCats, lambdaCats, fixGU, fixGCAU, integratestructure, unpairedmodel, maxfile, usecuda, maxbasepairdistance=maxbasepairdistance,gapcutoff=gapcutoff))

        model = ElasticGPE(15,
        mean = MeanConst(0.),         
        kernel = SEArd(zeros(Float64,15), 2.),
        capacity = 3000) 
        set_priors!(model.mean, [Normal(1, 2)])

        initializer_iterations = 100
        m = match(r"^([^\.]+)\..+$", basename(maxfile))
        cachefile = joinpath(dirname(maxfile),string("$(m[1]).cache.list",".sitecats$(siteCats).lambdaCats$(lambdaCats)"))
        if isfile(cachefile)
            x = Array{Float64,1}[]
            y = Float64[]
            fin = open(cachefile,"r")
            for line in readlines(fin)
                spl = split(line,"\t")
                if length(spl) == 2
                    push!(y, parse(Float64, spl[1]))
                    arrstr = split(spl[2][2:end-1], ", ")
                    push!(x, Float64[parse(Float64,v) for v in arrstr])
                end
            end
            close(fin)

            append!(model, hcat(x...), y)
            if length(x) > 5
                 initializer_iterations = 0
            end
        end

        # Optimize the hyperparameters of the GP using maximum a posteriori (MAP) estimates every 50 steps
        modeloptimizer = MAPGPOptimizer(every = 50,
        #kernbounds = [[-1, -1, 0], [4, 4, 10]],  # bounds of the 3 parameters GaussianProcesses.get_param_names(model.kernel)
        maxeval = 40)
        opt = BOpt(f,
        model,
        UpperConfidenceBound(),                # type of acquisition
        modeloptimizer,                        
        lower, upper,                  # lowerbounds, upperbounds         
        repetitions = 1,
        maxiterations = 500,                   # evaluate at 100 input positions
        sense = Max,
        verbosity = Progress,
        initializer_iterations = initializer_iterations)

        result = boptimize!(opt)
        println(result)
        exit()
    else
        
        rounds = 1
        if unpairedmodel
            rounds = 2
        end
        for roundnum=1:rounds
            localObjectiveFunction = ((param, grad) -> modellikelihood(param, dataset, paired, siteCats, lambdaCats, fixGU, fixGCAU, integratestructure, unpairedmodel, maxfile, usecuda, maxbasepairdistance=maxbasepairdistance, gapcutoff=gapcutoff))        
            opt = Opt(:LN_NELDERMEAD, 15)    
            if coevolutionparamsonly
                println("coevolutionparamsonly")
                opt = Opt(:GN_CRS2_LM, 15)
            end

            # unpaired optimization
            if rounds == 2 && roundnum == 1
                opt = Opt(:GN_CRS2_LM, 15)
            elseif rounds == 2 && roundnum == 2
                opt = Opt(:LN_NELDERMEAD, 15)
            end
            println(lower)
            println(upper)

            lower_bounds!(opt, lower)
            upper_bounds!(opt, upper)
            println("maxoptiter ", maxoptiter)
            maxeval!(opt, maxoptiter)
            max_objective!(opt, localObjectiveFunction)  
            xtol_rel!(opt,1e-4)
            

            (minf,minx,ret) = optimize(opt, initialparams)
            minf = max(minf,-1e20)
            println(minf,"\t",minx,"\t",ret)
            initialparams = minx

            if unpairedmodel
                if integratestructure
                    ll = getZ(minx, dataset, siteCats, lambdaCats, usecuda, true,gapcutoff=gapcutoff)
                    savemaximum(ll, minx, maxfile, true, string("Optimisation complete: ", ret))
                else
                    savemaximum(minf, minx, maxfile, true, string("Optimisation complete: ", ret))
                end
            end
        end
    end
    println("Complete")

    return getparams(minx, dataset, siteCats, lambdaCats,0,fixGU,fixGCAU)
end

#=
function optimizemodel(dataset::Dataset, pairedsites::Array{Int,1}, siteCats::Int=3, lambdaCats::Int=5, fixGU::Bool=false, fixGCAU::Bool=false, integratestructure::Bool=true, unpairedmodel::Bool=false,maxoptiter::Int=1000, initparams::ModelParameters=nothing, maxfile=nothing, usecuda::Bool=true; maxbasepairdistance::Int=1000000)
    global optiter
    optiter = 0
    global maxll
    maxll = -1e20
    fixLambdaWeight = false
    initialparams = 2.0*ones(Float64,20)
    if initialparams == nothing
        initialparams[1] = 1.5
        initialparams[2] = 1.5
        initialparams[3] = 1.5
        initialparams[12] = dataset.obsfreqs[1]
        initialparams[13] = dataset.obsfreqs[2]
        initialparams[14] = dataset.obsfreqs[3]
        initialparams[15] = 0.5
        initialparams[1] = 1.0
        initialparams[2] = 1.0
        initialparams[3] = 1.0
    else
        initialparams = getparamsvector(initparams)
    end

    paired = copy(pairedsites)
    if unpairedmodel
        fill!(paired, 0)
    end
    localObjectiveFunction = ((param, grad) -> modellikelihood(param, dataset, paired, siteCats, lambdaCats, fixGU, fixGCAU, integratestructure, unpairedmodel, maxfile, usecuda, maxbasepairdistance=maxbasepairdistance))
    opt = Opt(:LN_NELDERMEAD, 20)
    #opt = Opt(:LN_COBYLA, 20)

    #opt = Opt(:LN_NELDERMEAD, 15)
    #opt = Opt(:GN_ISRES, 20)
    lower = ones(Float64, 20)*0.0001
    lower[1] = 1.0
    lower[2] = 1.0
    lower[3] = 1.0
    #lower[10] = 0.0001
    lower[10] = 1.0 #GC
    lower[11] = 0.0001
    lower[12] = 0.01
    lower[13] = 0.01
    lower[14] = 0.01
    lower[15] = 0.0001
    if unpairedmodel
        lower[10] = 1.0
        lower[15] = 0.5
        #lower[11] = 1.0
    end
    lower[16] = 0.0001
    #lower[17] = 0.0001
    lower[17] = 1.0 # AT
    lower[18] = 0.0001
    lower[19] = 1.0 # GT
    #lower[19] = 0.0001
    lower[20] = 0.0001

    if fixLambdaWeight
        lower[15] = 0.5
    end
    lower_bounds!(opt, lower)

    upper = ones(Float64, 20)*50.0
    upper[1] = 1.0
    upper[2] = 1.0
    upper[3] = 1.0
    upper[12] = 2.0
    upper[13] = 2.0
    upper[14] = 2.0
    upper[15] = 0.9999
    upper[16] = 250.0
    if unpairedmodel
        upper[1] = 1.0
        upper[2] = 1.0
        upper[3] = 1.0
        upper[10] = 1.0
        upper[15] = 0.5
        #upper[11] = 1.0
    end

    if fixLambdaWeight
        upper[15] = 0.5
    end
    if fixGU
        upper[3] = 1.0
    end
    #ZZZZZZZZZZZZZZZZZZ
    upper[18] = 0.0001
    upper[20] = 0.0001
    upper_bounds!(opt, upper)

    xtol_rel!(opt,1e-4)
    maxeval!(opt, maxoptiter)

    max_objective!(opt, localObjectiveFunction)

    if unpairedmodel
        initialparams[1] = 1.0
        initialparams[2] = 1.0
        initialparams[3] = 1.0
        initialparams[10] = 1.0
    end
    if fixGU
        initialparams[3] = 1.0
    end
    if fixLambdaWeight
        initialparams[15] = 0.5
    end
    initialparams[15] = max(0.0001, initialparams[15])
    (minf,minx,ret) = optimize(opt, initialparams)
    println(minf,minx,ret)

    if unpairedmodel
        if integratestructure
            ll = getZ(minx, dataset, siteCats, lambdaCats, usecuda, true)
            savemaximum(ll, minx, maxfile, true, string("Optimisation complete: ", ret))
        else
            savemaximum(minf, minx, maxfile, true, string("Optimisation complete: ", ret))
        end
    end
    println("Complete")

    return getparams(minx, dataset, siteCats, lambdaCats,0,fixGU,fixGCAU)
end=#

function cleandataset(fastafile, outputfile)
    names = []
    sequences = []
    seqnametoindex  = Dict{AbstractString,Int}()
    FastaIO.FastaReader(fastafile) do fr
        seqindex = 1
        for (desc, seq) in fr
            len = length(seq)
            push!(names,desc)
            push!(sequences, seq)
            seqnametoindex[desc] = seqindex
            seqindex += 1
        end
    end

    #outfilename = string(outputdir, fastafile, ".norm")
    fout = open(outputfile, "w")
    seqindex = 1
    for seq in sequences
        write(fout, string(">seq", seqindex, "\n"))
        write(fout, string(uppercase(seq), "\n"))
        seqindex += 1
    end
    close(fout)
    return outputfile
end

function processmax(outputprefix, alignmentfile, maxfile::AbstractString, dataset::Dataset, params::ModelParameters, maxbasepairdistance::Int, usecuda::Bool, unpaired::Bool)
    siteloglikelihoods = nothing
    if 1 == 2
        if usecuda
            siteloglikelihoods = felsensteinposteriorsiterateconditionalscuda(unpairedposteriorprobs, pairedposteriorprobs,dataset, params, getmusespecificparamsarray(params))
        else
            siteloglikelihoods = getposteriorsiterates(unpairedposteriorprobs, pairedposteriorprobs, dataset, params, unpaired)
        end

        fout = open(string(maxfile, ".siterates.csv"),"w")
        pairingstr = ""
        if !run_args["unpaired"]
            pairingstr = string(",\"pairing probability\"")
        end
        println(fout, join(AbstractString[string("\"rate=", @sprintf("%.3f", siteRate), "\"") for siteRate in  currentparams.siteRates],","),",\"Mean\"", pairingstr)
        for i=1:dataset.numcols
            siteposteriorprobs = zeros(Float64, currentparams.siteCats)
            total = -Inf
            for siteCat=1:currentparams.siteCats
                siteposteriorprobs[siteCat] = siteloglikelihoods[siteCat, i]
                total = CommonUtils.logsumexp(total, siteposteriorprobs[siteCat])
            end
            siteposteriorprobs = exp.(siteposteriorprobs-total)
            pairingstr = ""
            if !run_args["unpaired"]
                pairingstr = string(",", 1.0 - unpairedposteriorprobs[i])
            end
            println(fout, join(AbstractString[string("\"", siteposteriorprobs[siteCat], "\"") for siteCat=1:currentparams.siteCats], ","),",", sum(siteposteriorprobs.*currentparams.siteRates), pairingstr)
        end
        close(fout)
    end

    if !unpaired
        unpairedlogprobs = computeunpairedlikelihoods(dataset, params)
        #pairedlogprobs, ret = coevolutionall(dataset,params,true,false,true,maxbasepairdistance,usecuda)
        println("G. maxbasepairdistance: ", maxbasepairdistance)
        pairedlogprobs, ret = coevolutionall(dataset,params,true,true,true,maxbasepairdistance,usecuda)
        inside = computeinsideKH99(unpairedlogprobs, pairedlogprobs, 1.0,false,usecuda)
        outside = computeoutsideKH99(inside,unpairedlogprobs, pairedlogprobs, usecuda)
        unpairedposteriorprobs,pairedposteriorprobs = computebasepairprobs(inside, outside, unpairedlogprobs, pairedlogprobs, KH99())
        fout = open(string(maxfile, ".posteriorunpaired"),"w")
        println(fout,unpairedposteriorprobs)
        close(fout)
        #=
        fout = open(string(maxfile, ".posteriorpaired"),"w")
        println(fout,printmatrix(pairedposteriorprobs))
        close(fout)
        =#

        ematrix = nothing
        smatrix = nothing
        consensus = getPosteriorDecodingConsensusStructure(pairedposteriorprobs, unpairedposteriorprobs, usecuda)
        fout = open(string(maxfile, ".consensus.dbn"),"w")
        println(fout, getdotbracketstring(consensus))
        close(fout)
        writectfile(consensus, repeat("N", length(consensus)), string(maxfile, ".consensus.ct"))

        #pairedlogprobs, ret = coevolutionall(dataset,params,true,true,true,maxbasepairdistance,usecuda)

        museparams = getmusespecificparamsarray(params)
        posteriorcoevolving = zeros(Float64, dataset.numcols, dataset.numcols)
        posteriormeanlambda = zeros(Float64, dataset.numcols, dataset.numcols)
        bayesfactorcoevolving = zeros(Float64, dataset.numcols, dataset.numcols)
        for i=1:dataset.numcols
            posteriorcoevolving[i,i] = 0.0
            for j=i+1:dataset.numcols
                posteriormeanlambda[i,j] = 0.0
                for k=1:length(museparams)
                    posteriormeanlambda[i,j] += museparams[k].lambdarate*exp(ret[k][i,j]-pairedlogprobs[i,j])
                    posteriormeanlambda[j,i] = posteriormeanlambda[i,j]
                    #println(k,"\t",ret[k][i,j],"\t",pairedlogprobs[i,j], "\t", exp(ret[k][i,j]-pairedlogprobs[i,j]))
                end
                postprob = exp(ret[1][i,j]-pairedlogprobs[i,j])
                bf = (postprob/(1.0-postprob))/(params.lambdazeroweight/(1.0-params.lambdazeroweight))
                bayesfactorcoevolving[i,j] = 1.0/bf
                bayesfactorcoevolving[j,i] = 1.0/bf
                posteriorcoevolving[i,j] = 1.0 - exp(ret[1][i,j]-pairedlogprobs[i,j])
                posteriorcoevolving[j,i] = posteriorcoevolving[i,j]
            end
        end
        fout = open(string(maxfile, ".consensus.csv"),"w")
        println(fout,"Position i,Paired with j,Unpaired probability,Pairing probability p(i^j),Post prob eta>0,Bayes fact eta>0,E[eta|i^j],E[eta],max(p[i,*])")
        for i=1:length(consensus)
            bp = "-,-,-,-,-,"
            if consensus[i] > 0
                #println(pairedposteriorprobs[i,consensus[i]])
                #println(posteriorcoevolving[i,consensus[i]])
                #println(bayesfactorcoevolving[i,consensus[i]])
                #println(posteriormeanlambda[i,consensus[i]])
                bp = string(@sprintf("%.4f", pairedposteriorprobs[i,consensus[i]]), ",", @sprintf("%.4f", posteriorcoevolving[i,consensus[i]]), ",", @sprintf("%.4e", bayesfactorcoevolving[i,consensus[i]]), ",", @sprintf("%.4f", posteriormeanlambda[i,consensus[i]]), ",", @sprintf("%.4f", pairedposteriorprobs[i,consensus[i]]*posteriormeanlambda[i,consensus[i]]))
            end
            println(fout, i,",", consensus[i], ",", @sprintf("%.4f", unpairedposteriorprobs[i]), ",", bp, ",", @sprintf("%.4f", maximum(pairedposteriorprobs[i,:])))
        end
        close(fout)
        lambdameans = zeros(Float64, length(consensus))
        postprobs = zeros(Float64, length(consensus))
        for i=1:length(consensus)
            if consensus[i] != 0
                lambdameans[i] = pairedposteriorprobs[i,consensus[i]]*posteriormeanlambda[i,consensus[i]]
                postprobs[i] = pairedposteriorprobs[i,consensus[i]]
            end
        end

        #println(lambdameans)
        try
            rankbycoevolution(outputprefix, alignmentfile, maxfile, consensus, Int[i for i=1:length(consensus)],dataset,lambdameans,postprobs)
        catch
            println("Error during ranking - structure is probably too short.")
        end

        cutoff = 0.001
        fout = open(string(maxfile,"_", cutoff, ".csv"),"w")
        println(fout, "\"i\",\"j\",\"paired(i)\",\"paired(j)\",\"paired(i,j)\",\"p(lambda>0)\",\"meanlambda\",\"bayesfactor(lambda>0)\"")
        for i=1:dataset.numcols
            for j=i+1:dataset.numcols
                if pairedposteriorprobs[i,j] >= cutoff
                    println(fout,"\"",i,"\",\"",j,"\",\"",1.0-unpairedposteriorprobs[i],"\",\"",1.0-unpairedposteriorprobs[j],"\",\"",pairedposteriorprobs[i,j],"\",\"",posteriorcoevolving[i,j],"\",\"",posteriormeanlambda[i,j],"\",\"",bayesfactorcoevolving[i,j],"\"")
                end
            end
        end
        close(fout)

        #=
        fout = open(string(maxfile, ".posteriorcoevolving"),"w")
        println(fout,printmatrix(posteriorcoevolving))
        close(fout)
        fout = open(string(maxfile, ".posteriormeanlambda"),"w")
        println(fout,printmatrix(posteriormeanlambda))
        close(fout)
        fout = open(string(maxfile, ".bayesfactorcoevolving"),"w")
        println(fout,printmatrix(bayesfactorcoevolving))
        close(fout)=#
    end
end

function processmax(outputprefix, alignmentfile, maxfile::AbstractString, dataset::Dataset, params::ModelParameters, paired::Array{Int,1}, mapping, usecuda::Bool)
    museparams = getmusespecificparamsarray(params)
    unpaired,pairedloglikelihoods,sampledstates,museconditionals = computeuncachedlikelihood(dataset, params, museparams, zeros(Int,1), paired, 1.0, false)
    lambdameans = zeros(Float64, length(paired))
    postprobs = zeros(Float64, length(paired))
    for i=1:length(paired)
        if paired[i] > i
            totalll = -Inf
            for musespecificparam in museparams
                lambdacat = musespecificparam.lambdacat
                totalll = CommonUtils.logsumexp(totalll, museconditionals[i,lambdacat])
                #
            end
            lambdamean = 0.0
            for musespecificparam in museparams
                lambdacat = musespecificparam.lambdacat
                if lambdacat == 1
                    postprobs[i] = 1.0-exp(museconditionals[i,1]-totalll)
                end
                lambdamean += exp.(museconditionals[i,lambdacat]-totalll)*musespecificparam.lambdarate
            end
            lambdameans[i] = lambdamean
            lambdameans[paired[i]] = lambdamean
        end
    end
    #println("Lambdas: ",lambdameans)
    fout = open(string(maxfile, ".mapped.dbn"),"w")
    println(fout, getdotbracketstring(paired))
    close(fout)
    rankbycoevolution(outputprefix, alignmentfile, maxfile, paired, mapping,dataset,lambdameans, postprobs)
    return lambdameans,postprobs
end

function computemaxstructure(dataset::Dataset, params::ModelParameters, maxbasepairdistance::Int=500, usecuda::Bool=true, deletegaps::Bool=false; gapcutoff::Float64=0.5)
    unpairedlogprobs = computeunpairedlikelihoods(dataset, params)
    println("A. maxbasepairdistance: ", maxbasepairdistance)
    pairedlogprobs, ret = coevolutionall(dataset,params,true,true,true,maxbasepairdistance,usecuda)

    if deletegaps
        mapping = zeros(Int,dataset.numcols)
        revmapping = zeros(Int, dataset.numcols)
        cutoff = 0.5
        index = 1
        for i=1:dataset.numcols
            if dataset.gapfrequency[i] < cutoff
                mapping[i] = index
                revmapping[index] = i
                index += 1
            end
        end
        newlen = index - 1

        unpairedlogprobs_trunc = ones(Float64, newlen)*-Inf
        pairedlogprobs_trunc = ones(Float64, newlen, newlen)*-Inf
        for i=1:dataset.numcols
            if mapping[i] > 0
                unpairedlogprobs_trunc[mapping[i]] = unpairedlogprobs[i]
                for j=1:dataset.numcols
                    if mapping[j] > 0
                        pairedlogprobs_trunc[mapping[i],mapping[j]] = pairedlogprobs[i,j]
                    end
                end
            end
        end
        inside = computeinsideKH99(unpairedlogprobs_trunc, pairedlogprobs_trunc, 1.0,false,usecuda)
        outside = computeoutsideKH99(inside,unpairedlogprobs_trunc, pairedlogprobs_trunc, usecuda)
        unpairedposteriorprobs_trunc,pairedposteriorprobs_trunc = computebasepairprobs(inside, outside, unpairedlogprobs_trunc, pairedlogprobs_trunc, KH99())
        paired_trunc = getPosteriorDecodingConsensusStructure(pairedposteriorprobs_trunc, unpairedposteriorprobs_trunc, usecuda)

        paired = zeros(Int, dataset.numcols)
        for i=1:dataset.numcols
            if mapping[i] > 0
                if paired_trunc[mapping[i]] > 0
                    paired[i] = revmapping[paired_trunc[mapping[i]]]
                end
            end
        end

        return paired
    else
        maskgapped!(pairedlogprobs,dataset.gapfrequency,gapcutoff,-Inf)
        inside = computeinsideKH99(unpairedlogprobs, pairedlogprobs, 1.0,false,usecuda)
        outside = computeoutsideKH99(inside,unpairedlogprobs, pairedlogprobs, usecuda)
        unpairedposteriorprobs,pairedposteriorprobs = computebasepairprobs(inside, outside, unpairedlogprobs, pairedlogprobs, KH99())
        return getPosteriorDecodingConsensusStructure(pairedposteriorprobs, unpairedposteriorprobs, usecuda)
    end
end

function likelihoodratiotests(outputprefix)
    try
        outputdir = dirname(outputprefix)
        csvfile = abspath(joinpath(outputdir,"likelihoodratiotests.csv"))
        fout = open(csvfile, "w")
        write(fout, printaictablescsv(String[outputprefix]))
        close(fout)

        template = read(open(joinpath(@__DIR__,"likelihoodratiotests_template.tex"),"r"), String)

        texfile = abspath(joinpath(outputdir,"likelihoodratiotests.tex"))
        fout = open(texfile,"w")
        write(fout, replace(template, "#INSERT#" => printaictables(String[outputprefix])))
        close(fout)
        run(`pdflatex -output-directory="$(abspath(joinpath(outputdir)))" $(texfile)`)
        rm(joinpath(outputdir,"likelihoodratiotests.aux"))
        rm(joinpath(outputdir,"likelihoodratiotests.log"))
    catch

    end
end

function main()
    parsed_args = parse_commandline()
    run_args_list = [deepcopy(parsed_args)]

    if parsed_args["full"]
        run_args_list = []

        args = deepcopy(parsed_args)
        args["unpaired"] = true
        args["fixgcau"] = false
        args["fixgu"] = false
        push!(run_args_list, args)

        args = deepcopy(parsed_args)
        args["unpaired"] = false
        args["fixgcau"] = true
        args["fixgu"] = true
        push!(run_args_list, args)

        args = deepcopy(parsed_args)
        args["unpaired"] = false
        args["fixgcau"] = true
        args["fixgu"] = false
        push!(run_args_list, args)

        args = deepcopy(parsed_args)
        args["unpaired"] = false
        args["fixgcau"] = false
        args["fixgu"] = true
        push!(run_args_list, args)

        args = deepcopy(parsed_args)
        args["unpaired"] = false
        args["fixgcau"] = false
        args["fixgu"] = false
        push!(run_args_list, args)
    end

    seed1 = parsed_args["seed"]
    rng = MersenneTwister(seed1)
    seed2 = rand(rng,1:typemax(Int64))
    Random.seed!(seed2)

    usecuda = !parsed_args["cpuonly"]

    alignmentfilein = parsed_args["alignment"]
    defaultname = split(basename(parsed_args["alignment"]), ".")[1]
    outputprefix = abspath(joinpath("results", defaultname, defaultname))
    if parsed_args["outputprefix"] != nothing
        outputprefix = abspath(string(parsed_args["outputprefix"]))
    end
    println("Output will be written to: \"", outputprefix, "\"")
    outputprefixold = string(outputprefix)
    if parsed_args["shuffle"]
        outputprefix = string(outputprefix, ".shuffle",seed1)
    end
    outputdir = dirname(outputprefix)
    mkpath(outputdir)
    newlog(string(outputprefix, ".log"))
    if parsed_args["shuffle"]
        alignmentfilein = shufflealignment(rng, alignmentfilein, string(outputprefix, ".fas"))
    end
    alignmentfile = cleandataset(alignmentfilein, string(outputprefix,".fas.norm"))

    if parsed_args["tree"] == nothing
        treefile = "$outputprefix.nwk"
        newickstring, treepath = Binaries.fasttreegtr(alignmentfile)
        fout = open(treefile,"w")
        println(fout, newickstring)
        close(fout)
    else
        treefile = parsed_args["tree"]
    end
    grammar = KH99()
    dataset = Dataset(getalignmentfromfastaandnewick(alignmentfile,treefile))
    
    #=
    println("Numcols = ", dataset.numcols)
    countbasepairs = 0
    totalbasepairs = 0
    for x=1:dataset.numcols
        println(x)
        for y=x+1:dataset.numcols
            if hascanonicalbasepair(dataset, x, y)
                countbasepairs += 1
            end
            totalbasepairs += 1
        end
    end
    println(countbasepairs,"\t",totalbasepairs,"\t",countbasepairs/totalbasepairs)
    =#

    maxbasepairprobs = ones(Float64, dataset.numcols, dataset.numcols)
    if parsed_args["expcanonicalconstraints"] || parsed_args["expcanonicalrnafold"] != -1000.0        
        if parsed_args["expcanonicalrnafold"] != -1000.0
            foldtemp  = parsed_args["expcanonicalrnafold"]
            maxbasepairprobs = zeros(Float64, dataset.numcols, dataset.numcols)
            for seq in dataset.sequences
                mapping = zeros(Int, length(seq))
                z = 1
                for i=1:length(seq)
                    if seq[i] != '-'
                        mapping[i] = z
                        z += 1
                    end
                end
                revmapping = CommonUtils.getrevmapping(mapping)

                I,J,V = RNABinaries.rnafold_basepairprobs(seq, foldtemp, false, "/media/michael/Sandisk500GB/Dropbox/dev/phylo/src/rnabox/rnafoldcache.h5")
                for (a,b,prob) in zip(I,J,V)
                    x = revmapping[a]
                    y = revmapping[b]
                    maxbasepairprobs[x,y] = max(maxbasepairprobs[x,y], prob)
                    maxbasepairprobs[y,x] = maxbasepairprobs[x,y]
                end
            end
        end

        countbasepairs = 0
        totalbasepairs = 0
        for x=1:dataset.numcols
            for y=x+1:dataset.numcols
                if parsed_args["expcanonicalconstraints"] && !hascanonicalbasepair(dataset, x, y)
                    maxbasepairprobs[x,y] = 0.0
                    maxbasepairprobs[y,x] = 0.0
                elseif parsed_args["expcanonicalrnafold"] != -1000.0 && maxbasepairprobs[x,y] < 1e-5
                    maxbasepairprobs[x,y] = 0.0
                    maxbasepairprobs[y,x] = 0.0
                else
                    countbasepairs += 1
                end
                totalbasepairs += 1
            end
        end        
        println(countbasepairs,"\t",totalbasepairs,"\t",countbasepairs/totalbasepairs)
    end
    dataset.maxbasepairprobs = sparse(maxbasepairprobs)

    for run_args in run_args_list
        mcmclogfile = string(outputprefix,".log")

        mode = MODE_IO_SAMPLE
        #mode = MODE_VIENNA_SAMPLE
        if run_args["structure"] != nothing
            sequence, pairedsites, mapping, revmapping = mapstructure(dataset,alignmentfile,run_args["structure"])
            mode = MODE_FIXED_STRUCTURE
            println("MODE_FIXED_STRUCTURE")
        end

        M = run_args["M"]

        lambdaCats = run_args["numlambdacats"]
        lambdarates = ones(Float64,lambdaCats)
        lambdaweights = ones(Float64,length(lambdarates)) / length(lambdarates)

        samplebranchlengths = false
        currentparams = ModelParameters(dataset.obsfreqs, 8.0, 4.0, 2.0, 1.0, 5.0, 1.0, 1.0, 5.0, 0.25, getnodelist(dataset.root), 0.5, lambdarates, lambdaweights)
        currentparams.lambdaCats = lambdaCats - 1
        #currentparams.lambdaweightsGC, currentparams.lambdaratesGC = discretizegamma2(currentparams.lambdazeroweight, currentparams.lambdaGammaShapeGC, currentparams.lambdaGammaScaleGC, currentparams.lambdaCats)
        #currentparams.lambdaweightsAT, currentparams.lambdaratesAT = discretizegamma2(currentparams.lambdazeroweight, currentparams.lambdaGammaShapeAT, currentparams.lambdaGammaScaleAT, currentparams.lambdaCats)
        #currentparams.lambdaweightsGT, currentparams.lambdaratesGT = discretizegamma2(currentparams.lambdazeroweight, currentparams.lambdaGammaShapeGT, currentparams.lambdaGammaScaleGT, currentparams.lambdaCats)
        currentparams.lambdaweightsGC, currentparams.lambdaratesGC, currentparams.lambdaweightsAT, currentparams.lambdaratesAT, currentparams.lambdaweightsGT, currentparams.lambdaratesGT = discretizegamma3(currentparams.lambdazeroweight, currentparams.lambdaGammaShapeGC, currentparams.lambdaGammaShapeAT, currentparams.lambdaGammaShapeGT, currentparams.lambdaGammaScaleGC, currentparams.lambdaCats)
        siteCats = run_args["numsitecats"]
        currentparams.siteCats = siteCats
        currentparams.siteWeights = ones(Float64,currentparams.siteCats)/currentparams.siteCats
        currentparams.siteRates = discretizegamma(currentparams.siteGammaShape, currentparams.siteGammaScale, currentparams.siteCats)
        currentparams.states = Int[rand(rng,1:currentparams.siteCats) for i=1:dataset.numcols]
        currentparams.pairedstates = Int[rand(rng,1:lambdaCats) for i=1:dataset.numcols]


        initialparams = Float64[2.0,2.0,2.0,1.7487740904105369,4.464074402974858,1.6941505179847807,0.4833108030758708,5.839004646491171,0.7168678100059017,1.0,0.6118067582467858,0.23307618715645315,0.2631203272837885,0.2430685428508905,0.5]
        initialparams[10] = 1.0


        maxbasepairdistance = run_args["maxbasepairdistance"]

        integratesiterates = true
        integratestructure = true

        maxfile = string(outputprefix, ".max",".sitecats$(siteCats).lambdaCats$(lambdaCats)")
        maxfilefixgu = string(outputprefix, ".fixgu.max",".sitecats$(siteCats).lambdaCats$(lambdaCats)")
        maxfileold = string(outputprefixold,".max",".sitecats$(siteCats).lambdaCats$(lambdaCats)")
        maxfileunpaired = string(outputprefix, ".max.unpaired",".sitecats$(siteCats)")
        if isfile(maxfile)
            jsondict = parsejsonfile(maxfile)
            initialparams = convert(Array{Float64,1}, jsondict["maxparams"])
        elseif isfile(maxfileold)
            jsondict = parsejsonfile(maxfileold)
            initialparams = convert(Array{Float64,1}, jsondict["maxparams"])
        elseif isfile(maxfilefixgu)
            jsondict = parsejsonfile(maxfilefixgu)
            initialparams = convert(Array{Float64,1}, jsondict["maxparams"])
            initialparams[3] = min(1.25, (initialparams[2]-1.0)*0.5 + 1.0)
        elseif isfile(maxfileunpaired)
            jsondict = parsejsonfile(maxfileunpaired)
            initialparams = convert(Array{Float64,1}, jsondict["maxparams"])
            initialparams[1] = 5.0
            initialparams[2] = 4.0
            initialparams[3] = 3.0
            initialparams[15] = 0.25
        end
        currentparams = getparams(initialparams,dataset,siteCats,lambdaCats,0)
      
       if run_args["fixgcau"] && run_args["fixgu"]
            maxfile = string(outputprefix, ".fixgcaugu.max",".sitecats$(siteCats).lambdaCats$(lambdaCats)")
            if isfile(maxfile)
                jsondict = parsejsonfile(maxfile)
                initialparams = convert(Array{Float64,1}, jsondict["maxparams"])
            end
            currentparams = getparams(initialparams,dataset,siteCats,lambdaCats,0)
            #currentparams.lambdazeroweight = max(currentparams.lambdazeroweight, 0.05)
        elseif run_args["fixgcau"]
            maxfile = string(outputprefix, ".fixgcau.max",".sitecats$(siteCats).lambdaCats$(lambdaCats)")
            if isfile(maxfile)
                jsondict = parsejsonfile(maxfile)
                initialparams = convert(Array{Float64,1}, jsondict["maxparams"])
            end
            currentparams = getparams(initialparams,dataset,siteCats,lambdaCats,0)
            #currentparams.lambdazeroweight = max(currentparams.lambdazeroweight, 0.05)
        elseif run_args["fixgu"]
            maxfile = string(outputprefix, ".fixgu.max",".sitecats$(siteCats).lambdaCats$(lambdaCats)")
            if isfile(maxfile)
                jsondict = parsejsonfile(maxfile)
                initialparams = convert(Array{Float64,1}, jsondict["maxparams"])
            end
            currentparams = getparams(initialparams,dataset,siteCats,lambdaCats,0)
            #currentparams.lambdazeroweight = max(currentparams.lambdazeroweight, 0.05)
        end

        if run_args["unpaired"]
            maxfile = string(outputprefix, ".max.unpaired",".sitecats$(siteCats)")
            if isfile(maxfile)
                jsondict = parsejsonfile(maxfile)
                initialparams = convert(Array{Float64,1}, jsondict["maxparams"])
            end
            currentparams = getparams(initialparams,dataset,siteCats,lambdaCats,0)
        end

        if run_args["processmax"]
            if run_args["structure"] != nothing
                processmax(outputprefix, alignmentfile, maxfile, dataset, currentparams, pairedsites, mapping, usecuda)
            else
                processmax(outputprefix, alignmentfile, maxfile, dataset, currentparams, maxbasepairdistance, usecuda, run_args["unpaired"])
            end
            exit()
        end

        if run_args["maxstructure"]
            println("Computing maximum likelihood structure")
            maxstructure_trunc = computemaxstructure(dataset, currentparams, maxbasepairdistance, usecuda,true, gapcutoff=run_args["gapcutoff"])
            writectfile(maxstructure_trunc, replace(dataset.sequences[1], "-" => "N"), string(outputprefix, ".maxstructuretrunc"))
            maxstructure = computemaxstructure(dataset, currentparams, maxbasepairdistance, usecuda, gapcutoff=run_args["gapcutoff"])
            writectfile(maxstructure, replace(dataset.sequences[1], "-" => "N"), string(outputprefix, ".maxstructure"))
            exit()
        end

        if run_args["simulate"] > 0
            simulatealignments(rng, dataset, currentparams, joinpath(outputdir,"simulations/"); samples=run_args["simulate"], usecuda=usecuda, maxbasepairdistance=maxbasepairdistance,fixGU=run_args["fixgu"],fixGCAU=run_args["fixgcau"],unpairedmodel=run_args["unpaired"])
            maxfiles = [f for f in filter(x -> match(r".+(.max(\..*)*)$", x) != nothing, readdir(outputdir))]
            
            for sim=1:run_args["simulate"]
                dst = joinpath(outputdir,"simulations/sim$(sim).nwk")
                if !isfile(dst)
                    cp(treefile,dst)
                end
            end
            for maxfile in maxfiles
                if !parsed_args["fixgu"] || (occursin(".fixgu.max", maxfile) || occursin(".max.unpaired",maxfile))
                    for sim=1:run_args["simulate"]
                        extension = match(r"(\..*)",maxfile)[1]
                        src = joinpath(outputdir,maxfile)
                        jsondict = JSON.parse(read(src, String))
                        jsondict["Z"] = -1.0e20
                        dst = joinpath(outputdir,"simulations/sim$(sim)$(extension)")
                        #println(dst)
                        if !isfile(dst)
                            open(dst,"w") do f
                                JSON.print(f, jsondict)
                            end
                        end
                    end
                end
            end
            exit()
        end

        if run_args["calcentropy"] || run_args["calcpriorentropy"]
            initialparams = convert(Array{Float64,1}, jsondict["maxparams"])
            currentparams = getparams(initialparams,dataset,siteCats,lambdaCats,0)
            println("PARAMS",initialparams)
            h, hmax = computeinformationentropy(dataset, currentparams, maxbasepairdistance, usecuda)
            #println("H=", h)
            H, Hstdev, Hstderr, Hmax, perc = estimateinformationentropy(rng, dataset, currentparams, maxbasepairdistance, usecuda, run_args["calcpriorentropy"])
            jsondict2 = Dict()
            jsondict2["params"] = initialparams
            jsondict2["H"] = H
            jsondict2["Hstdev"] = Hstdev
            jsondict2["Hstderr"] = Hstderr
            jsondict2["Hmax"] = Hmax
            jsondict2["percentage"] = perc
            jsondict2["length"] = dataset.numcols
            entropyfile =  string(outputprefix, ".entropy")
            if run_args["fixgu"]
                entropyfile = string(outputprefix, ".fixgu.entropy")
            end
            if run_args["calcpriorentropy"]
                entropyfile =  string(outputprefix, ".entropyprior")
            end
            out = open(entropyfile,"w")
            ret = replace(JSON.json(jsondict2),",\"" => ",\n\"")
            ret = replace(ret, "],[" => "],\n[")
            ret = replace(ret, "{" => "{\n")
            ret = replace(ret, "}" => "\n}")
            write(out,ret)
            close(out)

            println("Entropy = ", H)
            #println("Entropy est ", hest)
            println("Max. entropy = ", Hmax)
            println("Norm. entropy = ", perc)
            exit()
        end

        optimize = !run_args["mcmc"]
        if optimize
            fixGU = run_args["fixgu"]
            fixGCAU = run_args["fixgcau"]
            integratestructure = run_args["structure"] == nothing
            unpairedmodel = run_args["unpaired"]
            if integratestructure && unpairedmodel
                integratestructure = false
            end

            if integratestructure
                samplesperiter = run_args["samplesperiter"]
                if samplesperiter > 0
                    println("Using Monte Carlo EM")
                    currentparams = mcem(dataset, siteCats, lambdaCats, fixGU, fixGCAU, unpairedmodel,run_args["maxmcemiter"],samplesperiter,maxbasepairdistance,currentparams, maxfile, usecuda)
                else
                    currentparams = optimizemodel(dataset, (run_args["structure"] == nothing ? zeros(Int,dataset.numcols) : pairedsites), siteCats, lambdaCats, fixGU,fixGCAU, true, unpairedmodel,currentparams, maxfile, usecuda, maxbasepairdistance=maxbasepairdistance, gapcutoff=run_args["gapcutoff"], maxoptiter=run_args["maxoptiter"], coevolutionparamsonly=run_args["optimizecoevolutiononly"])
                end
            else
                currentparams = optimizemodel(dataset,  (run_args["structure"] == nothing ? zeros(Int,dataset.numcols) : pairedsites), siteCats,lambdaCats,fixGU,fixGCAU,integratestructure,unpairedmodel,currentparams, maxfile, usecuda, maxbasepairdistance=maxbasepairdistance, gapcutoff=run_args["gapcutoff"], maxoptiter=run_args["maxoptiter"], coevolutionparamsonly=run_args["optimizecoevolutiononly"])
            end
        end

        if run_args["mcmc"]
            #=
            parameterisation = 0
            if parameterisation == 1
                gammameanGC = currentparams.lambdaGammaShapeGC*currentparams.lambdaGammaScaleGC
                #currentparams.lambdaGammaScaleGC = gammameanGC/currentparams.lambdaGammaScaleGC
                currentparams.lambdaGammaShapeGC = gammameanGC

                gammameanAT = currentparams.lambdaGammaShapeAT*currentparams.lambdaGammaScaleGC
                #currentparams.lambdaGammaScaleAT = gammameanAT/currentparams.lambdaGammaScaleAT
                currentparams.lambdaGammaShapeAT = gammameanAT

                gammameanGT = currentparams.lambdaGammaShapeGT*currentparams.lambdaGammaScaleGC
                #currentparams.lambdaGammaScaleGT = gammameanGT/currentparams.lambdaGammaScaleGT
                currentparams.lambdaGammaShapeGT = gammameanGT
            end
            currentparams.lambdaweightsGC, currentparams.lambdaratesGC, currentparams.lambdaweightsAT, currentparams.lambdaratesAT, currentparams.lambdaweightsGT, currentparams.lambdaratesGT = discretizegamma3(currentparams.lambdazeroweight, currentparams.lambdaGammaShapeGC, currentparams.lambdaGammaShapeAT, currentparams.lambdaGammaShapeGT, currentparams.lambdaGammaScaleGC, currentparams.lambdaCats, parameterisation)
            =#
            proposedparams = deepcopy(currentparams)


            inside = zeros(Float64,1,1,1)
            if mode == MODE_IO_SAMPLE
                fout = open(string(outputprefix, ".calculations"), "w")
                computations = countcomputations(dataset, maxbasepairdistance)
                write(fout, string(computations, "\n"))
                write(fout, string(dataset.numseqs,"\t", dataset.numcols,"\n"))
                write(fout, string(@sprintf("%.3f", computations[6]), "\n"))
                write(fout, string(@sprintf("%.2f", 1.0/computations[6]), "\n"))
                close(fout)

                unpairedlogprobs = computeunpairedlikelihoods(dataset, currentparams)
                museparams = getmusespecificparamsarray(currentparams)
                #tic()
                println("B. maxbasepairdistance: ", maxbasepairdistance)
                pairedlogprobs, ret = coevolutionall(dataset,currentparams,true,false,true,maxbasepairdistance,usecuda)
                maskgapped!(pairedlogprobs,dataset.gapfrequency,run_args["gapcutoff"],-Inf)
                inside = computeinsideKH99(unpairedlogprobs, pairedlogprobs, 1.0,false,usecuda)
                savemaximum(inside[1,1,dataset.numcols], getparamsvector(currentparams), maxfile)
                #elapsed = toc();
                #println("MCMC initialised: ", elapsed)
            end

            thermodynamicsamples = Dict{Int, Array{Array{Int,1}}}()

            #Bs = [1.0]
            Bs = [run_args["bfactor"]^(z-1.0) for z=1:run_args["numchains"]]

            if samplebranchlengths
                currentparams.q6 = 1.0
            end
            chains = Chain[]
            id = 1
            for B in Bs
                if mode == MODE_IO_SAMPLE
                    chain = Chain(id, B, MersenneTwister(rand(rng,1:typemax(Int64))), 0.0, 0.0, currentparams, zeros(Int,dataset.numcols), copy(inside), pairedlogprobs, unpairedlogprobs, copy(inside), copy(pairedlogprobs), copy(unpairedlogprobs), deepcopy(currentparams))
                    chain.pairedlogprior = samplestructure(chain.rng, chain.inside, chain.pairedlogprobs, chain.unpairedlogprobs, 1, dataset.numcols, chain.paired, grammar, chain.B)
                elseif mode == MODE_VIENNA_SAMPLE
                    chain = Chain(id, B, MersenneTwister(rand(rng,1:typemax(Int64))), 0.0, 0.0, currentparams, zeros(Int,dataset.numcols), zeros(Float64,1,1,1), zeros(Float64,1,1), zeros(Float64,1))
                    chain.paired = samplethermodynamic(thermodynamicsamples, rng,dataset.sequences)
                else
                    chain = Chain(id, B, MersenneTwister(rand(rng,1:typemax(Int64))), 0.0, 0.0, currentparams, zeros(Int,dataset.numcols), zeros(Float64,1,1,1), zeros(Float64,1,1), zeros(Float64,1), zeros(Float64,1,1,1), zeros(Float64,1,1), zeros(Float64,1), deepcopy(currentparams))
                    chain.paired = copy(pairedsites)
                end
                chain.currentll = computetotallikelihood(chain.rng, dataset, chain.currentparams, chain.paired, samplebranchlengths,integratesiterates,integratestructure, M)
                chain.proposedll = chain.currentll
                push!(chains, chain)
                id += 1
            end


            burnins = Int[250, 500, 1000, 1500, 2000,2500, 3000, 4000,5000,6000, 7500, 9000, 10000]
            if samplebranchlengths
                burnins = Int[max(1500,dataset.numseqs*100),max(4500,dataset.numseqs*300)]
            end
            if length(burnins) > 0
                burnindata = zeros(Float64,burnins[end],10)
            end



            println("integratesiterates ", integratesiterates)
            println("integratestructure ", integratestructure)
            mcmcoptions = MCMCoptions(mode,M,samplebranchlengths,integratesiterates,integratestructure,maxbasepairdistance,usecuda)

            numChains =  length(Bs)
            swapacceptance = ones(Float64, numChains, numChains)*5e-11
            swaptotal = ones(Float64, numChains, numChains)*1e-10
            maxll = -Inf
            tuningvectors = Array{Float64,1}[]
            branchtuningvectors = Array{Float64,1}[]
            covmat = nothing
            for i=1:300
                for chain in chains
                    covmat = runchain(10,chain, dataset, grammar, outputprefix, string(outputprefix,".B",chain.B,".M",M), mcmcoptions, burnins, burnindata, thermodynamicsamples, tuningvectors, branchtuningvectors, covmat, siteCats=siteCats, lambdaCats=lambdaCats, fixGU=false, fixGCAU=false, gapcutoff=run_args["gapcutoff"])
                    if chain.currentll > maxll
                        maxll = chain.currentll
                    end
                end
                for k=1:length(chains)
                    sel = [i for i=1:numChains]
                    shuffle!(rng, sel)
                    if length(sel) > 1
                        a = sel[1]
                        b = sel[2]
                        S = (chains[b].B-chains[a].B)*(chains[a].currentll + chains[a].pairedlogprior - chains[b].currentll - chains[b].pairedlogprior)
                        if exp(S) > rand(rng)
                            chains[a].B, chains[b].B = chains[b].B, chains[a].B
                            chains[a].id, chains[b].id = chains[b].id, chains[a].id
                            chains[a].logger, chains[b].logger = chains[b].logger, chains[a].logger
                            chains[a].timings, chains[b].timings = chains[b].timings, chains[a].timings
                            swapacceptance[chains[a].id,chains[b].id] += 1.0
                            swapacceptance[chains[b].id,chains[a].id] += 1.0
                        end
                        swaptotal[chains[a].id,chains[b].id] += 1.0
                        swaptotal[chains[b].id,chains[a].id] += 1.0
                        println(swapacceptance./swaptotal)
                    end
                end

            end
        end
    end


    if parsed_args["full"]
        likelihoodratiotests(outputprefix)
    end
end

main()

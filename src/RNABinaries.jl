module RNABinaries
    using SHA
    using HDF5
    using Blosc
    using FastaIO


    push!(LOAD_PATH,joinpath(@__DIR__))
    #using RNATools
    using CommonUtils


    rnafold_binary = "RNAfold"
    rnafold_version = "Linux 2.3.5"
    if Sys.iswindows()
        rnafold_binary = "C:\\Program Files (x86)\\ViennaRNA Package\\RNAfold.exe"
        rnafold_version = "Windows 64bit 2.3.5"
    end

    function isrnafoldcached(inseq::AbstractString, temperature::Float64=37.0, circular::Bool=false, h5file::AbstractString="rnafoldcache.h5")
        seq = CommonUtils.canonicalise_sequence(inseq)
        key = CommonUtils.sha256base36(seq)
        tempstr = string(temperature)
        index = 1
        path = "$(key)/$(tempstr)/$(index)"
        iscached = h5open(h5file, "r") do file
            exists(file, path)
        end
        return iscached
    end

    function rnafold(inseq::AbstractString, temperature::Float64=37.0, circular::Bool=false, h5file::AbstractString="rnafoldcache.h5")
        seq = CommonUtils.canonicalise_sequence(inseq)
        key = CommonUtils.sha256base36(seq)
        tempstr = string(temperature)
        index = 1
        path = "$(key)/$(tempstr)/$(index)"
        iscached = h5open(h5file, "cw") do file
            exists(file, path)
        end

        if !iscached
            datalen = length(seq)
            tempfile, fout = mktemp()
            println(fout, ">seq")
            println(fout, seq)
            close(fout)
            iobuf = IOBuffer()
            output = ""
            conf = ""
            if circular
                output = readstring(pipeline(tempfile, `$(rnafold_binary) --temp=$(temperature) --partfunc --circ --id-prefix="$key"`))
            else
                output = readstring(pipeline(tempfile, `$(rnafold_binary) --temp=$(temperature) --partfunc --id-prefix="$key"`))
            end

            rm(tempfile)
            lines = split(output, r"[\r\n]+")
            fname = strip(lines[1])[2:end]
            fin = open("$(fname)_dp.ps", "r")
            I = Int[]
            J = Int[]
            V = Float64[]
            for ln in readlines(fin)
                line = strip(ln)
                if endswith(line, "ubox") && !startswith(line, "%")
                    spl = split(line)
                    push!(I, parse(Int, spl[1]))
                    push!(J, parse(Int, spl[2]))
                    push!(V, parse(Float64, spl[3])^2.0)
                end
            end
            probmatsparse = sparse(I,J,V,datalen,datalen)
            probarr = sum(Symmetric(probmatsparse),1)
            close(fin)
            rm("$(fname)_dp.ps")
            rm("$(fname)_ss.ps")
            splline3 = split(lines[3], r"\s+", limit=2)
            mfedbn = splline3[1]
            mfepaired = RNATools.parsedbn(mfedbn)
            free_energy = parse(Float64, strip(strip(splline3[2])[2:end-1]))
            splline6 = split(lines[6], ";", limit=2)
            ensembleprob = parse(Float64, split(strip(splline6[1]))[end])
            ensemblediversity = parse(Float64, split(strip(splline6[2]))[end])

            ret = h5open(h5file, "cw") do file
                while exists(file, path)
                    index += 1
                    path = "$(key)/$(tempstr)/$(index)"
                end
                h5write(h5file, "$(path)/sequence", seq, "blosc", 3)
                h5write(h5file, "$(path)/mfepaired", mfepaired, "blosc", 3)
                h5write(h5file, "$(path)/freeenergy", free_energy)
                h5write(h5file, "$(path)/ensembleprob", ensembleprob)
                h5write(h5file, "$(path)/ensemblediversity", ensemblediversity)
                h5write(h5file, "$(path)/probarr", probarr[:], "blosc", 3)
                h5write(h5file, "$(path)/I", I, "blosc", 3)
                h5write(h5file, "$(path)/J", J, "blosc", 3)
                h5write(h5file, "$(path)/V", V, "blosc", 3)
                h5write(h5file, "$(path)/version", rnafold_version)
                if circular
                    h5write(h5file, "$(path)/conformation", "circular")
                else
                    h5write(h5file, "$(path)/conformation", "linear")
                end
            end
        end
        return path
    end


    function rnafold2(inseq::AbstractString, temperature::Float64=37.0, circular::Bool=false, h5file::AbstractString="rnafoldcache.h5")
        seq = CommonUtils.canonicalise_sequence(inseq)
        key = CommonUtils.sha256base36(seq)
        tempstr = string(temperature)
        index = 1
        path = "$(key)/$(tempstr)/$(index)"
        iscached = h5open(h5file, "cw") do file
            exists(file, path)
        end

        if !iscached
            datalen = length(seq)
            tempfile, fout = mktemp()
            println(fout, ">seq")
            println(fout, seq)
            close(fout)
            iobuf = IOBuffer()
            output = ""
            conf = ""
            if circular
                output = readstring(pipeline(tempfile, `$(rnafold_binary) --temp=$(temperature) --partfunc --circ --id-prefix="$key"`))
            else
                output = readstring(pipeline(tempfile, `$(rnafold_binary) --temp=$(temperature) --partfunc --id-prefix="$key"`))
            end

            rm(tempfile)
            lines = split(output, r"[\r\n]+")
            fname = strip(lines[1])[2:end]
            fin = open("$(fname)_dp.ps", "r")
            I = Int[]
            J = Int[]
            V = Float64[]
            for ln in readlines(fin)
                line = strip(ln)
                if endswith(line, "ubox") && !startswith(line, "%")
                    spl = split(line)
                    push!(I, parse(Int, spl[1]))
                    push!(J, parse(Int, spl[2]))
                    push!(V, parse(Float64, spl[3])^2.0)
                end
            end
            probmatsparse = sparse(I,J,V,datalen,datalen)
            probarr = sum(Symmetric(probmatsparse),1)
            close(fin)
            rm("$(fname)_dp.ps")
            rm("$(fname)_ss.ps")
            splline3 = split(lines[3], r"\s+", limit=2)
            mfedbn = splline3[1]
            mfepaired = RNATools.parsedbn(mfedbn)
            free_energy = parse(Float64, strip(strip(splline3[2])[2:end-1]))
            splline6 = split(lines[6], ";", limit=2)
            ensembleprob = parse(Float64, split(strip(splline6[1]))[end])
            ensemblediversity = parse(Float64, split(strip(splline6[2]))[end])

            ret = h5open(h5file, "cw") do file
                while exists(file, path)
                    index += 1
                    path = "$(key)/$(tempstr)/$(index)"
                end
                h5write(h5file, "$(path)/sequence", seq, "blosc", 3)
                h5write(h5file, "$(path)/mfepaired", mfepaired, "blosc", 3)
                h5write(h5file, "$(path)/freeenergy", free_energy)
                h5write(h5file, "$(path)/ensembleprob", ensembleprob)
                h5write(h5file, "$(path)/ensemblediversity", ensemblediversity)
                h5write(h5file, "$(path)/probarr", probarr[:], "blosc", 3)
                h5write(h5file, "$(path)/I", I, "blosc", 3)
                h5write(h5file, "$(path)/J", J, "blosc", 3)
                h5write(h5file, "$(path)/V", V, "blosc", 3)
                h5write(h5file, "$(path)/version", rnafold_version)
                if circular
                    h5write(h5file, "$(path)/conformation", "circular")
                else
                    h5write(h5file, "$(path)/conformation", "linear")
                end
            end
        end

        ret = h5open(h5file, "r") do file
                   (convert(AbstractString, read(file["$(path)"]["sequence"])),
                    convert(Array{Int,1}, read(file["$(path)"]["mfepaired"])),
                    convert(Float64, read(file["$(path)"]["freeenergy"])),
                    convert(Float64, read(file["$(path)"]["ensembleprob"])),
                    convert(Float64, read(file["$(path)"]["ensemblediversity"])),
                    convert(Array{Float64,1}, read(file["$(path)"]["probarr"])[:]),
                    convert(Array{Int,1}, read(file["$(path)"]["I"])),
                    convert(Array{Int,1}, read(file["$(path)"]["J"])),
                    convert(Array{Float64,1}, read(file["$(path)"]["V"])),
                    convert(AbstractString, read(file["$(path)"]["version"])),
                    convert(AbstractString, read(file["$(path)"]["conformation"])))
        end
        return ret
    end

    function rnafold_probarr(inseq::AbstractString, temperature::Float64=37.0, circular::Bool=false, h5file::AbstractString="rnafoldcache.h5")
        path = rnafold(inseq, temperature, circular, h5file)
        probarr = h5open(h5file, "r") do file
            g = file["$(path)"]["probarr"]
            convert(Array{Float64,1}, read(g)[:])
        end
        return probarr
    end

    function rnafold_mfepaired(inseq::AbstractString, temperature::Float64=37.0, circular::Bool=false, h5file::AbstractString="rnafoldcache.h5")
        path = rnafold(inseq, temperature, circular, h5file)
        mfepaired = h5open(h5file, "r") do file
            g = file["$(path)"]["mfepaired"]
            convert(Array{Int,1}, read(g))
        end
        return mfepaired
    end

    function rnafold_basepairprobs(inseq::AbstractString, temperature::Float64=37.0, circular::Bool=false, h5file::AbstractString="rnafoldcache.h5")
        path = rnafold(inseq, temperature, circular, h5file)
        I,J,V = h5open(h5file, "r") do file
            len = length(convert(Array{Int,1}, read(file["$(path)"]["mfepaired"])))
            I = convert(Array{Int,1}, read(file["$(path)"]["I"]))
            J = convert(Array{Int,1}, read(file["$(path)"]["J"]))
            V = convert(Array{Float64,1}, read(file["$(path)"]["V"]))
            I,J,V
        end
        return I,J,V
    end

    function rnafold_posteriordecoding(inseq::AbstractString, temperature::Float64=37.0, circular::Bool=false, h5file::AbstractString="rnafoldcache.h5")
        seq = CommonUtils.canonicalise_sequence(inseq)
        len = length(seq)
        path = rnafold(inseq, temperature, circular, h5file)
        postdecoding = h5open(h5file, "cw") do file
            if exists(file, "$(path)/posteriordecodingpaired")
                return convert(Array{Int,1}, read(file["$(path)"]["posteriordecodingpaired"]))
            else
                I = convert(Array{Int,1}, read(file["$(path)"]["I"]))
                J = convert(Array{Int,1}, read(file["$(path)"]["J"]))
                V = convert(Array{Float64,1}, read(file["$(path)"]["V"]))
                basepairprobs = Array(sparse(I,J,V,len,len) + sparse(J,I,V,len,len))
                singleprobs = 1.0 - sum(basepairprobs,1)[:]
                postdecoding = RNATools.getPosteriorDecodingConsensusStructure(basepairprobs, singleprobs)
                h5write(h5file, "$(path)/posteriordecodingpaired", postdecoding, "blosc", 3)
                return postdecoding
            end
        end
        return postdecoding
    end
end

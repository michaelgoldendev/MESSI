push!(LOAD_PATH,joinpath(@__DIR__))
#using InsideCUDA
include("InsideCUDA.jl")


len = 7200
unpaired = zeros(Float32, len)
paired = zeros(Float32, len, len)
for z=1:1000
	println("$(z)")
	computeinsidecuda(unpaired, paired, KH99())
end
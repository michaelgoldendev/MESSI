using CUDAdrv, CUDAnative, CuArrays

function largeallocation()
	dev = CuDevice(0)
	ctx = CuContext(dev)
	largearray = CuArray(zeros(Float32, 56250000‬))
	finalize(largearray)
	GC.gc()
	CuArrays.pool_status()
	destroy!(dev)
end

for z=1:1000
	println("iter $z")
	largeallocation()
end
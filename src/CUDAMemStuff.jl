#using CUDAdrv
import CUDAdrv: @apicall, CuStream_t, CuDevice_t


function free(buf::Mem.Buffer)
    if buf.ptr != C_NULL
        @apicall(:cuMemFree, (Ptr{Cvoid},), buf.ptr)
    end
    return
end

"""
    upload!(dst::Mem.Buffer, src, nbytes::Integer, [stream=CuDefaultStream()]; async=false)

Upload `nbytes` memory from `src` at the host to `dst` on the device.
"""
function upload!(dst::Mem.Buffer, src::Ref, nbytes::Integer,
                 stream::CuStream=CuDefaultStream(); async::Bool=false)
    if async
        @apicall(:cuMemcpyHtoDAsync,
                 (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t, CuStream_t),
                 dst, src, nbytes, stream)
    else
        @assert stream==CuDefaultStream()
        @apicall(:cuMemcpyHtoD,
                 (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t),
                 dst, src, nbytes)
    end
end

"""
    download!(dst::Ref, src::Mem.Buffer, nbytes::Integer, [stream=CuDefaultStream()]; async=false)

Download `nbytes` memory from `src` on the device to `src` on the host.
"""
function download!(dst::Ref, src::Mem.Buffer, nbytes::Integer,
                   stream::CuStream=CuDefaultStream(); async::Bool=false)
    if async
        @apicall(:cuMemcpyDtoHAsync,
                 (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t, CuStream_t),
                 dst, src, nbytes, stream)
    else
        @assert stream==CuDefaultStream()
        @apicall(:cuMemcpyDtoH,
                 (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t),
                 dst, src, nbytes)
    end
end


## array based

"""
    alloc(src::AbstractArray)

Allocate space to store the contents of `src`.
"""
#=
function alloc(src::AbstractArray)
    return alloc(sizeof(src))
end=#

"""
    upload!(dst::Mem.Buffer, src::AbstractArray, [stream=CuDefaultStream()]; async=false)

Upload the contents of an array `src` to `dst`.
"""
function upload!(dst::Mem.Buffer, src::AbstractArray,
                 stream=CuDefaultStream(); async::Bool=false)
    upload!(dst, Ref(src, 1), sizeof(src), stream; async=async)
end

"""
    upload(src::AbstractArray)::Mem.Buffer

Allocates space for and uploads the contents of an array `src`, returning a Mem.Buffer.
Cannot be executed asynchronously due to the synchronous allocation.
"""
function upload(src::AbstractArray)
    dst = alloc(src)
    upload!(dst, src)
    return dst
end

"""
    download!(dst::AbstractArray, src::Mem.Buffer, [stream=CuDefaultStream()]; async=false)

Downloads memory from `src` to the array at `dst`. The amount of memory downloaded is
determined by calling `sizeof` on the array, so it needs to be properly preallocated.
"""
function download!(dst::AbstractArray, src::Mem.Buffer,
                   stream::CuStream=CuDefaultStream(); async::Bool=false)
    ref = Ref(dst, 1)
    download!(ref, src, sizeof(dst), stream; async=async)
    return
end

"""
    alloc(T::Type, [count::Integer=1])

Allocate space for `count` objects of type `T`.
"""
function alloc(::Type{T}, count::Integer=1) where {T}
    check_type(Mem.Buffer, T)

    return alloc(sizeof(T)*count)
end

"""
    download(::Type{T}, src::Mem.Buffer, [count::Integer=1], [stream=CuDefaultStream()]; async=false)::Vector{T}

Download `count` objects of type `T` from the device at `src`, returning a vector.
"""
function download(::Type{T}, src::Mem.Buffer, count::Integer=1,
                  stream::CuStream=CuDefaultStream(); async::Bool=false) where {T}
    dst = Vector{T}(undef, count)
    download!(dst, src, stream; async=async)
    return dst
end
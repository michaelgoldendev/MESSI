logout = nothing
logiter = 1

function newlog(newlogfilename)
  global logout
  closelogger()
  logfilename = newlogfilename
  logout = open(logfilename, "w")
end

function logtime(s::AbstractString, t::Float64, numcols::Int)
  global logout
  if logout != nothing
    write(logout, string(s, "\t", t, "\t", nprocs(), "\t", numcols,"\n"))
    global logiter
    if logiter % 25 == 0
      flush(logout)
    end
    logiter += 1
  end
end

function closelogger()
  global logout
  if logout != nothing
    close(logout)
    logout = nothing
  end
end

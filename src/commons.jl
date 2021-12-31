# commons.jl -- utility functions common in both implementations, i think

@inline function getbin(pos::StaticVector{N,T},
                base::BoundBox{N,T},
                cutoff::T;
                skipcheck=false) where {N,T}
    offset = pos .- base.start
    sz = size(base)
    skipcheck || @assert all(offset[i] <= sz[i], 1:N)
    return ceil.(Int, offset ./ cutoff) |> CartesianIndex
end

@inline function locflag(o::NTuple{N}, ref::NTuple{N}) where N
    ntuple(N) do i
        o[i] == 1 ? -1 : o[i] == ref[i] ? 1 : 0
    end
end

@inline function probeset(o::NTuple{N}) where N
    ntuple(N) do i
        o[i] == -1 ? (0:1) : o[i] == 1 ? (-1:0) : (-1:1)
    end |> CartesianIndices
end

@inline function wrapprobe(o::NTuple{N}, ref::NTuple{N}) where N
    ntuple(N) do i
        o[i] == 0 ? 1 : o[i] == ref[i]+1 ? 1 : o[i]
    end |> CartesianIndex
end

### THREADED MAP

"""
    tmap(f, src...; basesize, nt)

Threaded version of `map`. For use with relatively expensive `f`.
The function is allocating tasks for each collection element.

Keyword Arguments
===
+ `basesize`: length of array for when to fall back to serial. Choice
can depend on the load of `f`. Defaults to ``10000``.
+ `nt`: number of threads to be used. Defaults to `Threads.nthreads()`.
"""
function tmap(f, src...; basesize=10_000, nt=Threads.nthreads())
    ## if source is not large enough compared to basesize
    ## fall back to serial map
    length(src) < basesize && return map(f, src...)

    sem = Base.Semaphore(nt)
    return map(src...) do (x...)
        # set up tasks for each computation
        # basesize should depend on how heavy f is
        @async begin
            Base.acquire(sem)
            __t = Threads.@spawn f(x...)
            res = fetch(__t)
            Base.release(sem)
            return res
        end
    end .|> fetch
end

"""
    tmap!(f, dst, src...; basesize, nt)

Threaded version of `map!`. For use with relatively expensive `f` and
when `dst` can be provided (type is known upon compilation).

Keyword Arguments
===
+ `basesize`: length of array for when to fall back to serial. Choice
can depend on the load of `f`. Defaults to ``10000``.
+ `nt`: number of threads to be used. Defaults to `Threads.nthreads()`.

See also: [tmap][@ref]
"""
function tmap!(f, dst, src...; basesize=10_000, nt=Threads.nthreads())
    ## if source is not large enough compared to basesize
    ## fall back to serial map!
    length(dst) < basesize && return map!(f, dst, src...)

    sem = Base.Semaphore(nt)

    map!(dst, src...) do (x...)
        @async begin
            Base.acquire(sem)
            __t = Threads.@spawn f(x...)
            res = fetch(__t)
            Base.release(sem)
            return res
        end |> fetch
    end
end ## mostly untested so far; not safe

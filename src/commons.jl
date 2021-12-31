# commons.jl -- utility functions common in both implementations, i think

@inline function getbin(pos::StaticVector{N,T},
                base::BoundBox{N,T},
                cutoff::T;
                skipcheck=false) where {N,T}
    offset = pos .- base.start
    sz = size(base)
    skipcheck || @assert all(offset[i] <= sz[i], 1:N)
    return ceil.(Int, offset ./ cutoff).data |> CartesianIndex
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
    tmap(f, src, shaped=false; basesize, nt)

Threaded version of `map`. For use with relatively expensive `f`.
The function is allocating tasks per available thread.

# Keyword Arguments
+ `basesize`: length of array for when to fall back to serial. Choice
can depend on the load of `f`. Defaults to ``10000``.
+ `nt`: number of threads to be used. Defaults to `Threads.nthreads()`.
"""
function tmap(f, src, shaped=false; basesize=10_000, nt=Threads.nthreads())
    len = length(src)
    len < basesize && return map(f, src)
    sz = size(src)
    idcs = eachindex(src)

    tasks = map(1:nt) do tid
        ulen, rem = divrem(len, nt)
        # divide evenly and then dump into one by one
        fst = ulen * (tid-1) + first(idcs)
        lst = fst + ulen - 1
        if rem > 0
            if tid < rem
                fst += tid-1    # from dumping
                lst += tid      # from dumping
            else
                fst += rem      # from dumping
                lst += rem      # from dumping
            end
        end
        Threads.@spawn map(f, view(src, fst:lst))
    end

    unshaped = fetch.(tasks) |> q->vcat(q...)
    shaped && return reshape(unshaped, sz)
    return unshaped
end

"""
    tmap!(f, dst, src...; basesize, nt)

Threaded version of `map!`. For use with relatively expensive `f` and
when `dst` can be provided (type is known upon compilation).

# Keyword Arguments
+ `basesize`: length of array for when to fall back to serial. Choice
can depend on the load of `f`. Defaults to ``10000``.
+ `nt`: number of threads to be used. Defaults to `Threads.nthreads()`.

See also: [tmap][@ref]
"""
function tmap!(f, dst, src; basesize=10_000)
    ## if source is not large enough compared to basesize
    ## fall back to serial map!
    @assert length(dst) == length(src)
    len < basesize && return map!(f, dst, src)

    Threads.@threads for idx in eachindex(src)
        dst[idx] = f(src[idx])
    end
end ## mostly untested so far; not safe

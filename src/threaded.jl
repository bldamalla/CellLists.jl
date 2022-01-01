# threaded.jl -- threaded implementation of those in serial.jl

## the only things to be threaded are cacheing the 
## neighbor lists and probably CAL construction

function tCompressedAdjacencyList(pos::Vector{<:StaticVector{N}},
                                  bdict,
                                  base=getboundbox(pos);
                                  skipcheck=false,
                                  basesize=2000)
    ## get cache
    ## map and then reduce with vcat?
    ncache = tneighborcache(bdict)
    guarded = bdict isa GuardedBinDict
    peuc(a, b) = guarded ? peuclidean(a, b, bdict.guard) : nothing
    distf = guarded ? peuc : euclidean
    cutoff = bdict.cutoff

    tp = tmap(eachindex(pos); basesize=basesize) do idx
        cell = getbin(pos[idx], base, cutoff; skipcheck=skipcheck)
        filter(ncache[cell]) do q
            distf(pos[q], pos[idx]) < cutoff
        end
    end

    CompressedAdjacencyList(length.(tp), reduce(vcat, tp))
end

function tneighborcache(bdict::UnguardedBinDict, sort=false;
                        basesize=2000)
    sz = size(bdict); cidcs = CartesianIndices(sz)

    cache = tmap(cidcs; basesize=basesize) do idx
        prbset = locflag(idx.I, sz) |> probeset
        map(prbset) do probe
            bdict[CartesianIndex(probe.I .+ idx.I)]
        end |> q->vcat(q...)
    end

    if sort ## sort in place -- is there a better way to do this?
        sort!.(cache)
    end

    return cache
end

function tneighborcache(bdict::GuardedBinDict, sort=false;
                        basesize=2000)
    sz = size(bdict); cidcs = CartesianIndices(sz)

    around = let nones = ntuple(_->1, ndims(bdict))
        CartesianIndices(nones .* 3) .- CartesianIndex(nones .* 2)
    end

    cache = tmap(cidcs; basesize=basesize) do idx
        map(around) do probe
            bdict[wrapprobe(probe.I .+ idx.I, sz)]
        end |> q->vcat(q...)
    end

    if sort # sort inplace --- is there a better way to do this?
        sort!.(cache)
    end
    
    return cache
end

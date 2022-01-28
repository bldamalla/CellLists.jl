# serial.jl -- for serial implementation of binning

## TODO: given the positions, bounding box, and the bin dict
## update the bin dictionary

function fillbins!(bdict::BinDict{N},
                   pos::AbstractVector{StaticVector{N}},
                   base=getboundbox(pos);
                   skipcheck=false) where N
    # first is to check if bounding box of current positions
    # is applicable with current bdict
    cutoff = bdict.cutoff
    @assert size(bdict) == ceil.(Int, size(base) ./ cutoff)

    # next is iterate over the whole array and put in each bin
    @inbounds for idx in eachindex(pos)
        cell = getbin(pos[idx], base, cutoff; skipcheck=skipcheck)
        push!(bdict[cell], idx)
    end

    return bdict
end

function emptybins!(bdict::BinDict)
    for index in eachindex(bdict.bins)
       empty!(bdict.bins[index])
    end
    return bdict
end

function CompressedAdjacencyList(pos::Vector{<:StaticVector},
                                 bdict,
                                 base=getboundbox(pos);
                                 skipcheck=false)
    ## get cache
    ## map and then reduce with vcat
    ncache = neighborcache(bdict)
    guard = bdict isa GuardedBinDict ? bdict.guard : nothing
    peuc(a, b) = peuclidean(a, b, guard)    # never used when guard is nothing
    distf = bdict isa GuardedBinDict ? peuc : euclidean
    cutoff = bdict.cutoff

    tp = @inbounds map(eachindex(pos)) do idx
        cell = getbin(pos[idx], base, cutoff; skipcheck=skipcheck)
        filter(ncache[cell]) do q
            distf(pos[idx], pos[q]) < cutoff && idx != q
        end
    end

    return CompressedAdjacencyList(length.(tp), reduce(vcat, tp))
end

function neighborcache(bdict::UnguardedBinDict, sort=false)
    sz = size(bdict); cidcs = CartesianIndices(sz)

    cache = map(cidcs) do idx
        prbset = locflag(idx.I, sz) |> probeset
        map(prbset) do probe
            bdict[CartesianIndex(probe.I .+ idx.I)]
        end |> q->vcat(q...)
    end

    if sort # sort inplace --- is there a better way to do this?
        sort!.(cache)
    end

    return cache
end

function neighborcache(bdict::GuardedBinDict, sort=false)
    sz = size(bdict); cidcs = CartesianIndices(sz)

    around = let nones = ntuple(_->1, ndims(bdict))
        CartesianIndices(nones .* 3) .- CartesianIndex(nones .* 2)
    end

    cache = map(cidcs) do idx
        map(around) do probe
            bdict[wrapprobe(probe.I .+ idx.I, sz)]
        end |> q->vcat(q...)
    end

    if sort # sort inplace --- is there a better way to do this?
        sort!.(cache)
    end
    
    return cache
end

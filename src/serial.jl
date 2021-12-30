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

function neighborcache(bdict::UnguardedBinDict, sort=false)
    sz = size(bdict); cidcs = CartesianIndices(sz)

    cache = map(cidcs) do idx
        prbset = locflag(idx.I, sz) |> probeset
        vcat([bdict[CartesianIndex(prb.I .+ idx.I)] for prb in prbset]...)
    end

    if sort # sort inplace --- is there a better way to do this?
        map!(sort!, cache, cache)
    end

    return cache
end

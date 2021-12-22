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
    @assert bdict.dims == ceil.(Int, size(base) ./ cutoff)

    # next is iterate over the whole array and put in each bin
    @inbounds for idx in eachindex(pos)
        cell = getbin(pos[idx], base, cutoff; skipcheck=skipcheck)
        push!(bdict.bins[cell], idx)
    end

    return bdict
end

function emptybins!(bdict::BinDict)
    for key in keys(bdict.bins)
       empty!(bdict.bins[key])
    end
    return bdict
end

function getbin(pos::StaticVector{N,T},
                base::BoundBox{N,T},
                cutoff::T;
                skipcheck=false) where {N,T}
    offset = pos .- base.start
    sz = size(base)
    skipcheck || @assert all(offset[i] < sz[i], 1:N)
    return ceil.(Int, offset ./ cutoff)
end

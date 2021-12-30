# commons.jl -- utility functions common in both implementations, i think

@inline function getbin(pos::StaticVector{N,T},
                base::BoundBox{N,T},
                cutoff::T;
                skipcheck=false) where {N,T}
    offset = pos .- base.start
    sz = size(base)
    skipcheck || @assert all(offset[i] <= sz[i], 1:N)
    return ceil.(Int, offset ./ cutoff)
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

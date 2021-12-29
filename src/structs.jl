# structs.jl -- structures for handling the cell lists

## the primary structure I'll need is a set of bins arranged in a rectangular
## array.

abstract type BinDict{N} end

struct BoundBox{N,T}
    start::SVector{N,T}
    stop::SVector{N,T}
end
Base.size(bbox::BoundBox) = bbox.stop .- bbox.start

"""
    UnguardedBinDict{N,T}

An ``N``-dimensional binning dictionary for making cell-based neighbor lists.
Use this when there are no periodic boundary conditions to be followed.
"""
struct UnguardedBinDict{N,T} <: BinDict{N}
    bins::Array{Vector{Int},N}
    cutoff::T
end

"""
    GuardedBinDict{N}

An ``N``-dimensional binning dictionary for making cell-based neighbor lists.
This follows periodic boundary conditions as prescribed by a `guard` box.
"""
struct GuardedBinDict{N,T} <: BinDict{N}
    bins::Array{Vector{Int},N}
    cutoff::T
    guard::SVector{N,T}
end

Base.size(bdict::BinDict) = size(bdict.bins)

# TODO: something about the constructors of each
## you can have a base box if it is provided, but if it is not you can
## likely guess from the input information.

## allow construction from arrays of SVectors or MVectors
## use the abstract type StaticVector{N} where N

function createbindict(pos::AbstractVector{StaticVector{N}}, cutoff;
                       base=getboundbox(pos), guard=nothing) where N
    ## TODO: get the maximum and minimum position
    ## enclose this in a structure?
    ## get a temp grid size and create bins that way
    @assert base isa BoundBox "base should be a BoundBox"

    # from now get the grid
    sz = size(base)
    dims = ceil.(Int, sz ./ cutoff)

    # TODO: based on the grid dimensions create binning dictionary
    bdict = [Int[] for _ in CartesianIndices(dims)]

    if guard isa Nothing
        # create an unguarded dict
        return UnguardedBinDict{N}(bdict, cutoff)
    else
        # create a guarded dict
        return GuardedBinDict{N}(bdict, cutoff, guard)
    end
end

function getboundbox(pos::AbstractVector{StaticVector{N}}) where N
    T = eltype(eltype(pos))
    _mins = SVector{N,T}(minimum(p[n] for p in pos) for n in 1:N)
    _maxs = SVector{N,T}(maximum(p[n] for p in pos) for n in 1:N)
    return BoundBox(_mins, _maxs)
end

## define a getindex for the bounding boxes using Cartesian indices
getindex(bdict::BinDict{N}, idx) where N = bdict.bins[idx]

#### ADJACENCY LISTS ###########

"""
    CompressedAdjacencyList

Represents an iterable adjacency list. Contains information on the number of
neighbors of each element and their respective neighbors.
"""
struct CompressedAdjacencyList
    valences::Vector{Int}
    neighbors::Vector{Int}

    function CompressedAdjacencyList(valences, neighbors)
        @assert sum(valences) == length(neighbors)
        return new(valences, neighbors)
    end
end

# state can act like the indices of a matrix, the first is a placeholder
# for a total number of elements that have been passed
# the second element is the index of
# the current indexed element and the second is the nth neighbor
function Base.iterate(cal::CompressedAdjacencyList, state=(0,1,1))
    s, a, b = state
    s == length(cal.neighbors) && return nothing

    a_ = b == cal.valences[a] ? a+1 : a
    b_ = b == cal.valences[a] ? 1 : b+1

    return (a, cal.neighbors[s+1]), (s+1, a_, b_)
end

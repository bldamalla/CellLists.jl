# structs.jl -- structures for handling the cell lists

## the primary structure I'll need is a set of bins arranged in a rectangular
## array.

abstract type BinDict end

struct BoundBox{N,T}
    start::SVector{N,T}
    stop::SVector{N,T}
end
Base.size(bbox::BoundBox) = bbox.stop .- bbox.start

"""
    UnguardedBinDict{N}

An ``N``-dimensional binning dictionary for making cell-based neighbor lists.
Use this when there are no periodic boundary conditions to be followed.
"""
struct UnguardedBinDict{N,T} <: BinDict
    bins::Dict{NTuple{N,Int},Vector{Int}}
    cutoff::T
    dims::NTuple{N,Int}
end

"""
    GuardedBinDict{N}

An ``N``-dimensional binning dictionary for making cell-based neighbor lists.
This follows periodic boundary conditions as prescribed by a `guard` box.
"""
struct GuardedBinDict{N,T} <: BinDict
    bins::Dict{NTuple{N,Int},Vector{Int}}
    cutoff::T
    dims::NTuple{N,Int}
    guard::SVector{N,T}
end

# TODO: something about the constructors of each
## you can have a base box if it is provided, but if it is not you can
## likely guess from the input information.

## allow construction from arrays of SVectors or MVectors
## use the abstract type StaticVector{N} where N

function createbindict(pos::AbstractVector{StaticVector{N}}, cutoff;
                       base=nothing, guard=nothing) where N
    ## TODO: get the maximum and minimum position
    ## enclose this in a structure?
    ## get a temp grid size and create bins that way
    if base isa Nothing
        # get bounding box from the given input
        base = getboundbox(pos)
    end
    @assert base isa BoundBox "base should be a BoundBox"

    # from now get the grid
    sz = size(base)
    dims = sz ./ cutoff

    # TODO: based on the grid dimensions create binning dictionary
    bdict = Dict{NTuple{N,Int},Vector{Int}}()
    for index in CartesianIndices(dims)
        bdict[index.I] = Vector{Int}()
    end

    if guard isa Nothing
        # create an unguarded dict
        return UnguardedBinDict{N}(bdict, cutoff, dims)
    else
        # create a guarded dict
        return GuardedBinDict{N}(bdict, cutoff, dims, guard)
    end
end

function getboundbox(pos::AbstractVector{StaticVector{N}}) where N
    T = eltype(eltype(pos))
    _mins = SVector{N,T}(minimum(p[n] for p in pos) for n in 1:m)
    _maxs = SVector{N,T}(maximum(p[n] for p in pos) for n in 1:m)
    return BoundBox(_mins, _maxs)
end
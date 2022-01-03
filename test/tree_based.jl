## tree_based.jl -- for tree based distance finder

## can be used to check if the created cell based method is accurate
## i guess we can use a similar implementation in Molly.jl
## use the one without threading

using NearestNeighbors, StaticArrays, Distances

function treeadjlist(A::Vector{<:StaticVector{N}},
                     guard::SVector{N}=ntuple(_->Inf, N),
                     cutoff) where N
    btree = BallTree(A, PeriodicEuclidean(guard))
    nlist = Tuple{Int,Int}[]

    for i in eachindex(A)
        # get those within cutoff using the ball tree
        for j in inrange(btree, A[i], cutoff, true)
            if i > j
                push!(nlist, (i, j))
            end
        end
    end

    return nlist
end

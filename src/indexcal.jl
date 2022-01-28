# indexcal.jl -- tools for CAL indexing and threading?

@inline function viewsearch(a, src, i=firstindex(src), j=lastindex(src))
    # search for the largest index i of sorted vec src
    # where src[i] <= a --- use recursive binary search?
    i == j && return i
    1 <= a <= first(src) && return firstindex(src)
    mid = div(i+j, 2)

    @inbounds if a <= src[mid]
        return viewsearch(a, src, i, mid)
    else
        return viewsearch(a, src, mid+1, j)
    end
end

function Base.getindex(cal::CompressedAdjacencyList, i::Int)
    a_ix = viewsearch(i, cal.csvalences)
    return (a_ix, cal.neighbors[i])
end

struct CALview
    cal::CompressedAdjacencyList
    start::Int
    stop::Int
end
Base.firstindex(calv::CALview) = calv.start
Base.lastindex(calv::CALview) = calv.stop
Base.length(calv) = calv.stop - calv.start + 1

function Base.iterate(calv::CALview,
                      state=(calv.start,viewsearch(calv.start,
                                        calv.cal.csvalences,
                                        calv.start, calv.stop)))
    st, a = state
    st > calv.stop && return nothing
    vals = calv.cal.csvalences

    a_ = st == vals[a] ? a+1 : a

    return (a_, calv.cal.neighbors[st]), (st+1, a_)
end

function Base.getindex(cal::CompressedAdjacencyList, rng::UnitRange)
    return CALview(cal, rng.start, rng.stop)
end

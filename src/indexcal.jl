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

function Base.iterate(calv::CALview, state=(calv.start,0))
    st, a = state
    st > calv.stop && return nothing
    vals = calv.cal.csvalences
    _sta, _sto = calv.start, calv.stop

    a_ = a == 0 ? viewsearch(st, vals, _sta, _sto) : a
    a__ = st == vals[a_] ? a_+1 : a_

    return (a__, calv.cal.neighbors[st]), (st+1, a__)
end

function Base.getindex(cal::CompressedAdjacencyList, rng::UnitRange)
    return CALview(cal, rng.start, rng.stop)
end

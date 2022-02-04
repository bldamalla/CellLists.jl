## calinterface.jl --- interfacing with CAL structs

## something about further filtering --- i guess this is given when using `Filter`s
## build neighbors per index
## interface for collect --- i guess this is given

"""
    calcollect(cal[, sieve])

Return a `Vector{Tuple{Int,Int}}` of neighbor pairs in the CAL `cal`. Optional
argument `sieve` is for further filtering of the indices. Defaults to a "no-filter"
returning _all_ pairs in the CAL.

Possible filters:
+ `:full` - Default. Similar to `collect(cal)`
+ `:<` - Returns pairs `(i, j)` where `i < j`
+ `:>` - Returns pairs `(i, j)` where `i > j`

See also: [calfilter][@ref]
"""
calcollect(cal, sieve=:full) = collect(calfilter(cal, sieve))

"""
    calfilter(cal[, sieve])

Return a `Generator` of tuples of indices in the CAL `cal`. Optional argument
`sieve` is for further filtering of the indices. Defaults to a "no-filter" returning
_all_ pairs in the CAL. This can be chosen over `calcollect` if only iteration is
needed. This includes reduction and maps and excludes indexing operations.

Possible filters:
+ `:full` - Default. Similar to iterating over CAL without filter.
+ `:<` - Returns pairs `(i, j)` where `i < j`
+ `:>` - Returns pairs `(i, j)` where `i > j`

*Note*: If the last two filters are not chosen, the default behavior is done.
However, note that iteration over the result is the same as iteration over the CAL
itself.

See also: [calcollect][@ref]
"""
function calfilter(cal::CompressedAdjacencyList, sieve=:full)
    _filter = _getfilter(sieve)
    isnothing(_filter) && return ((i, j) for (i, j) in cal)
    return ((i, j) for (i, j) in cal if _filter(i, j))
end

"""
    neighborlist(cal, i::Int[, sieve])

Return a `Vector{Int}` containing indices in the CAL `cal` that are considered
neighbors of the point with index `i`. An optional filter (`sieve`) can be added.
Equivalent to `collect(neighboritr(cal, i[, sieve]))`.

Possible filters:
+ `:full` - Default. Includes indices greater and lower than `i`.
+ `:<` - Includes indices lower than `i`.
+ `:>` - Includes indices above `i`.

*Note*: If the last two filters are not chosen, the default behavior is done.

See also: [neighboritr][@ref]
"""
neighborlist(cal, i, sieve=:full) = map(identity, neighboritr(cal, i, sieve))

"""
    neighboritr(cal, i::Int[, sieve])

Return a `Generator` of indices in the CAL `cal` that are considered
neighbors of the point with index `i`. An optional filter (`sieve`) can be added.
Use this instead of `neighborlist` if only iteration over the list is needed.
This includes reduction and maps and excludes indexing operations.

Possible filters:
+ `:full` - Default. Includes indices greater and lower than `i`.
+ `:<` - Includes indices lower than `i`.
+ `:>` - Includes indices above `i`.

*Note*: If the last two filters are not chosen, the default behavior is done.

See also: [neighborlist][@ref]
"""
function neighboritr(cal, i::Int, sieve=:full)
    start = cal.csvalences[i-1]+1
    stop = cal.csvalences[i]
    _filter = _getfilter(sieve)
    isnothing(_filter) && return (cal.neighbors[q] for q in start:stop)
    return (cal.neighbors[q] for q in start:stop if _filter(q, i))
end

@inline _getfilter(sieve) = ifelse(sieve == :<, <, ifelse(sieve == :>, >, nothing))



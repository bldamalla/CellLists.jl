module CellLists

using Distances
using StaticArrays

include("structs.jl")
include("commons.jl")
include("serial.jl")

include("threaded.jl")

end # module

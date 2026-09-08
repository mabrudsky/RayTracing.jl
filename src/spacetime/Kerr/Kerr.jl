module kerr
using StaticArrays
using Parameters
using ForwardDiff: Dual, Partials, partials

include(joinpath(@__DIR__, "..", "spacetime_types.jl"))


#------------------------------- spacetime --------------------------------- 
include("KerrSchildCoordinates.jl")

export metric, metric_inverse, christoffel

end # module kerr   
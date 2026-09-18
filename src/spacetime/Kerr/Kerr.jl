module Kerr
using StaticArrays
using Parameters
using ..RayTracing: AbstractSpacetime, AbstractKerrSpacetime, AbstractChristoffel, ChristoffelTag, _Γcomp, _christoffel

#------------------------------- spacetime --------------------------------- 
include("KerrSchildCoordinates.jl")

export metric, metric_inverse, christoffel

end # module kerr   
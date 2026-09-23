module kerr
using StaticArrays
using Parameters
using ..RayTracing: AbstractSpacetime, AbstractKerrSpacetime, AbstractCoordinate, CartesianCoordinate, SphericalCoordinate
using ..RayTracing: AbstractChristoffel, ChristoffelTag, _Γcomp, _christoffel

import ..RayTracing: metric, metric_inverse, christoffel # Bring in "metric", "metric_inverse" and "christoffel" so this module adds a new case to each instead of accidentally creating a separate, disconnected function.
#------------------------------- spacetime --------------------------------- 
include("KerrSchildCoordinates.jl")

export metric, metric_inverse, christoffel

end # module kerr   
module RayTracing
using StaticArrays
using Parameters
import ForwardDiff: Dual, Partials, partials

include("spacetime/spacetime.jl")

# Export spacetime module
export Kerr

# Export spacetime functions
export metric, metric_inverse, christoffel

end
module RayTracing
import StaticArrays: @SArray, SVector
import Parameters: @with_kw
import ForwardDiff: Dual, Partials, partials

include("spacetime/spacetime.jl")

# Export spacetime module
export Kerr

# Export spacetime functions
export metric, metric_inverse, christoffel

end
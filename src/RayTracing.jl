module RayTracing
import StaticArrays: @SArray, SVector
import Parameters: @with_kw
import ForwardDiff: Dual, Partials, partials
using KernelAbstractions

include("RayTracing_types.jl")
include("configuration/configuration.jl")
include("spacetime/spacetime.jl")
include("propagate/propagate.jl")   

# Export spacetime module
export Kerr

# Export spacetime functions
export metric, metric_inverse, christoffel

end
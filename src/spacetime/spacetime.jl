include("spacetime_types.jl")
include("christoffel.jl")

#= Declare "metric" and "metric_inverse" here first, with no
   implementation yet, so that every spacetime module (Kerr,
   Schwarzschild, etc.) adds to these SAME functions instead of
   each one creating its own.
=#
function metric end
function metric_inverse end

include("Kerr/Kerr.jl")

using .kerr
abstract type AbstractRayTracing end

abstract type AbstractPropagate <: AbstractRayTracing end
abstract type AbstractEnsembleProblem <: AbstractPropagate end
abstract type AbstractSpacetime <: AbstractRayTracing end
abstract type AbstractConfiguration <: AbstractRayTracing end
abstract type AbstractPhysicalModel <: AbstractRayTracing end
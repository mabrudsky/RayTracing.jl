# spacetime types
abstract type AbstractKerrSpacetime <: AbstractSpacetime end


# topology types
abstract type AbstractCoordinate <: AbstractRayTracing end
struct CartesianCoordinate <: AbstractCoordinate end
struct SphericalCoordinate <: AbstractCoordinate end
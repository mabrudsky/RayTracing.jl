abstract type AbstractConfiguration end
abstract type AbstractParameter <: AbstractConfiguration end

#-------------------------------------------------------------------------
@kwdef struct Configuration{N,T<:Union{Float32,Float64}} <: AbstractConfiguration 
    spacetime :: AbstractSpacetime
    model     :: AbstractModel
    U0        :: SVector{N,T} 
    p         :: SVector{N,T}
    tspan     :: SVector{2,T}
end 
@kwdef struct Configuration{N,M,T<:Union{Float32,Float64}} <: AbstractConfiguration 
    spacetime :: AbstractSpacetime
    model     :: AbstractEnsembleProblem
    backend   :: KernelAbstractions.Backend
    u0        :: SVector{N,T} 
    p         :: SVector{M,T}
    tspan     :: Tuple{T,T}
end 
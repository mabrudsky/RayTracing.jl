@kwdef struct IntegrationParameters{S <: AbstractSpacetime, M <: AbstractPhysicalModel} <: AbstractConfiguration
    spacetime::S
    model::M
end

@kwdef struct Configuration{N, T <: Union{Float32,Float64}, S <: AbstractSpacetime, E <: AbstractEnsembleProblem, B <: KernelAbstractions.Backend, M <: AbstractPhysicalModel} <: AbstractConfiguration
    spacetime::S
    equations::E
    backend::B
    u0::SVector{N,T}
    model::M
    tspan::Tuple{T,T}
end

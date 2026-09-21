struct GeodesicEquationsVacuum <: AbstractEnsembleProblem end
struct ParallelTransportEquationsVacuum <: AbstractEnsembleProblem end
#-------------------------------------------------------------------------
function Ensemble_Problem(configuration::AbstractConfiguration{N,M,T}, initial_data::AbstractMatrix)  where {N,M,T<:Union{Float32,Float64}} 
    size(initial_data, 1) == N || throw(DimensionMismatch("Each column of initial_data must be an initial condition. Expected size $N, got $(size(initial_data, 1))"))    
    return _Ensemble_Problem(configuration, configuration.model, initial_data)
end

function _Ensemble_Problem(configuration::AbstractConfiguration{N,M,T}, ::GeodesicEquationsVacuum, initial_data::AbstractMatrix) where {N,M,T<:Union{Float32,Float64}}
    (; u0, p, tspan) = configuration
    prob         = ODEProblem{false}(geodesic_equations_Vacuum, u0, tspan, p)
    prob_func    = (prob, i, repeat) -> remake(prob, u0 = SVector{8, T}(initial_data[:,i]))
    output_func  = (sol, i) -> (sol[end], false)
    return EnsembleProblem(prob, prob_func=prob_func, output_func=output_func, safetycopy=false) 
end

function _Ensemble_Problem(configuration::AbstractConfiguration{N,M,T}, ::ParallelTransportEquationsVacuum, initial_data::AbstractMatrix) where {N,M,T<:Union{Float32,Float64}}
    (; u0, p, tspan) = configuration
    prob         = ODEProblem{false}(parallel_transport_equations_Vacuum, u0, tspan, p)
    prob_func    = (prob, i, repeat) -> remake(prob, u0 = SVector{16, T}(initial_data[:,i]))
    output_func  = (sol, i) -> (sol[end], false)
    return EnsembleProblem(prob, prob_func=prob_func, output_func=output_func, safetycopy=false) 
end

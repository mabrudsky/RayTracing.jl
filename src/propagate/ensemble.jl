struct GeodesicEquationsVacuum <: AbstractEquationSet end
struct ParallelTransportEquationsVacuum <: AbstractEquationSet end
#-------------------------------------------------------------------------
function Ensemble_Problem(configuration::Configuration{N,T,S,E,B,M}, initial_data::AbstractMatrix) where {N, T <: Union{Float32,Float64}, S <: AbstractSpacetime, E <: AbstractEquationSet, B <: KernelAbstractions.Backend, M <: AbstractPhysicalModel}
    size(initial_data, 1) == N || throw(DimensionMismatch("Each column of initial_data must be an initial condition. Expected size $N, got $(size(initial_data, 1))"))    
    eltype(initial_data) === T || throw(ArgumentError("initial_data must have element type $T"))
    return _Ensemble_Problem(configuration, configuration.equations, initial_data)
end

function _Ensemble_Problem(configuration::Configuration{N,T,S,E,B,M}, ::GeodesicEquationsVacuum, initial_data::AbstractMatrix) where {N, T <: Union{Float32,Float64}, S <: AbstractSpacetime, E <: AbstractEquationSet, B <: KernelAbstractions.Backend, M <: AbstractPhysicalModel}
    (; u0, spacetime, model, tspan) = configuration
    p = IntegrationParameters(spacetime, model)
    prob        = ODEProblem{false}(geodesic_equations_Vacuum, u0, tspan, p)
    prob_func   = (prob, i, repeat) -> remake(prob, u0 = SVector{N, T}(@view initial_data[:,i]))
    output_func = (sol, i) -> (sol[end], false)
    return EnsembleProblem(prob, prob_func=prob_func, output_func=output_func, safetycopy=false) 
end

function _Ensemble_Problem(configuration::Configuration{N,T,S,E,B,M}, ::ParallelTransportEquationsVacuum, initial_data::AbstractMatrix) where {N, T <: Union{Float32,Float64}, S <: AbstractSpacetime, E <: AbstractEquationSet, B <: KernelAbstractions.Backend, M <: AbstractPhysicalModel}
    (; u0, spacetime, model, tspan) = configuration
    p = IntegrationParameters(spacetime, model)
    prob        = ODEProblem{false}(parallel_transport_equations_Vacuum, u0, tspan, p)
    prob_func   = (prob, i, repeat) -> remake(prob, u0 = SVector{N, T}(@view initial_data[:,i]))
    output_func = (sol, i) -> (sol[end], false)
    return EnsembleProblem(prob, prob_func=prob_func, output_func=output_func, safetycopy=false) 
end

abstract type AbstractEnsembleProblem end
#abstract type AbstractEnsembleCPUBackend end
#abstract type AbstractEnsembleGPUBackend end

function ensemble_problem(initial_data::AbstractInitialData, configuration::AbstractConfiguration) #; kwargs...)
    u0 = initial_data.u0
    p = configuration.parameters
end
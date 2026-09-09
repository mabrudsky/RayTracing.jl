module RayTracing
include("spacetime/spacetime.jl")

# Export spacetime module
export Kerr

# Export spacetime functions
export metric, metric_inverse, christoffel

end


#= Esto 

a = @SVector [- sum(Γ[μ, ν, λ] * v[μ] * v[ν] for ν in 1:4, μ in 1:4) for λ in 1:4]

es equivalente a 

# El bucle interno sobre μ lee la memoria de forma estrictamente lineal
for λ in 1:4
    s = zero(T)
    for ν in 1:4
        @simd for μ in 1:4  # Dimensión 1 -> Lectura lineal contigua (Stride 1)
            s += Γ[μ, ν, λ] * v[μ] * v[ν]
        end
    end
    a[λ] = -s
end

=#
using RayTracing.Kerr: KerrSchildCoordinates, metric, christoffel
using ForwardDiff: Dual, Partials, partials
using StaticArrays

function covariant_derivative_metric(k, u::SVector{4,T}) where T
    Γ = christoffel(k, u)
    g = metric(k, u)

    xd = SVector{4, Dual{Nothing,T,4}}(ntuple(i -> Dual{Nothing,T,4}(u[i], Partials(ntuple(j -> T(i == j), 4))), 4))
    gd = metric(k, xd)
    ∂g = partials.(gd)               # ∂g[μ,ν][σ] = ∂g_{μν}/∂x^σ

    ∇g = zeros(T, 4, 4, 4)
    for σ in 1:4, μ in 1:4, ν in 1:4
        s = ∂g[μ, ν][σ]
        for ρ in 1:4
            s -= Γ[σ, μ, ρ] * g[ρ, ν] # Γ^ρ_{σ μ} * g_{ρ ν}
            s -= Γ[σ, ν, ρ] * g[μ, ρ] # Γ^ρ_{σ ν} * g_{μ ρ}
        end
        ∇g[σ, μ, ν] = s
    end

    return ∇g
end

k  = KerrSchildCoordinates(M=1.0, a=0.5)
u  = SA[0.0, 3.0, 0.0, 1.0]

∇g = covariant_derivative_metric(k, u)
@show size(∇g)
println("máximo |∇g| = ", maximum(abs.(∇g)))
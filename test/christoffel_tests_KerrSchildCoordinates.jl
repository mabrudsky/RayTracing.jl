using Test
using RayTracing
using RayTracing.Kerr: KerrSchildCoordinates, metric, christoffel
using ForwardDiff: Dual, Partials, partials
using StaticArrays

@testset "Símbolos de Christoffel - KerrSchildCoordinates" begin

    function covariant_derivative_metric(k, u::SVector{4,T}) where T
        Γ = christoffel(k, u)
        g = metric(k, u)

        xd = SVector{4, Dual{Nothing,T,4}}(ntuple(i -> Dual{Nothing,T,4}(u[i], Partials(ntuple(j -> T(i == j), 4))), 4))
        gd = metric(k, xd)
        ∂g = partials.(gd)

        ∇g = zeros(T, 4, 4, 4)
        for σ in 1:4, μ in 1:4, ν in 1:4
            s = ∂g[μ, ν][σ]
            for ρ in 1:4
                s -= Γ[σ, μ, ρ] * g[ρ, ν]
                s -= Γ[σ, ν, ρ] * g[μ, ρ]
            end
            ∇g[σ, μ, ν] = s
        end

        return ∇g
    end

    # Caso de prueba
    k = KerrSchildCoordinates(M=1.0, a=0.5)
    u = SA[0.0, 3.0, 0.0, 1.0]

    ∇g = covariant_derivative_metric(k, u)

    # Verifica que el valor máximo de |∇g| sea cero dentro de tolerancia flotante
    @test maximum(abs.(∇g)) ≈ 0 atol=1e-12
end
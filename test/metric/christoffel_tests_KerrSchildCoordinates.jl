using Test
using RayTracing
import ForwardDiff: Dual, Partials, partials
import StaticArrays: SVector, @SVector, @SArray, @MArray

@testset "Christoffel symbols - KerrSchildCoordinates" begin

    function covariant_derivative_metric(k, u::SVector{4,T}) where T
        Γ = christoffel(k, u)
        g = metric(k, u)

        xd = SVector{4, Dual{Nothing,T,4}}(ntuple(i -> Dual{Nothing,T,4}(u[i], Partials(ntuple(j -> T(i == j), 4))), 4))
        gd = metric(k, xd)
        ∂g = partials.(gd)

        ∇g = @MArray zeros(T, 4, 4, 4)
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

    # Test case: Verify that the covariant derivative of the metric is zero for KerrSchildCoordinates
    spacetime = Kerr.KerrSchildCoordinates(M=1.0, a=0.5)
    u = @SArray [0.1, 3.0, 0.05, 1.0]

    ∇g = covariant_derivative_metric(spacetime, u)

    # Verify that the maximum value of |∇g| is zero within floating-point tolerance
    @test maximum(abs.(∇g)) ≈ 0 atol=1e-12  # Equivalent to @test isapprox(maximum(abs.(∇g)), 0; atol=1e-12)
    println("Maximum value of |∇g|: ", maximum(abs.(∇g)))
end

    
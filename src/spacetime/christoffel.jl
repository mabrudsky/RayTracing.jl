abstract type AbstractChristoffel end
struct ChristoffelTag <: AbstractChristoffel end
#struct Christoffel <: AbstractChristoffel end  
#============================================================================================================#


@inline function _Γcomp(ginv::SMatrix{4,4,T}, ∂g, λ::Int, μ::Int, ν::Int) where T
    s = zero(T)
    @inbounds for α in 1:4
        # ∂g[α,ν][μ] is ∂(g_αν)/∂x^μ
        s += ginv[λ, α] * (∂g[α, ν][μ] + ∂g[α, μ][ν] - ∂g[μ, ν][α])
    end
    return T(0.5) * s
end


@inline function _christoffel(k::AbstractSpacetime, position::SVector{N,T}) where {N,T}
    xd = SVector{4, Dual{ChristoffelTag,T,4}}(ntuple(i -> Dual{ChristoffelTag,T,4}(position[i], Partials(ntuple(j -> T(i == j), 4))), 4))
    
    gd   = metric(k, xd)         # SMatrix{4,4,Dual}
    ∂g   = partials.(gd)         # ∂g[μ,ν][σ] = ∂g_{μν}/∂x^σ
    ginv = metric_inverse(k, position)

    return @SArray [_Γcomp(ginv, ∂g, λ, μ, ν) for μ in 1:4, ν in 1:4, λ in 1:4]
end
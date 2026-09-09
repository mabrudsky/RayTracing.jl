abstract type AbstractChristoffel end
struct ChristoffelTag end
#struct Christoffel <: AbstractChristoffel end  
#============================================================================================================#
# @inline function _Γcomp(ginv::SMatrix{4,4,T}, ∂g, λ::Int, μ::Int, ν::Int) where T
    # s  = ginv[λ,1] * (∂g[1,ν][μ] + ∂g[1,μ][ν] - ∂g[μ,ν][1])
    # s += ginv[λ,2] * (∂g[2,ν][μ] + ∂g[2,μ][ν] - ∂g[μ,ν][2])
    # s += ginv[λ,3] * (∂g[3,ν][μ] + ∂g[3,μ][ν] - ∂g[μ,ν][3])
    # s += ginv[λ,4] * (∂g[4,ν][μ] + ∂g[4,μ][ν] - ∂g[μ,ν][4])
    # return T(0.5) * s
# end


@inline function _Γcomp(ginv::SMatrix{4,4,T}, ∂g, λ::Int, μ::Int, ν::Int) where T
    s = zero(T)
    @inbounds for α in 1:4
        # ∂g[α,ν][μ] is ∂(g_αν)/∂x^μ
        s += ginv[λ, α] * (∂g[α, ν][μ] + ∂g[α, μ][ν] - ∂g[μ, ν][α])
    end
    return T(0.5) * s
end


@inline function _christoffel(k::AbstractSpacetime, u) 
    T  = eltype(u)

    x  = SVector{4, T}(u[1], u[2], u[3], u[4])
    xd = SVector{4, Dual{ChristoffelTag,T,4}}(ntuple(i -> Dual{ChristoffelTag,T,4}(x[i], Partials(ntuple(j -> T(i == j), 4))), 4))
    
    gd   = metric(k, xd)         # SMatrix{4,4,Dual}
    ∂g   = partials.(gd)         # ∂g[μ,ν][σ] = ∂g_{μν}/∂x^σ
    ginv = metric_inverse(k, x)

    return @SArray [_Γcomp(ginv, ∂g, λ, μ, ν) for μ in 1:4, ν in 1:4, λ in 1:4]
end
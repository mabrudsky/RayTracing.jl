abstract type AbstractChristoffel end
struct TimeIndependentChristoffel <: AbstractChristoffel end  # static or stationary: ∂g/∂t = 0
struct TimeDependentChristoffel   <: AbstractChristoffel end  # metric with explicit dependence on t
#============================================================================================================#

@inline function _Γcomp(ginv::SMatrix{4,4,T}, ∂g, λ::Int, μ::Int, ν::Int, ::TimeDependentChristoffel) where T
    s  = ginv[λ,1] * (∂g[1,ν][μ] + ∂g[1,μ][ν] - ∂g[μ,ν][1])
    s += ginv[λ,2] * (∂g[2,ν][μ] + ∂g[2,μ][ν] - ∂g[μ,ν][2])
    s += ginv[λ,3] * (∂g[3,ν][μ] + ∂g[3,μ][ν] - ∂g[μ,ν][3])
    s += ginv[λ,4] * (∂g[4,ν][μ] + ∂g[4,μ][ν] - ∂g[μ,ν][4])
    return T(0.5) * s
end


# ∂g[μ,ν][σ-1] con σ=1 (t) igual a cero, sin recalcular la derivada.
@inline function _Γcomp(ginv::SMatrix{4,4,T}, ∂g, λ::Int, μ::Int, ν::Int, ::TimeIndependentChristoffel) where T
    s  = ginv[λ,2] * (∂g[1,ν][μ] + ∂g[1,μ][ν] - ∂g[μ,ν][1])   # σ=2 (x), índice 1 en ∂g (3 comps)
    s += ginv[λ,3] * (∂g[2,ν][μ] + ∂g[2,μ][ν] - ∂g[μ,ν][2])   # σ=3 (y)
    s += ginv[λ,4] * (∂g[3,ν][μ] + ∂g[3,μ][ν] - ∂g[μ,ν][3])   # σ=4 (z)
    # σ=1 (t): ∂g_{μν}/∂t = 0 en todo el término → se omite directamente
    return T(0.5) * s
end

#============================================================================================================#

@inline function _christoffel(::TimeDependentChristoffel, k::AbstractSpacetime, u)
    x = SVector{4}(u[1], u[2], u[3], u[4])
    T = eltype(x)

    xd = SVector{4}(ntuple(i -> Dual{Nothing,T,4}(x[i], Partials(ntuple(j -> T(i == j), 4))), 4))

    gd   = metric(k, xd)         # SMatrix{4,4,Dual}
    ∂g   = partials.(gd)         # ∂g[μ,ν][σ] = ∂g_{μν}/∂x^σ
    ginv = metric_inv(k, x)

    return @SArray [_Γcomp(ginv, ∂g, λ, min(μ,ν), max(μ,ν), TimeDependentChristoffel()) for λ in 1:4, μ in 1:4, ν in 1:4]
end


@inline function _christoffel(::TimeIndependentChristoffel, k::AbstractSpacetime, u) 
    x = SVector{4}(u[1], u[2], u[3], u[4])
    T = eltype(x)

    # Sólo 3 seeds: x, y, z. t queda como Dual sin derivada (partials = 0).
    xd = SVector{4}(Dual{Nothing,T,3}(x[1]),
                    Dual{Nothing,T,3}(x[2], Partials((one(T), zero(T), zero(T)))),
                    Dual{Nothing,T,3}(x[3], Partials((zero(T), one(T), zero(T)))),
                    Dual{Nothing,T,3}(x[4], Partials((zero(T), zero(T), one(T)))),
                   )

    gd   = metric(k, xd)         # SMatrix{4,4,Dual}
    ∂g   = partials.(gd)         # ∂g[μ,ν][σ'] = ∂g_{μν}/∂x^σ', σ'∈{1,2,3}↔{x,y,z}
    ginv = metric_inv(k, x)

    return @SArray [_Γcomp(ginv, ∂g, λ, min(μ,ν), max(μ,ν), TimeIndependentChristoffel()) for λ in 1:4, μ in 1:4, ν in 1:4]
end
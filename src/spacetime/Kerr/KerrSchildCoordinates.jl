# =====================================================================
# Parámetros de Kerr (métrica Kerr-Schild,
# coords cartesianas q = (t, x, y, z))
# =====================================================================
@with_kw struct KerrSchildCoordinates{T<:Union{Float32,Float64}} <: AbstractKerrSpacetime
    M::T
    a::T

    function KerrSchildCoordinates{T}(M, a) where T<:Union{Float32,Float64}
        M >= zero(T) || throw(ArgumentError("M must be non-negative"))
        abs(a) <= M  || throw(ArgumentError("|a| must be smaller than M"))
        new{T}(M, a)
    end
end

KerrSchildCoordinates(M::T, a::T) where T<:Union{Float32,Float64} = KerrSchildCoordinates{T}(M, a)

# ---------------------------------------------------------------------
# auxiliar interno (no exportado): calcula el escalar H y el vector nulo
# covariante l_μ, compartidos entre metric() y metric_inv().
# ---------------------------------------------------------------------
@inline function _H_l(k::KerrSchildCoordinates, x::SVector{4,T}) where T
    _, x1, y1, z1 = x
    
    a = k.a
    M = k.M

    A2   = a^2
    rho2 = x1^2 + y1^2 + z1^2
    r2   = T(0.5) * (rho2 - A2) + sqrt(T(0.25) * (rho2 - A2)^2 + A2 * z1^2)
    r    = sqrt(r2)
    H    = 2*M*r / (r2 + A2 * z1^2 / r2)

    l1 = one(T)
    l2 = (r*x1 + a*y1) / (r2 + A2)
    l3 = (r*y1 - a*x1) / (r2 + A2)
    l4 = z1 / r

    return H, SVector(l1, l2, l3, l4)
end

# =====================================================================
# 1) metric(k, x) -> g_{μν}   (x = (t, x, y, z))
#    Genérica en T: acepta Float32/Float64 o Dual (para christoffel).
# =====================================================================
@inline function metric(k::KerrSchildCoordinates, x::SVector{4,T}) where T
    H, l = _H_l(k, x)

    # términos cruzados calculados una sola vez y reusados abajo,
    # igual que en la versión mutable original (g[2,1]=g[1,2], etc.)
    # -> 10 multiplicaciones garantizadas por construcción, no por
    # confiar en que el optimizador deduplique H*l[i]*l[j] vs H*l[j]*l[i].
    Hl12 = H*l[1]*l[2]; Hl13 = H*l[1]*l[3]; Hl14 = H*l[1]*l[4]
    Hl23 = H*l[2]*l[3]; Hl24 = H*l[2]*l[4]; Hl34 = H*l[3]*l[4]
    Hl11 = H*l[1]*l[1]; Hl22 = H*l[2]*l[2]; Hl33 = H*l[3]*l[3]; Hl44 = H*l[4]*l[4]

    return @SMatrix T[
        -one(T)+Hl11   Hl12          Hl13          Hl14        ;
         Hl12          one(T)+Hl22   Hl23          Hl24        ;
         Hl13           Hl23         one(T)+Hl33   Hl34        ;
         Hl14           Hl24          Hl34         one(T)+Hl44
    ]
end

# =====================================================================
# 2) metric_inv(k, x) -> g^{μν}
#    Kerr-Schild admite inversa analítica exacta (Sherman-Morrison,
#    ya que l es nulo respecto de η): g^{μν} = η^{μν} - H l^μ l^ν,
#    con l^μ = η^{μν} l_ν = (-l_t, l_x, l_y, l_z).
# =====================================================================
@inline function metric_inv(k::KerrSchildCoordinates, x::SVector{4,T}) where T
    H, l = _H_l(k, x)
    lup = SVector(-l[1], l[2], l[3], l[4])

    Hl12 = H*lup[1]*lup[2]; Hl13 = H*lup[1]*lup[3]; Hl14 = H*lup[1]*lup[4]
    Hl23 = H*lup[2]*lup[3]; Hl24 = H*lup[2]*lup[4]; Hl34 = H*lup[3]*lup[4]
    Hl11 = H*lup[1]*lup[1]; Hl22 = H*lup[2]*lup[2]; Hl33 = H*lup[3]*lup[3]; Hl44 = H*lup[4]*lup[4]

    return @SMatrix T[
        -one(T)-Hl11   -Hl12          -Hl13          -Hl14        ;
        -Hl12           one(T)-Hl22   -Hl23          -Hl24        ;
        -Hl13          -Hl23           one(T)-Hl33   -Hl34        ;
        -Hl14          -Hl24          -Hl34           one(T)-Hl44
    ]
end

# =====================================================================
# 3) christoffel(k, u) -> Γ^λ_{μν}
#    Deriva metric() en tiempo de ejecución vía Dual "multi-seed"
#    (isbits -> compatible con kernels CUDA).
# =====================================================================
@inline function _Γcomp(ginv::SMatrix{4,4,T}, ∂g, λ::Int, μ::Int, ν::Int) where T
    s  = ginv[λ,1] * (∂g[1,ν][μ] + ∂g[1,μ][ν] - ∂g[μ,ν][1])
    s += ginv[λ,2] * (∂g[2,ν][μ] + ∂g[2,μ][ν] - ∂g[μ,ν][2])
    s += ginv[λ,3] * (∂g[3,ν][μ] + ∂g[3,μ][ν] - ∂g[μ,ν][3])
    s += ginv[λ,4] * (∂g[4,ν][μ] + ∂g[4,μ][ν] - ∂g[μ,ν][4])
    return T(0.5) * s
end

@inline function christoffel(k::KerrSchildCoordinates, u)
    x = SVector{4}(u[1], u[2], u[3], u[4])
    T = eltype(x)

    xd = SVector{4}(ntuple(
        i -> Dual{Nothing,T,4}(x[i], Partials(ntuple(j -> T(i == j), 4))),
        4))
    gd  = metric(k, xd)         # SMatrix{4,4,Dual}
    ∂g  = partials.(gd)         # ∂g[μ,ν][σ] = ∂g_{μν}/∂x^σ
    ginv = metric_inv(k, x)

    return @SArray [_Γcomp(ginv, ∂g, λ, μ, ν) for λ in 1:4, μ in 1:4, ν in 1:4]
end


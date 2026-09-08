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

# =====================================================================
# auxiliar interno (no exportado): calcula el escalar H y el vector
# nulo l_μ, compartidos entre metric() y metric_inv().
# =====================================================================
@inline function _H_l(k::KerrSchildCoordinates, x::SVector{4,T}) where T
    _, x, y, z = x
    
    a = k.a
    M = k.M

    A2   = a^2
    rho2 = x^2 + y^2 + z^2
    r2   = T(0.5) * (rho2 - A2) + sqrt(T(0.25) * (rho2 - A2)^2 + A2 * z^2)
    r    = sqrt(r2)
    H2   = 2*M*r / (r2 + A2 * z^2 / r2)

    l1 = one(T)
    l2 = (r*x + a*y) / (r2 + A2)
    l3 = (r*y - a*x) / (r2 + A2)
    l4 = z / r

    return H2, SVector(l1, l2, l3, l4)
end
   
# =====================================================================
# 1) metric(k, x) -> g_{μν}   (x = (t, x, y, z))
#    Genérica en T: acepta Float32/Float64 o Dual (para christoffel).
# =====================================================================
@inline function metric(k::KerrSchildCoordinates, x::SVector{4,T}) where T
    H2, l = _H_l(k, x)
    
    g11 = - one(T) + H2 * l[1]*l[1]
    g12 = H2 * l[1]*l[2]
    g13 = H2 * l[1]*l[3]
    g14 = H2 * l[1]*l[4]
    
    g21 = g12
    g22 = one(T) + H2 * l[2]*l[2]
    g23 = H2 * l[2]*l[3]
    g24 = H2 * l[2]*l[4]
    
    g31 = g13
    g32 = g23
    g33 = one(T) + H2 * l[3]*l[3]
    g34 = H2 * l[3]*l[4]
    
    g41 = g14
    g42 = g24
    g43 = g34
    g44 = one(T) + H2 * l[4]*l[4]    

    return SMatrix{4, 4, T, 16}(g11,  g12,  g13,  g14, 
                                g21,  g22,  g23,  g24, 
                                g31,  g32,  g33,  g34, 
                                g41,  g42,  g43,  g44)
end

# =====================================================================
# 2) metric_inv(k, x) -> g^{μν}
#    Kerr-Schild admite inversa analítica exacta (Sherman-Morrison,
#    ya que l es nulo respecto de η): g^{μν} = η^{μν} - H l^μ l^ν,
#    con l^μ = η^{μν} l_ν = (-l_t, l_x, l_y, l_z).
# =====================================================================
@inline function metric_inverse(k::KerrSchildCoordinates, x::SVector{4,T}) where T
    H, l = _H_l(k, x)

    g11 = -one(T) - H2 * l[1] * l[1]
    g12 = - H2 * l[1] * l[2]
    g13 = - H2 * l[1] * l[3]
    g14 = - H2 * l[1] * l[4]
    
    g21 = g12
    g22 = one(T) - H2 * l[2] * l[2]
    g23 = - H2 * l[2] * l[3]
    g24 = - H2 * l[2] * l[4]
    
    g31 = g13
    g32 = g23
    g33 = one(T) - H2 * l[3] * l[3]
    g34 = - H2 * l[3] * l[4]
    
    g41 = g14
    g42 = g24
    g43 = g34
    g44 = one(T) - H2 * l[4] * l[4]

    return SMatrix{4, 4, T, 16}(g11,  g12,  g13,  g14, 
                                g21,  g22,  g23,  g24, 
                                g31,  g32,  g33,  g34, 
                                g41,  g42,  g43,  g44)
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

    xd = SVector{4}(ntuple(i -> Dual{Nothing,T,4}(x[i], Partials(ntuple(j -> T(i == j), 4))),4) )

    gd  = metric(k, xd)         # SMatrix{4,4,Dual}
    ∂g  = partials.(gd)         # ∂g[μ,ν][σ] = ∂g_{μν}/∂x^σ
    ginv = metric_inv(k, x)

    return @SArray [_Γcomp(ginv, ∂g, λ, μ, ν) for λ in 1:4, μ in 1:4, ν in 1:4]
end
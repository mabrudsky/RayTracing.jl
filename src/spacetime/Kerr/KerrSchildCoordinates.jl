# =====================================================================
# Parámetros de Kerr (métrica Kerr-Schild,
# coords cartesianas q = (t, x, y, z))
# =====================================================================
@with_kw struct KerrSchildCoordinates{T<:Union{Float32,Float64}} <: AbstractKerrSpacetime
    M::T
    a::T
    @assert M >= zero(T) "M must be non-negative"
    @assert abs(a) <= M "|a| must be smaller than M"
end

KerrSchildCoordinates(M::Real, a::Real) = throw(ArgumentError("Spacetime parameters must be Float32 or Float64 (both the same type); got M::$(typeof(M)), a::$(typeof(a))."))

# =====================================================================
# auxiliar interno (no exportado): calcula el escalar H y el vector
# nulo l_μ, compartidos entre metric() y metric_inv().
# =====================================================================
@inline function _H_l(k::KerrSchildCoordinates, point::SVector{4,T}) where T
    t, x, y, z = point
    
    a = k.a
    M = k.M

    A2   = a^2
    rho2 = x^2 + y^2 + z^2
    r2   = T(0.5) * (rho2 - A2) + sqrt(T(0.25) * (rho2 - A2)^2 + A2 * z^2)
    r    = sqrt(r2)
    _2H  = 2*M*r / (r2 + A2 * z^2 / r2)

    l1 = one(T)
    l2 = (r*x + a*y) / (r2 + A2)
    l3 = (r*y - a*x) / (r2 + A2)
    l4 = z / r

    return _2H, SVector(l1, l2, l3, l4)
end
   
# =====================================================================
# 1) metric(k, x) -> g_{μν}   (x = (t, x, y, z))
#    Genérica en T: acepta Float32/Float64 o Dual (para christoffel).
# =====================================================================
@inline function metric(k::KerrSchildCoordinates, point::SVector{4,T}) where T
    _2H, l = _H_l(k, point)
    
    g11 = - one(T) + _2H * l[1]*l[1]
    g12 = _2H * l[1]*l[2]
    g13 = _2H * l[1]*l[3]
    g14 = _2H * l[1]*l[4]
    
    g21 = g12
    g22 = one(T) + _2H * l[2]*l[2]
    g23 = _2H * l[2]*l[3]
    g24 = _2H * l[2]*l[4]
    
    g31 = g13
    g32 = g23
    g33 = one(T) + _2H * l[3]*l[3]
    g34 = _2H * l[3]*l[4]
    
    g41 = g14
    g42 = g24
    g43 = g34
    g44 = one(T) + _2H * l[4]*l[4]    

    return SMatrix{4, 4, T, 16}(g11, g21, g31, g41, 
                                g12, g22, g32, g42, 
                                g13, g23, g33, g43, 
                                g14, g24, g34, g44)
end

# ===============================================================#
# 2) metric_inverse(k, x) -> g^{μν} = η^{μν} - 2H l^μ l^ν        #
#  with l^μ = η^{μν} l_ν = (-1, l^i) (η^{μν} = diag(-1, 1, 1, 1))#
# ===============================================================#
@inline function metric_inverse(k::KerrSchildCoordinates, point::SVector{4,T}) where T
    _2H, l = _H_l(k, point)

    g11 = -one(T) - _2H * l[1] * l[1]
    g12 =  _2H * l[1] * l[2]
    g13 =  _2H * l[1] * l[3]
    g14 =  _2H * l[1] * l[4]
    
    g21 = g12
    g22 = one(T) - _2H * l[2] * l[2]
    g23 = - _2H * l[2] * l[3]
    g24 = - _2H * l[2] * l[4]
    
    g31 = g13
    g32 = g23
    g33 = one(T) - _2H * l[3] * l[3]
    g34 = - _2H * l[3] * l[4]
    
    g41 = g14
    g42 = g24
    g43 = g34
    g44 = one(T) - _2H * l[4] * l[4]


    return SMatrix{4, 4, T, 16}(g11, g21, g31, g41, 
                                g12, g22, g32, g42, 
                                g13, g23, g33, g43, 
                                g14, g24, g34, g44)
end

#========================================================#
# 3) christoffel(k, u) -> Γ^λ_{μν}                       #  
#    Differentiates metric() at runtime via "multi-seed" #
#    Duals (isbits -> compatible with CUDA kernels).     #  
#========================================================#
christoffel(spacetime::KerrSchildCoordinates, position::SVector{4,T}) where T = _christoffel(spacetime, position)
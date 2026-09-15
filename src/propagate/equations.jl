#=
function condition(u,p)
    @inbounds R_NS, d = p[3], p[4]

        Rmin = R_NS       # Inner radius cut
        Rmax = d*1.1f0    # Outer radius cut
       Rmin2 = Rmin*Rmin
       Rmax2 = Rmax*Rmax
          r2 = u[2] * u[2] + u[3] * u[3] + u[4] * u[4] 
          
   return (Rmin2 <= r2) * (Rmax2 >= r2) * 1.0f0
end
=#
function geodesic_equations(u::SVector{N,T}, p, t) where {N,T}
    position = @SVector T[u[1], u[2], u[3], u[4]]
    k = @SVector T[u[5], u[6], u[7], u[8]]  

    kkT = k * k'   
     
    Γ = christoffel(p.spacetime, position) 
    a = @SVector T[-sum(Γ[:, :, λ] .* kkT) for λ in 1:4]    # Geodesic equation kᵃ∇ₐkᵇ=0 where k is the momentum of the photon.

    #cond = 1.0 # condition(u,p)
    
    #du = cond .* vcat(k, a)
    du = vcat(k, a)
    return du
end 

function parallel_transport_equations(u::SVector{N,T}, p, t) where {N,T}
    position = @SVector T[u[1], u[2], u[3], u[4]]
    k = @SVector T[u[5], u[6], u[7], u[8]]  

    ex = @SVector T[u[9], u[10], u[11], u[12]]
    ey = @SVector T[u[13], u[14], u[15], u[16]]

    kkT = k * k'    
    k_ex = k * ex'
    k_ey = k * ey'
    
    Γ = christoffel(p.spacetime, position)
    
    a   = @SVector T[-sum(Γ[:, :, λ] .* kkT) for λ in 1:4]   # Geodesic equation kᵃ∇ₐkᵇ=0 where k is the momentum of the photon. 
    dex = @SVector T[-sum(Γ[:, :, λ] .* k_ex) for λ in 1:4]  # Parallel transport kᵃ∇ₐexᵇ=0 where ex is the rectangular coordinates of the image plane.            
    dey = @SVector T[-sum(Γ[:, :, λ] .* k_ey) for λ in 1:4]  # Parallel transport kᵃ∇ₐeyᵇ=0 where ey is the rectangular coordinates of the image plane. 

    du = vcat(k, a, dex, dey)
    return du
end
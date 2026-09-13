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
    momentum = @SVector T[u[5], u[6], u[7], u[8]]

    mmT = momentum * momentum'   
     
    Γ = christoffel(p.spacetime, position) 
    a = @SVector T[-sum(Γ[:, :, λ] .* mmT) for λ in 1:4]

    #cond = 1.0 # condition(u,p)
    
    #du1 = cond * u[5] 
    #du2 = cond * u[6]
    #du3 = cond * u[7]
    #du4 = cond * u[8]
    #du5 = cond * a[1]
    #du6 = cond * a[2]
    #du7 = cond * a[3]
    #du8 = cond * a[4]
    
    du = vcat(momentum, a)
    return du
end 

function parallel_transport_equations(u::SVector{N,T}, p, t) where {N,T}
    position = @SVector T[u[1], u[2], u[3], u[4]]
    momentum = @SVector T[u[5], u[6], u[7], u[8]]

    ex = @SVector T[u[9], u[10], u[11], u[12]]
    ey = @SVector T[u[13], u[14], u[15], u[16]]

    mmT = momentum * momentum'    
    momentum_ex = momentum * ex'
    momentum_ey = momentum * ey'
    
    Γ = christoffel(p.spacetime, position)
    
    a   = @SVector T[-sum(Γ[:, :, λ] .* mmT) for λ in 1:4]
    dex = @SVector T[-sum(Γ[:, :, λ] .* momentum_ex) for λ in 1:4]
    dey = @SVector T[-sum(Γ[:, :, λ] .* momentum_ey) for λ in 1:4]

    du = vcat(momentum, a, dex, dey)
    return du
end
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
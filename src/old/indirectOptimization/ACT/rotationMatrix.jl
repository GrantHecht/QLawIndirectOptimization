
crossProductMatrix(a) = SMatrix{3,3}(0.0, a[3], -a[2], -a[3], 0.0, a[1], a[2], -a[1], 0.0)

function rotationMatrix(cart, ps)
    # Compute dv
    r       = norm(cart[1:3])
    dv      = SVector(-ps.μ*cart[1] / r^3,
                      -ps.μ*cart[2] / r^3,
                      -ps.μ*cart[3] / r^3)

    # Compute required cross product matricies
    @views rx = crossProductMatrix(cart[1:3])
    @views vx = crossProductMatrix(cart[4:6])

    # Compute drhatdr
    rinv      = 1.0 / r
    rinv3     = rinv^3
    drhatdr   = SMatrix{3,3}(rinv - cart[1]^2*rinv3,
                             -cart[1]*cart[2]*rinv3,
                             -cart[1]*cart[3]*rinv3,
                             -cart[1]*cart[2]*rinv3,
                             rinv - cart[2]^2*rinv3,
                             -cart[2]*cart[3]*rinv3,
                             -cart[1]*cart[3]*rinv3,
                             -cart[2]*cart[3]*rinv3,
                             rinv - cart[3]^2*rinv3)

    # Compute dhhatdh
    hvec      = SVector(cart[2]*cart[6] - cart[3]*cart[5],
                        cart[3]*cart[4] - cart[1]*cart[6],
                        cart[1]*cart[5] - cart[2]*cart[4])
    hx        = crossProductMatrix(hvec)
    h         = hvec / norm(hvec)
    hinv      = 1.0 / h
    hinv3     = hinv^3
    dhhatdh   = SMatrix{3,3}(hinv - hvec[1]^2*hinv3,
                             -hvec[1]*hvec[2]*hinv3,
                             -hvec[1]*hvec[3]*hinv3,
                             -hvec[1]*hvec[2]*hinv3,
                             hinv - hvec[2]^2*hinv3,
                             -hvec[2]*hvec[3]*hinv3,
                             -hvec[1]*hvec[3]*hinv3,
                             -hvec[2]*hvec[3]*hinv3,
                             hinv - hvec[3]^2*hinv3)

    # Compute dthatdt
    tvec      = SVector(hvec[2]*cart[3] - hvec[3]*cart[2],
                        hvec[3]*cart[1] - hvec[1]*cart[3],
                        hvec[1]*cart[2] - hvec[2]*cart[1])
    t         = tvec / norm(tvec)
    tinv      = 1.0 / t
    tinv3     = tinv^3
    dthatdt   = SMatrix{3,3}(tinv - tvec[1]^2*tinv3,
                             -tvec[1]*tvec[2]*tinv3,
                             -tvec[1]*tvec[3]*tinv3,
                             -tvec[1]*tvec[2]*tinv3,
                             tinv - tvec[2]^2*tinv3,
                             -tvec[2]*tvec[3]*tinv3,
                             -tvec[1]*tvec[3]*tinv3,
                             -tvec[2]*tvec[3]*tinv3,
                             tinv - tvec[3]^2*tinv3)   

    # Compute unit vectors
    u1        = SVector(cart[1] * rinv, cart[2] * rinv, cart[3] * rinv)
    u2        = SVector(tvec[1] * tinv, tvec[2] * tinv, tvec[3] * tinv)
    u3        = SVector(hvec[1] * hinv, hvec[2] * hinv, hvec[3] * hinv)

    # Compute unit vector time derivatives
    du1       = drhatdr*cart[4:6]
    du2       = dthatdt*((rx*vx + hx)*cart[4:6] - rx*rx*dv)
    du3       = dhhatdh*(rx*dv - vx*cart[4:6])

    # Return matricies
    return (hcat(u1,u2,u3), hcat(du1, du2, du3))
end
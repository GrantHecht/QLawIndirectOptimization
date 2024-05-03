
function get_third_body_perturbation(cart, t, ps::qLawParams)
    return AstroEOMs.computeThirdBodyPerterbations(cart, t, ps.meePs)
end

function get_j2_perturbation(mee, ps::qLawParams)
    return AstroEOMs.computeJ2Perterbations(mee, ps.meePs)
end

function find_c2c3(ψ)
    if ψ > 1e-6
        c2 = (1.0 - cos(sqrt(ψ))) / ψ
        c3 = (sqrt(ψ) - sin(sqrt(ψ))) / sqrt(ψ^3)
    else
        if ψ < -1e-6
            c2 = (1.0 - cosh(sqrt(-ψ))) / ψ
            c3 = (sinh(sqrt(-ψ)) - sqrt(-ψ)) / sqrt((-ψ)^3)
        else
            c2 = 0.5
            c3 = 1.0 / 6.0
        end
    end
    return c2, c3
end

# See valado for details
function kepler(r0_vec, v0_vec, Δt, μ; tol = 1e-6, max_iters = 100)
    # Compute commonly used terms
    r0 = norm(r0_vec)
    v0 = norm(v0_vec)
    dotr0v0 = dot(r0_vec, v0_vec)

    r0inv = 1.0 / r0
    v0squared = v0^2

    sqrtmu = sqrt(μ)
    sqrtmuinv = 1.0 / sqrtmu

    # Find initial guess for χ
    ξ  = 0.5*v0squared - μ * r0inv
    α  = -v0squared / μ + 2.0 * r0inv
    if α >= 1e-6 # Circule or ellipse
        χ₀ = sqrtmu * Δt * α
    elseif abs(α) < 1e-6 # Parabola
        h_vec = cross(r0_vec, v0_vec)
        p  = norm(h_vec)^2 / μ
        s  = 0.5*acot(3.0*sqrt(μ / p)*Δt)
        w  = atan(tan(s)^(1.0/3.0))
        χ₀ = 2.0*sqrt(p)*cot(2.0*w)
    else # Hyperbola
        a   = 1.0 / α
        lnt = -2.0*μ*α*Δt / (dotr0v0 + sign(Δt)*sqrt(-μ*a)*(1.0 - r0*α))
        χ₀  = sign(Δt)*sqrt(-a)*log(lnt)
    end

    # Begin iteration
    χₙ   = χ₀
    i    = 0
    done = false
    local r
    local c2
    local c3
    local ψ
    while !done && i <= max_iters
        # Requirements
        χₙsquared = χₙ^2
        ψ = χₙsquared * α

        # Compute c2 and c3
        c2, c3 = find_c2c3(ψ)

        # Compute r
        r = χₙsquared*c2 + dotr0v0*χₙ*(1.0 - ψ*c3)*sqrtmuinv + r0*(1.0 - ψ*c2)

        # Update χ
        update = (sqrtmu*Δt - χₙ*χₙsquared*c3 - dotr0v0*χₙsquared*c2*sqrtmuinv - r0*χₙ*(1.0 - ψ*c3)) / r
        χₙ₊₁ = χₙ + update
        if χₙ₊₁ < 0.0
            χₙ₊₁ = 0.5*(χₙ + χₙ₊₁)
        end
        χₙ₋₁ = χₙ
        χₙ = χₙ₊₁

        # Check convergence
        if abs(χₙ - χₙ₋₁) < tol
            done = true
        end

        # Increment counter
        i += 1
    end

    # Compute f and g functions
    χₙsquared = χₙ^2
    f = 1.0 - χₙsquared*c2*r0inv
    g = Δt - χₙ*χₙsquared*c3*sqrtmuinv
    dg = 1.0 - χₙsquared*c2 / r
    df = sqrtmu*χₙ*(ψ*c3 - 1.0)*r0inv / r

    # Compute r and v
    r_vec = f*r0_vec + g*v0_vec
    v_vec = df*r0_vec + dg*v0_vec
    flag  = abs(f*dg - df*g - 1.0) < 1e-10
    return r_vec, v_vec, flag
end
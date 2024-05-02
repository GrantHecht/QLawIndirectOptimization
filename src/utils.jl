
function get_third_body_perturbation(cart, t, ps::qLawParams)
    return AstroEOMs.computeThirdBodyPerterbations(cart, t, ps.meePs)
end

function get_j2_perturbation(mee, ps::qLawParams)
    return AstroEOMs.computeJ2Perterbations(mee, ps.meePs)
end
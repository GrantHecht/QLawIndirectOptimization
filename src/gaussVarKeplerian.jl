
function gaussVarKeplerian(sma, e, inc, ran, ape, tru, ps)
    # Compute requirements
    cta     = cos(tru)
    sta     = sin(tru)
    s_ta_ω  = sin(tru + ape)

    p       = sma*(1.0 - e*e)
    h       = sqrt(ps.μ*p)
    r       = p / (1.0 + e*cta)

    hinv    = 1.0 / h
    einv    = 1.0 / e

    # Compute variational equations (without ta)
    A = SA[
        2*sma*sma*hinv*e*sta    2*sma*sma*hinv*p / r        0.0;
        hinv*p*sta              hinv*((p + r)*cta + r*e)    0.0;                   
        0.0                     0.0                         r*hinv*cos(tru + ape);
        0.0                     0.0                         r*hinv*s_ta_ω;
        -p*hinv*cta*einv        hinv*(p + r)*sta*einv       -(r*hinv*s_ta_ω*cot(inc));
        p*hinv*cta*einv         -(p + r)*hinv*sta*einv      0.0;
    ]
    B = SA[0.0, 0.0, 0.0, 0.0, 0.0, h / (r*r)]

    return (A,B)
end

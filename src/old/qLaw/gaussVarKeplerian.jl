
function gaussVarKeplerian(sma, e, inc, ran, ape, tru, ps)
    # Compute requirements
    p       = sma*(1.0 - e*e)
    h       = sqrt(ps.μ*p)
    r       = p / (1.0 + e*cos(tru))
    hinv    = 1.0 / h

    # Compute variational equations (without ta)
    A       = SMatrix{6,3}(2*sma*sma*hinv*e*sin(tru),       # dadfr
                           hinv*p*sin(tru),                 # dedfr
                           0.0,                             # didfr
                           0.0,                             # dΩdfr
                           -p*hinv*cos(tru) / e,            # dΩdfr
                           p*hinv*cos(tru) / e,             # dθdfr
                           2*sma*sma*hinv*p / r,            # dadfθ
                           hinv*((p + r)*cos(tru) + r*e),   # dedfθ
                           0.0,                             # didfθ 
                           0.0,                             # dΩdfθ
                           hinv*(p + r)*sin(tru) / e,       # dωdfθ
                           -(p + r)*hinv*sin(tru) / e,      # dθdfθ
                           0.0,                             # dadfh
                           0.0,                             # dedfh
                           r*hinv*cos(tru + ape),           # didfh
                           r*hinv*sin(tru + ape) / sin(inc),# dΩdfh
                           -(r*hinv*sin(tru + ape)*cos(inc) / sin(inc)),
                           0.0)

    B       = SVector(0.0, 0.0, 0.0, 0.0, 0.0, h / (r*r))

    return (A,B)
end

function qLaw(ps::qLawParams)
    # Compute p and true long.
    L0      = ps.oei[4] + ps.oei[5] + ps.oei[6]
    p       = ps.oei[1]*(1.0 - ps.oei[2]*ps.oei[2])

    # Prepare for simulation
    Lf              = L0 + 2.0*pi*ps.maxRevs
    L_bin           = 0.0
    nseg            = ceil(Int, (Lf - L0) / ps.step)
    X_next          = zeros(7)
    X_curr          = [p, ps.oei[2], ps.oei[3], ps.oei[4], 
                        ps.oei[5], ps.m0, ps.t0]
    Xf              = zeros(7)
    sma_curr        = ps.oei[1]
    L_curr          = L0
    percentComplete = 0.0
    accumulate      = 1.0
    converged       = false
    writebin        = fill(NaN, nseg + 1, 8)
    writerow        = 1

    # Write initial row to storage
    writebin[writerow, 1:7] .= X_curr
    writebin[writerow, 8]   = L_curr
    writerow += 1

    # Begin fixed step integration
    print("Simulation Progress: \n")
    for j = 1:nseg
        # If we run out of mass, break from sim loop 
        if X_curr[6] <= 0.0
            print("Out of mass.\n")
            break
        end

        # Take an integration step

    end
end
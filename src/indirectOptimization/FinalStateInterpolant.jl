
struct FinalStateInterpolant{TPF,FPF,GPF,HPF,KPF,TFPF}
    pf::TPF
    ff::FPF
    gf::GPF
    hf::HPF
    kf::KPF
    tf::TFPF
end

function FinalStateInterpolant(meeAtSteps, startingIdx)
    # Compute number of steps
    n   = size(meeAtSteps,1) - startingIdx + 1

    # Create independent variable vector
    κ   = range(0.0, 1.0, length = n)

    # Create interpolants
    pf  = linear_interpolation(κ, meeAtSteps[startingIdx:end, 1])
    ff  = linear_interpolation(κ, meeAtSteps[startingIdx:end, 2])
    gf  = linear_interpolation(κ, meeAtSteps[startingIdx:end, 3])
    hf  = linear_interpolation(κ, meeAtSteps[startingIdx:end, 4])
    kf  = linear_interpolation(κ, meeAtSteps[startingIdx:end, 5])
    tf  = linear_interpolation(κ, meeAtSteps[startingIdx:end, 8])

    # Contruct and return object
    return FinalStateInterpolant(pf,ff,gf,hf,kf,tf)
end

# Get desired final state and time for give κ
function getFinalStateAndTime(interp::FinalStateInterpolant, κ)
    return SVector(interp.pf(κ),
                   interp.ff(κ),
                   interp.gf(κ),
                   interp.hf(κ),
                   interp.kf(κ),
                   interp.tf(κ))
end


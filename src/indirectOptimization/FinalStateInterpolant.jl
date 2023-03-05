
# struct FinalStateInterpolant{TPF,FPF,GPF,HPF,KPF,TFPF}
#     pf::TPF
#     ff::FPF
#     gf::GPF
#     hf::HPF
#     kf::KPF
#     tf::TFPF
# end

# function FinalStateInterpolant(meeAtSteps, startingIdx)
#     # Compute number of steps
#     n   = size(meeAtSteps,1) - startingIdx + 1

#     # Create independent variable vector
#     κ   = range(0.0, 1.0, length = n)

#     # Create interpolants
#     pf  = linear_interpolation(κ, meeAtSteps[startingIdx:end, 1])
#     ff  = linear_interpolation(κ, meeAtSteps[startingIdx:end, 2])
#     gf  = linear_interpolation(κ, meeAtSteps[startingIdx:end, 3])
#     hf  = linear_interpolation(κ, meeAtSteps[startingIdx:end, 4])
#     kf  = linear_interpolation(κ, meeAtSteps[startingIdx:end, 5])
#     tf  = linear_interpolation(κ, meeAtSteps[startingIdx:end, 8])

#     # Contruct and return object
#     return FinalStateInterpolant(pf,ff,gf,hf,kf,tf)
# end

# # Get desired final state and time for give κ
# function getFinalStateAndTime(interp::FinalStateInterpolant, κ)
#     return SVector(interp.pf(κ),
#                    interp.ff(κ),
#                    interp.gf(κ),
#                    interp.hf(κ),
#                    interp.kf(κ),
#                    interp.tf(κ))
# end

struct FinalStateInterpolant{KT}
    # Independant variables
    κs::KT

    # Matrix of values to interpolate
    vals::Matrix{Float64}
end

function FinalStateInterpolant(meeAtSteps, startingIdx)
    # Compute number of steps
    n   = size(meeAtSteps,1) - startingIdx + 1

    # Create independant variable
    κs  = range(0.0, 1.0, length = n)

    # Prune data
    vals    = meeAtSteps[startingIdx:end,:]

    # Construct and return
    return FinalStateInterpolant(κs, vals)
end

# Get desired final state and time for given κ
function getFinalStateAndTime(interp::FinalStateInterpolant, κ)
    if κ < 0.0 # Linear extrapolation
        κ1  = interp.κs[1]
        κ2  = interp.κs[2]
        τ1  = (κ2 - κ) / (κ2 - κ1)
        τ2  = (κ - κ1) / (κ2 - κ1)
        return SVector(interp.vals[1,1]*τ1 + interp.vals[2,1]*τ2,
                       interp.vals[1,2]*τ1 + interp.vals[2,2]*τ2,
                       interp.vals[1,3]*τ1 + interp.vals[2,3]*τ2,
                       interp.vals[1,4]*τ1 + interp.vals[2,4]*τ2,
                       interp.vals[1,5]*τ1 + interp.vals[2,5]*τ2,
                       interp.vals[1,8]*τ1 + interp.vals[2,8]*τ2)
    elseif κ > 1.0 # Linear extrapolation
        κ1  = interp.κs[end - 1]
        κ2  = interp.κs[end]
        τ1  = (κ2 - κ) / (κ2 - κ1)
        τ2  = (κ - κ1) / (κ2 - κ1)
        return SVector(interp.vals[end - 1,1]*τ1 + interp.vals[end,1]*τ2,
                       interp.vals[end - 1,2]*τ1 + interp.vals[end,2]*τ2,
                       interp.vals[end - 1,3]*τ1 + interp.vals[end,3]*τ2,
                       interp.vals[end - 1,4]*τ1 + interp.vals[end,4]*τ2,
                       interp.vals[end - 1,5]*τ1 + interp.vals[end,5]*τ2,
                       interp.vals[end - 1,8]*τ1 + interp.vals[end,8]*τ2)
    end
    @inbounds for i in 1:length(interp.κs) - 1
        if κ >= interp.κs[i] && κ <= interp.κs[i + 1]
            κ1  = interp.κs[i]
            κ2  = interp.κs[i + 1]
            τ1  = (κ2 - κ)/(κ2 - κ1) 
            τ2  = (κ - κ1)/(κ2 - κ1)
            return SVector(interp.vals[i,1]*τ1 + interp.vals[i + 1,1]*τ2,
                           interp.vals[i,2]*τ1 + interp.vals[i + 1,2]*τ2,
                           interp.vals[i,3]*τ1 + interp.vals[i + 1,3]*τ2,
                           interp.vals[i,4]*τ1 + interp.vals[i + 1,4]*τ2,
                           interp.vals[i,5]*τ1 + interp.vals[i + 1,5]*τ2,
                           interp.vals[i,8]*τ1 + interp.vals[i + 1,8]*τ2)
            break
        end
    end
end
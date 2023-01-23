using DrWatson
@quickactivate "QLawIndirectOptimization"

# This is not working quite as well as I hoped....
# Think we're going to have to wait for Symbolics to improve

using Symbolics

# Flag to set if simplifying final expressions
simplifyExp = false

# Define our variables
@variables μ a e i ω Ω θ f α β Wp Wa We Wi Wω WΩ rpmin k at et it ωt Ωt fr fθ fh;

# Define some differential operators
Da = Differential(a)
De = Differential(e)
Di = Differential(i)
Dω = Differential(ω)
DΩ = Differential(Ω)
Dfr = Differential(fr)
Dfθ = Differential(fθ)
Dfh = Differential(fh)


# Set some constants which can be changed for tuning
m       = 3
n       = 4
rt      = 2
b       = 0.01

# Define some required expressions
p    = a*(1.0 - e)
r    = p / (1.0 + e*cos(θ))
h    = sqrt(p*μ)
rp   = a*(1.0 - e)
P    = exp(k*(1.0 - rp/rpmin))
Sa   = (1 + ((a - at)/(m*at))^n)^(1/rt)
Se   = 1
Si   = 1
Sω   = 1
SΩ   = 1
dat  = a - at
det  = e - et
dit  = i - it
dωt  = acos(cos(ω - ωt))
dΩt  = acos(cos(Ω - Ωt))

# Gauss's form of the variational equations
dΩ = r*sin(θ + ω)*fh/(h*sin(i))
di = r*cos(θ + ω)*fh/h
dω = (1.0/(e*h))*(-p*cos(θ)*fr + (p+r)*sin(θ)*fθ) - r*sin(θ+ω)*cos(i)*fh/(h*sin(i))
da = 2*a^2*(e*sin(θ)*fr + p*fθ/r)/h
de = (1.0/h)*(p*sin(θ)*fr + ((p+r)*cos(θ) + r*e)*fθ)
dθ = h/r^2 + (p*cos(θ)*fr - (p+r)*sin(θ)*fθ)/(e*h)

# Maximum rate of changes
daxx  = 2*f*sqrt(a^3*(1 + e)/(μ*(1 - e)))
dexx  = 2*p*f/h
dixx  = p*f/(h*(sqrt(1 - e^2*sin(ω)^2) - e*abs(cos(ω))))
dΩxx  = p*f/(h*sin(i)*(sqrt(1 - e^2*cos(ω)^2) - e*abs(sin(ω))))

cθxx  = ((1 - e^2)/(2*e^3) + sqrt((1/4)*((1 - e^2)/e^3)^2 + 1/27))^(1/3) - 
        (-(1-e^2)/(2*e^3) + sqrt((1/4)*((1-e^2)/e^3)^2 + 1/27))^(1/3) - 1/e
rxx   = p/(1 + e*cθxx)
dωxxi = f*sqrt(p^2*cθxx^2 + (p + rxx)^2*(1 - cθxx^2))/(e*h)
dωxxo = dΩxx*abs(cos(i))
dωxx  = (dωxxi + b*dωxxo) / (1 + b) 

# Define the QLaw function
Q = (1 + Wp*P)*(Wa*Sa*(dat / daxx)^2 + We*Se*(det / dexx)^2 + Wi*Si*(dit / dixx)^2 + Wω*Sω*(dωt / dωxx)^2 + WΩ*SΩ*(dΩt / dΩxx)^2)

# Compute the D terms
D1 = expand_derivatives(Da(Q)*Dfθ(da) + De(Q)*Dfθ(de) + Di(Q)*Dfθ(di) + Dω(Q)*Dfθ(dω) + DΩ(Q)*Dfθ(dΩ))
D2 = expand_derivatives(Da(Q)*Dfr(da) + De(Q)*Dfr(de) + Di(Q)*Dfr(di) + Dω(Q)*Dfr(dω) + DΩ(Q)*Dfr(dΩ))
D3 = expand_derivatives(Da(Q)*Dfh(da) + De(Q)*Dfh(de) + Di(Q)*Dfh(di) + Dω(Q)*Dfh(dω) + DΩ(Q)*Dfh(dΩ))

# Compute the QLaw thrust directions
αs = simplifyExp ? simplify(atan(-D2,-D1); threaded = true) : atan(-D2,-D1)
βs = simplifyExp ? simplify(atan(-D3 / sqrt(D1^2 + D2^2)); threaded = true) : atan(-D3 / sqrt(D1^2 + D2^2))
as = [αs, βs]

# Compute Qdot_n
dQ = expand_derivatives(Da(Q)*da + De(Q)*de + Di(Q)*di + Dω(Q)*dω + DΩ(Q)*dΩ)
dQ = substitute(dQ, Dict([fr => f*cos(β)*sin(α), fθ => f*cos(β)*cos(α), fh => f*sin(β)]))
dQ = simplifyExp ? simplify(dQ) : dQ

# Build Julia functions
asFuncName      = "qLawThrustAngles"
vars            = [a,e,i,ω,Ω,θ,at,et,it,ωt,Ωt,μ,f,Wp,Wa,We,Wi,Wω,WΩ,rpmin,k]
asFunc          = Symbolics._build_function(Symbolics.JuliaTarget(), as, vars...;
                        parallel = Symbolics.SerialForm(),
                        linenumbers = false,
                        force_SA = true)
                
qnFuncName      = "Qn"
vars            = [a,e,i,ω,Ω,θ,at,et,it,ωt,Ωt,μ,f,α,β,Wp,Wa,We,Wi,Wω,WΩ,rpmin,k]
qnFunc          = Symbolics._build_function(Symbolics.JuliaTarget(), dQ, vars...;
                        parallel = Symbolics.SerialForm(),
                        linenumbers = false)

# Write functions to source directory
asFuncFileName = srcdir(asFuncName * ".jl")
touch(asFuncFileName)
asFuncFile = open(asFuncFileName, "w")
write(asFuncFile, string(asFunc[1]))
close(asFuncFile)

qnFuncFileName = srcdir(qnFuncName * ".jl")
touch(qnFuncFileName)
qnFuncFile = open(qnFuncFileName, "w")
write(qnFuncFile, string(qnFunc))
close(qnFuncFile)
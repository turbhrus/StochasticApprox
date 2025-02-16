
################################################################################
# Clenshaw stage
struct SGClenshawStage <: StochasticGate
    mul::SGMulBipolar
    mac::SG_SDP
end
SGClenshawStage(coeffs=[1,2,-1]) = SGClenshawStage(SGMulBipolar(), SG_SDP(coeffs))
function proc(cs::SGClenshawStage, x::SCBit, aₙ::SCBit, bₙ₊₁::SCBit, bₙ₊₂::SCBit)
    xbₙ₊₁ = proc(cs.mul, x, bₙ₊₁)       # x*bₙ₊₁
    return proc(cs.mac, aₙ, xbₙ₊₁, bₙ₊₂)
end


###############################################################################
# Clenshaw alg implementation
struct SGClenshaw <: StochasticGate
    n::Int                         # number of terms
    coeff::Vector{Float64}         # in unipolar [0,1] domain
    b::Vector{SCBit}
    stage::Vector{SGClenshawStage}
    noise::SCNoise
end
# constructor
function SGClenshaw(coeff::Vector{Float64}, ϵ=0.0)
    n = length(coeff)
    b = BitVector(undef, n+2)
    stage = Vector{SGClenshawStage}(undef, n)
    cv(c, a) = [c[1]*(a!=0.5); c[2:end]]   # if a==0.5 (x==0) then ignore it
    stage[1]   = SGClenshawStage(cv([1,1,-1], coeff[1]))
    stage[n-1] = SGClenshawStage(cv([1,2,0],  coeff[n-1]))
    stage[n]   = SGClenshawStage(cv([1,0,0],  coeff[n]))
    for i in 2:n-2
        stage[i] = SGClenshawStage(cv([1,2,-1],coeff[i]))
    end
    noise = SCNoise(ϵ/n)
    return SGClenshaw(n, coeff, b, stage, noise)
end
Base.broadcastable(g::SGClenshaw) = Ref(g)
# processing
function proc(c::SGClenshaw, x::SCBit)
    ismissing(x) && return missing
    for i in c.n:-1:1
        a = bernoulli_samp(c.coeff[i])
        c.b[i] = proc(c.stage[i], c.noise, x, a, c.b[i+1], c.b[i+2])
    end
    return c.b[1]
end




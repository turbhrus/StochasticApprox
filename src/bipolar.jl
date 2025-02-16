const OPBP_BUFLEN = 2  

################################################################################
# Identity: x = x
struct SGIdentity <: StochasticGate
    info::String          # description of the function implemented by this gate
    nin::Int
    nout::Int
    memoryless::Bool      # output sequence uncorrelated with itself?
    SGIdentity() = new("x", 1, 1, true)
end
Base.broadcastable(g::SGIdentity) = Ref(g)
proc(g::SGIdentity, b::SCBit) = b

################################################################################
# Bipolar multiplier: x₁x₂/L
struct SGMulBipolar <: StochasticGate
    info::String          # description of the function implemented by this gate
    nin::Int
    nout::Int
    memoryless::Bool      # output sequence uncorrelated with itself?
    SGMulBipolar() = new("x₁x₂/L", 2, 1, true)
end
Base.broadcastable(g::SGMulBipolar) = Ref(g)
proc(g::SGMulBipolar, b1::SCBit, b2::SCBit) = !(b1 ⊻ b2)


################################################################################
# Stochastic dot-product (SDP)
# Requires an interleaver at the output of minimum length `M`
mutable struct SG_SDP <: StochasticGate
    M::Int
    cs::Vector{Int}       # scale factors
    m::Int
    σ₀::Int
    σ::Int
    π₁::Int               # number of "1"s to output in buffer cycle
    globalctr::UInt
    t::Int                # accumulator
    info::String          # description of the function implemented by this gate
    nin::Int              # number of input bits
    nout::Int             # number of output bits
    memoryless::Bool      # output sequence uncorrelated with itself?
end
function SG_SDP(coeffs::Vector{<:Integer}, M=OPBP_BUFLEN)
    M = M
    cs = [Int(c) for c in coeffs]
    m = globalctr = t = 0
    σ₀ = (1-sum(cs)) * (M >> 1)
    σ = σ₀
    π₁ = 0
    info = "sign(Σ(cᵢxᵢ))⋅min(L,abs(Σ(cᵢxᵢ)))"
    nin = length(coeffs)
    nout = 1
    memoryless = false
    SG_SDP(M, cs, m, σ₀, σ, π₁, globalctr, t,
                info, nin, nout, memoryless)
end 
Base.broadcastable(g::SG_SDP) = Ref(g)
################################################################################
function proc(g::SG_SDP, bs::Vararg{SCBit,N}) where N
    N == g.nin || throw(DimensionMismatch("Expected $(g.nin) inputs, got $N."))
    none(ismissing, bs...) || return missing
    g.m += 1
    bout = g.m <= g.π₁
    # bout = bernoulli_samp(g.π₁ / g.M)
    for (i, bi) in enumerate(bs)
        bi && (g.σ += g.cs[i])
    end
    if g.m == g.M
       π̅₁ = g.t + g.σ                   # number of "1"s required to represent output
        g.π₁  = clamp(π̅₁, 0, g.M)       # limit to an implementable number
        g.t =π̅₁ - g.π₁                  # store the excess "1"s for later
        g.m = 0
        g.σ = g.σ₀
    end
    return (g.globalctr += 1) <= g.M ? missing : bout
end

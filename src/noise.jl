struct SCNoise
    ϵ::Float64   # bit flip probability
    function SCNoise(ϵ)
        0 <= ϵ <= 1 || throw(DomainError(ϵ))
        new(ϵ)
    end
end
Base.broadcastable(n::SCNoise) = Ref(n)

################################################################################
noisybit(noise::SCNoise, b) = rand() < noise.ϵ ? bitflip(b) : b

################################################################################
function proc(g::StochasticGate, noise::SCNoise, bs...)
    bout = proc(g, bs...)
    noise.ϵ == 0 && return bout
    noisybit(noise, bout)
end

################################################################################
# direct version
function apply_noise!(b::BitVector, noise::SCNoise)
    noise.ϵ == 0 && return b
    for (i, bit) ∈ enumerate(b)
        b[i] = noisybit(noise, bit)
    end
end
################################################################################
# return a new BitVector
function apply_noise(b::BitVector, noise)
    bout = similar(b)
    apply_noise!(bout, noise)
    bout
end

################################################################################
# The following samplers are a lot faster than Distributions.Bernoulli
bernoulli_samp(λ::AbstractFloat) = rand() < λ
bernoulli_samp(n::Integer, λ::AbstractFloat) = BitVector(rand() < λ for _ in 1:n)
function bernoulli_samp!(bv::BitVector, λ::AbstractFloat)
    for i in 1:length(bv)
        bv[i] = rand() < λ
    end
end

################################################################################
function get_test_input(N, λ, n::SCNoise)
    b = StochasticApprox.bernoulli_samp(N, float(λ))
    n.ϵ > 0 && apply_noise!(b, n)
    return b
end
get_test_input(N, λ) = get_test_input(N, λ, SCNoise(0))
################################################################################
function get_test_input(N, n::SCNoise)
    λ = rand()
    b = get_test_input(N, λ, n)
    return λ, b
end
get_test_input(N) = get_test_input(N, SCNoise(0))


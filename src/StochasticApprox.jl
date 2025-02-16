module StochasticApprox

using Statistics, Random

abstract type StochasticGate end

SCBit = Union{Missing, Bool}

include("noise.jl")
include("source.jl")
include("utils.jl")
include("bipolar.jl")
include("clenshaw.jl")

end # module

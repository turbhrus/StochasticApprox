
none(pred, bs::Vararg{SCBit,N}) where {N} = !any(pred, bs)

################################################################################
bitflip(b::Bool) = b ⊻ true

################################################################################
evalbits(g::StochasticGate, M::Function, bs...) = StochasticApprox.proc.(g, bs...) |> skipmissing |> mean |> M
evalbits(g::StochasticGate, bs...) = evalbits(g, identity, bs...)


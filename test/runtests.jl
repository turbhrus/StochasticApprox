using Test
using StochasticApprox
using Statistics


################################################################################
@testset "bipolar" begin
    N = 100_000
    #################################################
    L = 10.0
    M = λ -> L*(2λ-1)    # bipolar mapping
    Minv = x -> (x+L)/2L   # inverse map

    #################################################
    # product
    g = StochasticApprox.SGMulBipolar()
    λ1, b1 = StochasticApprox.get_test_input(N)
    λ2, b2 = StochasticApprox.get_test_input(N)
    x1 = M(λ1)
    x2 = M(λ2)
    @test isapprox(x1*x2/L, StochasticApprox.evalbits(g, M, b1, b2), atol=0.02L)
    #################################################
    cs = [1, 2, -3]
    g = StochasticApprox.SG_SDP(cs)
    b1 = StochasticApprox.get_test_input(N, Minv(0.9))
    b2 = StochasticApprox.get_test_input(N, Minv(0.4))
    b3 = StochasticApprox.get_test_input(N, Minv(0.6))
    @test isapprox(-0.1, StochasticApprox.evalbits(g, M, b1, b2, b3), atol=0.02*3L)
    ##################
    cs = rand(-5:5, 6)
    g = StochasticApprox.SG_SDP(cs)
    λ1, b1 = StochasticApprox.get_test_input(N); x1 = M(λ1)
    λ2, b2 = StochasticApprox.get_test_input(N); x2 = M(λ2)
    λ3, b3 = StochasticApprox.get_test_input(N); x3 = M(λ3)
    λ4, b4 = StochasticApprox.get_test_input(N); x4 = M(λ4)
    λ5, b5 = StochasticApprox.get_test_input(N); x5 = M(λ5)
    λ6, b6 = StochasticApprox.get_test_input(N); x6 = M(λ6)
    s = cs[1]*x1 + cs[2]*x2 + cs[3]*x3 + cs[4]*x4 + cs[5]*x5 + cs[6]*x6
    xout = StochasticApprox.evalbits(g, M, b1, b2, b3, b4, b5, b6)
    @test isapprox(sign(s)*min(L,abs(s)), xout, atol=0.02*6L)
end


################################################################################
@testset "clenshaw" begin
    N = 1_000_000
    M = λ -> L*(2λ-1)   # bipolar mapping
    Minv = x -> (x+L)/2L   # inverse map
    L = 4

    f(x) = 0.4*(cos(11x^4)) * (2x-0.5)
    a = [-0.04204922469426184,-0.12020408562146889,0.14420049219925812,-0.1983127465052685,-0.045044118946623836,0.05055150116936228,0.019768368361942682,-0.16935628776505465,0.06490977552058465,-0.037176590047742604,-0.04632148049671333,0.3049355595022372,-0.1061462992544053,0.29759692441872465,-0.042652162954956975,0.04041205958706344,0.02244613316142525,-0.1010001840598668,0.02805395886850814,-0.07292563865327949]
    g = StochasticApprox.SGClenshaw(Minv.(a))
    x = 2*rand()-1
    b = StochasticApprox.get_test_input(N, Minv(L*x))
    @test isapprox(f(x), StochasticApprox.evalbits(g, M, b), atol=0.02*L)
end
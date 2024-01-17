using ElementalBalancing
using Test
using Revise
using LinearAlgebra

A = [1 0 1 1;
    0 -0.286 -0.286 0.014;
    0 0.572 0.572 -0.028;
    0 0.858 0.858 -0.042]

E = [1 1 1 0;
    4.113 4 0 -4]
i_known = [2,3,4]
i_unknown = [1]
Em = E[:,i_known]
Ec = E[:,i_unknown]

reduced(redundancy(Ec, Em)) |> display

prob = EBALProblem(Em, Ec, Matrix{Float64}(LinearAlgebra.I, length(i_known), length(i_known)))
reduced(redundancy(Ec, Em))

solve(prob, [0.1, 0.1, 0.1]) |> display


@testset "ElementalBalancing.jl" begin
    # Write your tests here.
end




redundancy(Ec, Em) |> reduced

Ec_star=(inv(Ec'*Ec))*Ec'
R=Em-Ec*Ec_star*Em
U,S,V=svd(R)
Sconv=[1 0]
C=Sconv*S
K=C*S'*U'
Rred=K*R
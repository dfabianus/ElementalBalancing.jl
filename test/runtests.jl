using ElementalBalancing
using Test
using Revise

@testset "ElementalBalancing.jl" begin
    include("integration_test_1.jl")
    include("soft_sensing_test_1.jl")
end
# This type contains the data for the ElementalBalancing
mutable struct EBALProblem
    E::Matrix{Float64}
    M::Vector{Float64}
    i_known::Vector{Int}
end

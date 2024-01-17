mutable struct EBALProblem
    E::Matrix{Float64}
    M::Vector{Float64}
    i_known::Vector{Int}
    i_unknown::Vector{Int}
end

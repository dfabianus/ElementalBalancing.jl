using ElementalBalancing
using Test
using LinearAlgebra

A = [1 0 1 1;
    0 -0.286 -0.286 0.014;
    0 0.572 0.572 -0.028;
    0 0.858 0.858 -0.042]
E = [1 1 1 0;
    4.113 4 0 -4]
i_known = [2,3,4] # i_unknown = [1], automatically inferred
M = [26.5, 30, 44, 32]
F = diagm([0.1, 0.1, 0.1])

prob = EBALProblem(E, M, i_known)
reduced(redundancy(prob))
sol = solve(prob, [0.1, 0.1, 0.1], F)
sol.h

@test sol.h > 0.2
@test sol.h < 0.3
@test length(sol.rm_best_mol) == length(i_known)
@test length(sol.rc_best_mol) == 1
@test sol.rm_best_g == sol.rm_best_mol .* M[i_known]
@test sol.rc_best_g == sol.rc_best_mol .* M[1]
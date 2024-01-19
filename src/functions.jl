i_unknown(prob::EBALProblem) = setdiff(collect(range(1,size(prob.E,2))), prob.i_known)

# Splitting E and M into known and unknown parts
Em(prob::EBALProblem) = prob.E[:,prob.i_known]
Ec(prob::EBALProblem) = prob.E[:,i_unknown(prob)]
Mm(prob::EBALProblem) = prob.M[prob.i_known]
Mc(prob::EBALProblem) = prob.M[i_unknown(prob)]
## todo: function to recombine from known and unknown parts

# retranslation between mole and mass basis
translate_rg(r, M) = r .* M
translate_rm(rg, M) = rg ./ M

# Equations from the book chapter: Gregory N., Material balances and data consistency [1]:
# https://edisciplinas.usp.br/pluginfile.php/5381165/mod_resource/content/1/Chap4.pdf
solve_rc(Ec, Em, rm) = -pinv(Ec) * Em * rm              # eq. 4.14
redundancy(Ec, Em) = Em-Ec*pinv(Ec)*Em                  # eq. 4.17
eps(Rred, rm) = Rred * rm                               # eq. 4.20
P(Rred, F) = Rred * F *Rred'                            # eq. 4.24
delta(Rred, rm, F, P) = (F*Rred'*inv(P) * Rred)* rm     # eq. 4.26
rm_best(rm, delta) = rm-delta                           # eq. 4.27
h(eps, P) = eps' * inv(P) * eps                         # eq. 4.29

# Transfering the equations to the EBALProblem type
solve_rc(prob::EBALProblem, rm)     = solve_rc(Ec(prob), Em(prob), rm)
redundancy(prob::EBALProblem)       = redundancy(Ec(prob), Em(prob))
eps(prob::EBALProblem, rm)          = eps(reduced(redundancy(prob)), rm)
P(prob::EBALProblem, F)             = P(reduced(redundancy(prob)), F)
delta(prob::EBALProblem, rm, F)     = delta(reduced(redundancy(prob)), rm, F, P(prob, F))
rm_best(prob::EBALProblem, rm, F)   = rm_best(rm, delta(prob, rm, F))
h(prob::EBALProblem, rm, F)         = h(eps(prob, rm), P(prob, F))

# Reducing a matrix to its independent rows, 
# important for the reduced redundancy matrix
function reduced(A)
    Q, R = qr(A)
    tolerance = 1e-10
    independent_rows = findall(x -> abs(x) > tolerance, diag(R))
    return A[independent_rows, :]
end

# The main function of the package
function solve(prob::EBALProblem, rm, F)
    out_rm_best_mol = rm_best(prob, rm, F)
    out_rc_best_mol = solve_rc(prob, out_rm_best_mol)
    out_rm_best_g = translate_rg(out_rm_best_mol, Mm(prob))
    out_rc_best_g = translate_rg(out_rc_best_mol, Mc(prob))
    out_r_best_mol = zeros(size(prob.E,2))
    out_r_best_mol[prob.i_known] = out_rm_best_mol
    out_r_best_mol[i_unknown(prob)] = out_rc_best_mol
    out_r_best_g = translate_rg(out_r_best_mol, prob.M)
    return (
        rm_best_mol = out_rm_best_mol,
        rc_best_mol = out_rc_best_mol,
        r_best_mol = out_r_best_mol,
        rm_best_g = out_rm_best_g,
        rc_best_g = out_rc_best_g,
        r_best_g = out_r_best_g,
        h = h(prob, rm, F),
    )
end
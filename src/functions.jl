

function solve(prob::EBALProblem)
    
end

function solve_rc(Ec, Em, rm)
    return -pinv(Ec) * Em * rm
end

function redundancy(Ec, Em)
    return Em-Ec*pinv(Ec)*Em
end

function reduced(A)
    Q, R = qr(A)
    tolerance = 1e-10
    independent_rows = findall(x -> abs(x) > tolerance, diag(R))
    return A[independent_rows, :]
end

# eps: residual vector
function eps(Rred, rm)
    return Rred * rm
end

# P: Residual variance covariance matrix
function P(Rred, sigma)
    return Rred * sigma *Rred'
end

# Reconciliation of measured and calculated rates
function delta(Rred, rm, sigma)
    return (sigma*Rred'*inv(P) * Rred)* rm
end

function rm_best(rm, delta)
    return rm-delta
end

# Sum of weighted squares of residuals
function h(eps, P)
    return eps' * inv(P) * eps
end

solve_rc(prob::EBALProblem, rm) = solve_rc(prob.Ec, prob.Em, rm)
redundancy(prob::EBALProblem) = redundancy(prob.Ec, prob.Em)
eps(prob::EBALProblem, rm) = eps(reduced(redundancy(prob)), rm)
P(prob::EBALProblem) = P(reduced(redundancy(prob)), prob.sigma)
delta(prob::EBALProblem, rm) = delta(reduced(redundancy(prob)), rm, prob.sigma)
rm_best(prob::EBALProblem, rm) = rm_best(rm, delta(prob, rm))
h(prob::EBALProblem, rm) = h(eps(prob, rm), P(prob))

A = [1 0 1 1;
    0 -0.286 -0.286 0.014;
    0 0.572 0.572 -0.028;
    0 0.858 0.858 -0.042]

reduced(A)

E = [1 1 1 0;
    4.113 4 0 -4]
i_known = [2,3,4]
i_unknown = [1]
Em = E[:,i_known]
Ec = E[:,i_unknown]

redundancy(Ec, Em) |> reduced

Ec_star=(inv(Ec'*Ec))*Ec'
R=Em-Ec*Ec_star*Em
U,S,V=svd(R)
Sconv=[1 0]
C=Sconv*S
K=C*S'*U'
Rred=K*R
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
function delta(Rred, rm, sigma, P)
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
delta(prob::EBALProblem, rm) = delta(reduced(redundancy(prob)), rm, prob.sigma, P(prob))
rm_best(prob::EBALProblem, rm) = rm_best(rm, delta(prob, rm))
h(prob::EBALProblem, rm) = h(eps(prob, rm), P(prob))

function solve(prob::EBALProblem, rm)
    return (
        rm_best = rm_best(prob, rm),
    )
end
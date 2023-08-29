module ElementalBalancing

import LinearAlgebra
import BioprocessKinetics
# Write your package code here.
function K2S1(r, E, i_known, i_unknown)
    # Q is a vector of the supply rates and exhaust rates measured at the bioreactor
    # Units are mol/h
    sigma = Matrix{Float64}(LinearAlgebra.I, length(i_known), length(i_known))
    rm = r[i_known]
    Em = E[:,i_known]
    Ec = E[:,i_unknown]

    #Ec_star: Moore Penrose pseudo inverse of Ec
    Ec_star=(inv(Ec'*Ec))*Ec'
    # R: Redundancy matrix
    R=Em-Ec*Ec_star*Em
    # Rred: Reduced redundancy matrix (containing just the independent rows of R
    U,S,V=svd(R)
    Sconv=[1 0]
    C=Sconv*S
    K=C*S'*U'
    Rred=K*R
    # eps: residual vector
    eps = Rred * rm
    # P: Residual variance covariance matrix
    P = Rred * sigma *Rred'
    # Reconciliation of measured and calculated rates
    delta = (sigma*Rred'*inv(P) * Rred)* rm
    rm_best = rm-delta
    xc_best = -Ec_star*Em*rm_best
    # Sum of weighted squares of residuals
    h = eps' * inv(P) * eps
    # Calculate the function outputs
    return (r=vcat(xc_best, rm_best), h=h)
end

function K2S1m!(dx,x,p,t)

    i_known = [2,3,4]
    i_unknown = [1]

    M = [26.5, 30, 44, 32]

    E = [1 1 1 0;
         4.113 4 0 -4]
    
    c_S = x[2] * p.V_L(t)
    qS = p.qSmax(t) * c_S / (p.kS(t) + c_S)
    rS = -qS * x[1]

    # constructing Q and r in mol/h and g/h
    Q_mol = vcat(0, p.Q_S(t), p.Q_CO2(t), p.Q_O2(t))
    Q_g = Q_mol .* M
    r_g = vcat(0, rS, Q_g[3,4])
    r_mol = r_g ./ M

    r_hat_mol, h = K2S1(r_mol, E, i_known, i_unknown)
    r_hat_g = r_hat_mol .* M
    # if x[1] < 7
    #     #r_hat[1] = -0.37 * r_hat[2]
    #     r_hat[1] = -0.45 * r_hat[2]
    # end
    
    #dx[:] = [Qᵢ+rᵢ for (Qᵢ,rᵢ) in zip(Q_g,r_hat_g)]
    dx[:] = Q_g .+ r_hat_g
end

end

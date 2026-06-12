# Method-of-moments drop-in for RTIGER's BetaBinomial M-step.
#
# Redefines emissionUpdateState (defined in Main by sourceJulia()) so the
# dispersion tau = alpha+beta is estimated in closed form instead of by the
# numerical optimize() MLE. Same signature, same return (a_i, b_i, mi, tau_i),
# same model (BetaBinomial) and same clamps -> a drop-in faster estimator.
#
# Estimator (Williams / Pearson moment estimator for the intra-class
# correlation rho = 1/(tau+1), with varying cluster sizes n_p and gamma
# weights w_p):
#   mi = sumk/sumn                                   (weighted mean, as before)
#   For BetaBinomial:  E[(k - n*mi)^2] = n*mi*(1-mi)*(1 + (n-1)*rho)
#   So with  S = Σ_p w_p (k_p - n_p mi)^2 / (n_p mi(1-mi)),
#            W = Σ_p w_p,   Dd = Σ_p w_p (n_p - 1):
#       rho_hat = (S - W) / Dd        (closed form, one pass over distinct pairs)
#       tau = 1/rho - 1
# rho<=0 (underdispersion, common at low coverage) clamps to tau=100, matching
# the MLE path's tau cap. Dd<=0 (all n==1) -> no info on rho -> keep prior tau.

function emissionUpdateState(i, ks, ns, ws, sumk, sumn, alpha_old, beta_old)
    mi = sumk / sumn
    if isnan(mi)
        mi = alpha_old[i] / (alpha_old[i] + beta_old[i])
        mi = clamp(mi, 1e-12, 1 - 1e-12)        # guard 0/0 if a prior shape is 0
        a_i = alpha_old[i]
        b_i = beta_old[i]
        tau_i = (alpha_old[i] / mi + beta_old[i] / (1 - mi))
        return a_i, b_i, mi, tau_i
    end
    mi < 0.01 && (mi = 0.01)
    mi > 0.99 && (mi = 0.99)
    v = mi * (1 - mi)

    S = 0.0
    W = 0.0
    Dd = 0.0
    @inbounds for p = 1:length(ws)
        np = ns[p]
        np == 0 && continue
        resid = ks[p] - np * mi
        S += ws[p] * resid * resid / (np * v)
        W += ws[p]
        Dd += ws[p] * (np - 1)
    end

    if Dd <= 0.0
        tau_i = (alpha_old[i] / mi + beta_old[i] / (1 - mi))   # no replication -> prior
    else
        rho = (S - W) / Dd
        rho = clamp(rho, 1.0 / 101.0, 0.999)                   # tau in (~0, 100]
        tau_i = 1.0 / rho - 1.0
    end
    tau_i > 100 && (tau_i = 100.0)

    a_i = tau_i * mi
    b_i = tau_i * (1 - mi)
    return a_i, b_i, mi, tau_i
end

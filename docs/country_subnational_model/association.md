# Feasible TBI and medical factor association

## Within age joint distribution

For each disjoint age band b, let p = P(TBI), q = P(any included non-recency factor), and J = P(TBI and any factor). Both marginals lie in [0,1]. The Fréchet bounds are L = max(0,p+q−1) and U = min(p,q). Independence has J0 = p*q.

For the requested normalized association alpha in [−1,1]:

    alpha < 0: J = J0 + alpha * (J0 − L)
    alpha >= 0: J = J0 + alpha * (U − J0)

Thus −1 gives L, zero gives J0, and +1 gives U. This is a piecewise linear overlap scale, not Pearson correlation. It remains defined when one marginal is zero or one, though then the overlap is fixed and alpha has no effect.

| | Any factor | No factor |
| --- | ---: | ---: |
| TBI | J | p−J |
| No TBI | q−J | 1−p−q+J |

For nondegenerate marginals, achieved phi is:

    phi = (J − p*q) / sqrt(p*(1−p)*q*(1−q))

Its attainable limits are the same expression evaluated at L and U, generally not −1 and +1. For degenerate marginals, return null with reason undefined_degenerate_marginal, never a fabricated zero. The joint proportion with both conditions is J; the proportion among TBI with any factor is J/p when p > 0.

Example p = 0.30, q = 0.20: L = 0, J0 = 0.06, U = 0.20. The requested scale endpoints achieve phi approximately −0.3273 and +0.7638. No attainable joint has phi = +1 here.

Reject nonfinite values, booleans, out-of-range marginals/alpha, negative joint cells, overlaps outside [L,U] and distributions not summing to one. Do not repair invalid inputs by clipping. The pure prototype uses exact rational arithmetic from JSON-number decimal representations, with decimal square roots for phi, before conversion to JSON numbers. This avoids tiny negative endpoint cells from binary cancellation.

## Age driven and overall association

Let w_b be demographic weights summing to one. Overall P = sum(w_b*p_b), Q = sum(w_b*q_b), and J_all = sum(w_b*J_b). Compute overall cells and phi from P,Q,J_all; never average phi or alpha across bands. If all within-age alpha values are zero, overall phi can still be nonzero because age induces association.

The covariance decomposition makes the distinction explicit:

    J_all − P*Q = sum(w_b*(J_b − p_b*q_b))
                 + [sum(w_b*p_b*q_b) − P*Q]

The first term is residual within-age covariance; the second is age-driven covariance. Report both. With sex strata, use their demographic weights and disclose whether the control applies within age only or within age and sex; introducing sex-specific control changes the estimand.

Every report includes age band, weight, requested alpha, p, q, L, U, achieved J, four cells, phi or null reason, and TBI-with-any-factor proportion. The overall report has the aggregate quantities and the vector of requested age controls; there is no unique overall requested alpha.

## Individual factors and any factor prevalence

The aggregate joint does not determine individual-factor overlap, within-TBI composition or multimorbidity. Recency is excluded from q. Included factor marginals p_j,b are by age. Missing values cannot be assumed zero; exclusion and non-applicability require reasons.

For a parsimonious first implementation, generate the unconditioned factor pattern from independent Bernoulli p_j,b within age. This defines q_b = 1 − product_j(1−p_j,b). Then couple TBI to its any-factor indicator R using:

    P(TBI | R=1,b) = J_b / q_b                     if q_b > 0
    P(TBI | R=0,b) = (p_b−J_b) / (1−q_b)          if q_b < 1

Equivalently generate the aggregate TBI/R cells first, then draw a nonempty factor pattern with probability product_j p_j^f_j (1−p_j)^(1−f_j) / q_b in either R-positive TBI stratum. R-negative people receive no factors. This truncation induces conditional dependence; it must not be called independence given age and TBI. It preserves every unconditional age-specific factor prevalence and the selected q and J in expectation. Do not sample independent factors within R-positive people without truncation, since some would have no factors and break the requested union.

If q is separately observed, require max_j p_j <= q <= min(1,sum_j p_j) and additional joint feasibility. The independent-union formula cannot silently replace observed q. Reconcile using a reviewed pattern model or report incompatibility. Exact finite populations may not attain the requested continuous probabilities: store expected and realized cells, seed, weights and sampling deviation. Weighted deterministic populations can retain exact expected proportions; never silently round percentages and claim exact agreement.

## Alternative allocation methods

| Method | Inputs and strengths | Limitations and feasibility checks |
| --- | --- | --- |
| Conditional independence given age and TBI | p_j conditional on TBI and no TBI, or approved factor-specific associations | Marginal factors and aggregate alpha alone do not identify these values. Unions in each TBI stratum must agree with the aggregate joint; naive independent draws need not do so. |
| User pairwise associations | Reviewed factor-pair data with age, definition and conditioning | Pairwise Fréchet feasibility and positive-semidefinite correlation are necessary but not sufficient for a multivariate Bernoulli joint. Many full joints fit the same pairs. |
| Gaussian copula or latent variable | Marginal thresholds and a reviewed latent correlation structure; scalable sampling | Latent correlations differ from observed binary phi. Calibration and positive-semidefinite constraints can fail; tail dependence and aggregate union may be wrong. |
| Iterative proportional fitting | Seed contingency table plus compatible marginal/selected joint constraints | Sparse high-dimensional tables and structural zeros cause instability or nonexistence. Pairwise margins leave higher-order interactions unidentified; convergence does not establish scientific validity. |

Start with the explicit independent age-pattern plus aggregate coupling policy. Validate against available multimorbidity data and compare selected plausible dependent scenarios before adding many unidentifiable parameters. Factor association with recency within total TBI is an additional unobserved assumption; initially use the same conditional infection-time distribution within each age/TBI stratum and disclose it.

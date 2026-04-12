# ─────────────────────────────────────────────────────────────────────────────
# Mean-field (Landau) universality class — exact critical exponents
#
# Mean-field exponents are the baseline reference for all universality
# classes. They are exact above the upper critical dimension d_c (= 4
# for standard Ising-like transitions) and serve as the starting point
# for ε-expansion (d = d_c − ε) corrections.
#
# References:
#   L. D. Landau, "On the theory of phase transitions",
#     Zh. Eksp. Teor. Fiz. 7, 19 (1937).
#   J. Zinn-Justin, "Quantum Field Theory and Critical Phenomena",
#     Oxford University Press (2002), Ch. 23.
# ─────────────────────────────────────────────────────────────────────────────

"""
    MeanField

Dispatch tag for the mean-field (Landau) universality class.

These exponents are exact for spatial dimension d ≥ 4 (the upper
critical dimension for Ising-like transitions). Below d = 4, critical
fluctuations renormalize the exponents (see `Ising2D` for d = 2).
"""
struct MeanField end

"""
    fetch(::MeanField, ::CriticalExponents) -> NamedTuple

Exact critical exponents of the mean-field universality class:

| Symbol | Value | Physical meaning                              |
|--------|-------|-----------------------------------------------|
| β      | 1/2   | Order parameter: M ∼ (T_c − T)^β             |
| ν      | 1/2   | Correlation length: ξ ∼ |T − T_c|^{−ν}       |
| γ      | 1     | Susceptibility: χ ∼ |T − T_c|^{−γ}           |
| η      | 0     | Anomalous dimension: G(r) ∼ r^{−(d−2+η)}     |
| δ      | 3     | Critical isotherm: M ∼ h^{1/δ} at T = T_c    |
| α      | 0     | Specific heat: discontinuity (step function)  |

All values are `Rational{Int}` for exact representation.
Satisfies the standard scaling relations:
  α + 2β + γ = 2,   γ = β(δ−1),   γ = ν(2−η).

# References
    L. D. Landau, Zh. Eksp. Teor. Fiz. 7, 19 (1937).
"""
function fetch(::MeanField, ::CriticalExponents; kwargs...)
    return (
        β=1 // 2,    # order parameter exponent
        ν=1 // 2,    # correlation length exponent
        γ=1 // 1,    # susceptibility exponent
        η=0 // 1,    # anomalous dimension
        δ=3 // 1,    # critical isotherm exponent
        α=0 // 1,    # specific heat (discontinuity)
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# 2D Ising universality class — exact critical exponents
#
# The 2D Ising model is the archetype of a second-order phase transition
# in two dimensions. Its universality class is characterized by the
# Virasoro minimal model M(3,4) with central charge c = 1/2.
#
# All critical exponents are known exactly (rational numbers).
#
# References:
#   L. Onsager, Phys. Rev. 65, 117 (1944).
#   C. N. Yang, Phys. Rev. 85, 808 (1952).
#   A. A. Belavin, A. M. Polyakov, A. B. Zamolodchikov,
#     Nucl. Phys. B 241, 333 (1984) — conformal field theory.
# ─────────────────────────────────────────────────────────────────────────────

"""
    Ising2D

Dispatch tag for the 2D Ising universality class.

Central charge c = 1/2 (Virasoro minimal model M(3,4)).
Applies to: 2D classical Ising model, 1+1D TFIM at criticality,
and any system in the same universality class.
"""
struct Ising2D end

"""
    CriticalExponents

Dispatch tag for the set of critical exponents of a universality class.
Returns a `NamedTuple` of exact values (typically `Rational{Int}`).
"""
struct CriticalExponents end

"""
    fetch(::Ising2D, ::CriticalExponents) -> NamedTuple

Exact critical exponents of the 2D Ising universality class:

| Symbol | Value | Physical meaning                              |
|--------|-------|-----------------------------------------------|
| β      | 1/8   | Order parameter: M ∼ (T_c − T)^β             |
| ν      | 1     | Correlation length: ξ ∼ |T − T_c|^{−ν}       |
| γ      | 7/4   | Susceptibility: χ ∼ |T − T_c|^{−γ}           |
| η      | 1/4   | Anomalous dimension: G(r) ∼ r^{−(d−2+η)}     |
| δ      | 15    | Critical isotherm: M ∼ h^{1/δ} at T = T_c    |
| α      | 0     | Specific heat: C ∼ ln|T − T_c| (logarithmic) |
| c      | 1/2   | Central charge (CFT)                          |

All values are `Rational{Int}` for exact representation.

# References
    A. A. Belavin, A. M. Polyakov, A. B. Zamolodchikov,
      Nucl. Phys. B 241, 333 (1984).
"""
function fetch(::Ising2D, ::CriticalExponents; kwargs...)
    return (
        β=1 // 8,    # order parameter exponent
        ν=1 // 1,    # correlation length exponent
        γ=7 // 4,    # susceptibility exponent
        η=1 // 4,    # anomalous dimension
        δ=15 // 1,   # critical isotherm exponent
        α=0 // 1,    # specific heat (logarithmic divergence)
        c=1 // 2,    # central charge
    )
end

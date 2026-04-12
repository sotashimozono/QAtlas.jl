# ─────────────────────────────────────────────────────────────────────────────
# Potts universality classes (q-state, d=2 exact)
#
# The q-state Potts model generalizes the Ising model (q=2).
# For d=2, exact exponents are known for q ≤ 4 via Coulomb gas / CFT.
#
# References:
#   R. J. Baxter, "Exactly Solved Models in Statistical Mechanics" (1982).
#   B. Nienhuis, J. Stat. Phys. 34, 731 (1984).
# ─────────────────────────────────────────────────────────────────────────────

"""
    fetch(::Universality{:Potts3}, ::CriticalExponents; d=2) -> NamedTuple

Exact critical exponents of the 3-state Potts model in d=2 (S₃ symmetry).
"""
function fetch(::Universality{:Potts3}, ::CriticalExponents; d::Int=2, kwargs...)
    if d == 2
        return (
            α=1 // 3,
            β=1 // 9,
            γ=13 // 9,
            δ=14 // 1,
            ν=5 // 6,
            η=4 // 15,
        )
    end
    return error("Potts3 universality: only d=2 implemented; got d=$d.")
end

"""
    fetch(::Universality{:Potts4}, ::CriticalExponents; d=2) -> NamedTuple

Exact critical exponents of the 4-state Potts model (Ashkin–Teller point)
in d=2 (S₄ symmetry). This is the marginal case: the transition is
second-order with logarithmic corrections.
"""
function fetch(::Universality{:Potts4}, ::CriticalExponents; d::Int=2, kwargs...)
    if d == 2
        return (
            α=2 // 3,
            β=1 // 12,
            γ=7 // 6,
            δ=15 // 1,
            ν=2 // 3,
            η=1 // 4,
        )
    end
    return error("Potts4 universality: only d=2 implemented; got d=$d.")
end

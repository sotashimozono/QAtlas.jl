# ─────────────────────────────────────────────────────────────────────────────
# Percolation universality class
#
# References:
#   d=2: Nienhuis, Phys. Rev. Lett. 49, 1062 (1982) — exact via Coulomb gas.
#   d=3: Wang, Zhou, Zhang, Garoni, Deng, Phys. Rev. E 87, 052107 (2013).
# ─────────────────────────────────────────────────────────────────────────────

"""
    fetch(::Universality{:Percolation}, ::CriticalExponents; d) -> NamedTuple

Critical exponents of the ordinary (site/bond) percolation universality class.

- **d = 2**: Exact rational values (Coulomb gas mapping).
- **d = 3**: Numerical (Monte Carlo).
- **d ≥ 6**: Mean-field.
"""
function fetch(::Universality{:Percolation}, ::CriticalExponents; d::Int, kwargs...)
    if d == 2
        return (
            α=-2 // 3,
            β=5 // 36,
            γ=43 // 18,
            δ=91 // 5,
            ν=4 // 3,
            η=5 // 24,
        )
    elseif d == 3
        return (
            α=-0.625,
            α_err=0.003,
            β=0.4181,
            β_err=0.0008,
            γ=1.793,
            γ_err=0.003,
            δ=5.29,
            δ_err=0.06,
            ν=0.87619,
            ν_err=0.00012,
            η=0.46,
            η_err=0.08,
        )
    elseif d >= 6
        # Percolation upper critical dimension is 6
        return (α=-1 // 1, β=1 // 1, γ=1 // 1, δ=2 // 1, ν=1 // 2, η=0 // 1)
    end
    return error("Percolation universality: d=$d not supported (d ∈ {2, 3, ≥6}).")
end

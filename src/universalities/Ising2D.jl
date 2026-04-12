# ─────────────────────────────────────────────────────────────────────────────
# Ising universality class — exact (d=2) and numerical (d=3) exponents
#
# References:
#   d=2: Belavin, Polyakov, Zamolodchikov, Nucl. Phys. B 241, 333 (1984).
#   d=3: Kos, Poland, Simmons-Duffin, Vichi, JHEP 08, 036 (2016)
#         — conformal bootstrap bounds.
# ─────────────────────────────────────────────────────────────────────────────

# Backward-compatible alias — delegates to Universality{:Ising} with d=2
struct Ising2D end

"""
    fetch(::Ising2D, ::CriticalExponents) -> NamedTuple

Backward-compatible alias for `fetch(Universality(:Ising), CriticalExponents(); d=2)`.
"""
function fetch(::Ising2D, ::CriticalExponents; kwargs...)
    return fetch(Universality(:Ising), CriticalExponents(); d=2, kwargs...)
end

"""
    fetch(::Universality{:Ising}, ::CriticalExponents; d) -> NamedTuple

Critical exponents of the Ising universality class (Z₂ symmetry).

- **d = 2**: Exact rational values (CFT minimal model M(3,4), c = 1/2).
- **d = 3**: High-precision numerical estimates from the conformal
  bootstrap (Kos et al. 2016). Fields `α_err`, `β_err`, … give the
  uncertainty in the last digits.
- **d ≥ 4**: Mean-field (Landau) exponents (upper critical dimension).
"""
function fetch(::Universality{:Ising}, ::CriticalExponents; d::Int, kwargs...)
    if d == 2
        return (α=0 // 1, β=1 // 8, γ=7 // 4, δ=15 // 1, ν=1 // 1, η=1 // 4, c=1 // 2)
    elseif d == 3
        # Conformal bootstrap: Kos, Poland, Simmons-Duffin, Vichi (2016)
        return (
            α=0.11009,
            α_err=0.00001,
            β=0.32642,
            β_err=0.00001,
            γ=1.23708,
            γ_err=0.00001,
            δ=4.78984,
            δ_err=0.00001,
            ν=0.62997,
            ν_err=0.00001,
            η=0.03630,
            η_err=0.00005,
        )
    elseif d >= 4
        return fetch(MeanField(), CriticalExponents())
    end
    return error("Ising universality: d=$d not supported (d ∈ {2, 3, ≥4}).")
end

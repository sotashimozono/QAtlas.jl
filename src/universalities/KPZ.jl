# ─────────────────────────────────────────────────────────────────────────────
# KPZ (Kardar–Parisi–Zhang) universality class
#
# The KPZ equation is a non-equilibrium growth model. Its scaling
# exponents differ from equilibrium critical exponents and are accessed
# via the `GrowthExponents` tag rather than `CriticalExponents`.
#
# References:
#   M. Kardar, G. Parisi, Y.-C. Zhang, Phys. Rev. Lett. 56, 889 (1986).
#   T. Sasamoto, H. Spohn, Nucl. Phys. B 834, 523 (2010).
# ─────────────────────────────────────────────────────────────────────────────

"""
    KPZ1D <: AbstractQAtlasModel

Dispatch tag for the **1+1-dimensional** KPZ (Kardar–Parisi–Zhang)
universality class. Acts as a convenience wrapper over
`Universality(:KPZ)` with the dimension pinned to `d = 1`:

```julia
QAtlas.fetch(KPZ1D(), CriticalExponents())
# ≡ QAtlas.fetch(Universality(:KPZ), GrowthExponents(); d = 1)
```

The struct is kept as a distinct nominal type (rather than
`const KPZ1D = Universality{:KPZ}`) because the dimension is fixed
here — a type-level alias would defeat that fix. Subtyped to
`AbstractQAtlasModel` so dispatch composes with the canonical
`fetch(::AbstractQAtlasModel, ::AbstractQuantity, ::BoundaryCondition)`
signature.

See also: [`Universality`](@ref), [`GrowthExponents`](@ref).
"""
struct KPZ1D <: AbstractQAtlasModel end
function fetch(::KPZ1D, ::CriticalExponents; kwargs...)
    return fetch(Universality(:KPZ), GrowthExponents(); d=1, kwargs...)
end

"""
    fetch(::Universality{:KPZ}, ::GrowthExponents; d) -> NamedTuple

Scaling exponents of the KPZ universality class.

- **d = 1** (1+1D): All three exponents are exact (Rational{Int}).
  Scaling relations: α + z = 2 (Galilean invariance), β = α / z.

Higher dimensions (d ≥ 2) are not exactly known and are not yet
included.
"""
function fetch(::Universality{:KPZ}, ::GrowthExponents; d::Int, kwargs...)
    if d == 1
        return (β_growth=1 // 3, α_rough=1 // 2, z=3 // 2)
    end
    return error("KPZ growth exponents: only d=1 (1+1D) implemented; got d=$d.")
end

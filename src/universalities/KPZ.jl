# ─────────────────────────────────────────────────────────────────────────────
# KPZ (Kardar–Parisi–Zhang) universality class — exact scaling exponents
#
# The KPZ equation describes the nonequilibrium growth of interfaces:
#
#   ∂h/∂t = ν ∇²h + (λ/2)(∇h)² + η(x,t)
#
# In 1+1 dimensions all three scaling exponents are known exactly
# (rigorous, not conjectured) from the connection to directed polymers
# and the replica method.
#
# References:
#   M. Kardar, G. Parisi, Y.-C. Zhang, "Dynamic Scaling of Growing
#     Interfaces", Phys. Rev. Lett. 56, 889 (1986).
#   M. Prähofer, H. Spohn, "Exact Scaling Functions for One-Dimensional
#     Stationary KPZ Growth", J. Stat. Phys. 115, 255 (2004) — exact
#     scaling function via Tracy–Widom distribution.
#   T. Sasamoto, H. Spohn, "Exact height distributions for the KPZ
#     equation with narrow wedge initial condition",
#     Nucl. Phys. B 834, 523 (2010).
# ─────────────────────────────────────────────────────────────────────────────

"""
    KPZ1D

Dispatch tag for the KPZ universality class in 1+1 dimensions.

The 1+1D KPZ equation has exact scaling exponents connected via
the Galilean invariance relation α + z = 2 and the scaling relation
β = α / z.
"""
struct KPZ1D end

"""
    fetch(::KPZ1D, ::CriticalExponents) -> NamedTuple

Exact scaling exponents of the 1+1D KPZ universality class:

| Symbol    | Value | Physical meaning                           |
|-----------|-------|--------------------------------------------|
| β_growth  | 1/3   | Growth: width ∼ t^{β}                      |
| α_rough   | 1/2   | Roughness: width ∼ L^{α} (saturated)       |
| z         | 3/2   | Dynamic: t_sat ∼ L^{z}                     |

Scaling relations: α + z = 2 (Galilean invariance), β = α / z.

All values are `Rational{Int}` for exact representation.

# References
    M. Kardar, G. Parisi, Y.-C. Zhang, Phys. Rev. Lett. 56, 889 (1986).
"""
function fetch(::KPZ1D, ::CriticalExponents; kwargs...)
    return (
        β_growth=1 // 3,   # growth exponent (width ~ t^β)
        α_rough=1 // 2,    # roughness exponent (width ~ L^α)
        z=3 // 2,          # dynamic exponent (t_sat ~ L^z)
    )
end

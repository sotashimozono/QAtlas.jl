# test/util/thermodynamic_identities.jl — model-agnostic self-validation harness.
#
# Issue #117: every QAtlas implementation that exposes (Energy, FreeEnergy,
# ThermalEntropy, SpecificHeat, ...) at a given BC should obey universal
# thermodynamic identities like
#
#   ε(β) = f(β) + T·s(β)               (Gibbs)
#   c_v(β) = -β² ∂ε/∂β                 (specific heat from energy)
#   m_α(h) = -∂f/∂h_α                  (magnetisation from free energy)
#   χ_αα   = ∂m_α/∂h_α = β·Var(M_α)/N  (linear response from variance)
#
# These hold *between the model's own outputs* — we are not comparing to
# literature here; we are checking internal consistency of the dispatch
# layer.  Self-validation catches per-site/total drift, sign errors,
# missing field-dependence, and ForwardDiff-incompatible kwargs.
#
# This file deliberately lives in `test/util/`: the harness is a
# debugging/verification tool, not part of QAtlas's public surface.  If
# downstream packages need it later, lift it to `src/verification/` as a
# weakdep extension on ForwardDiff.

using ForwardDiff
using QAtlas:
    fetch,
    AbstractQAtlasModel,
    AbstractQuantity,
    BoundaryCondition,
    Energy,
    FreeEnergy,
    ThermalEntropy,
    SpecificHeat,
    OBC,
    PBC,
    Infinite

"""
    ThermoIdentity(name, requires, check)

A single self-validation rule.

- `name::String`: human-readable label, surfaced in `IdentityCheckResult`.
- `requires::Vector{Type}`: quantity types the identity needs to be
  evaluable on the target `(model, bc)`.  The harness checks dispatch
  existence (via `which(fetch, ...)`) and skips the identity if any
  required type lacks a non-catch-all method.
- `check::Function`: `(model, bc, params::NamedTuple) -> (lhs, rhs)`.
  The two sides should be equal up to the harness's `(rtol, atol)`.
"""
struct ThermoIdentity
    name::String
    requires::Vector{Type}
    check::Function
end

"""
    IdentityCheckResult

One row per `(identity, params)` evaluated by
[`verify_thermodynamic_identities`](@ref).

`status` is `:pass`, `:fail`, or `:skipped` (the latter when one of
`identity.requires` is not dispatchable for `(model, bc)`).  For
`:skipped` rows, the numerical fields are `NaN`.
"""
struct IdentityCheckResult
    model::Any
    bc::Any
    identity::String
    params::NamedTuple
    lhs::Float64
    rhs::Float64
    abs_err::Float64
    rel_err::Float64
    status::Symbol
end

# ──────────────────────────────────────────────────────────────────────
# Default identity set
# ──────────────────────────────────────────────────────────────────────

"""
    GIBBS_RELATION

Self-validation of `ε = f + T·s` (no AutoDiff — three independent
fetches whose values must reconcile).  Catches per-site/total mix-ups,
sign errors in entropy, and missing temperature factors.
"""
const GIBBS_RELATION = ThermoIdentity(
    "Gibbs ε = f + T·s",
    Type[Energy{:per_site}, FreeEnergy, ThermalEntropy],
    function (model, bc, params)
        β = params.β
        T = 1 / β
        ε = fetch(model, Energy(:per_site), bc; beta=β)
        f = fetch(model, FreeEnergy(), bc; beta=β)
        s = fetch(model, ThermalEntropy(), bc; beta=β)
        return Float64(ε), Float64(f + T * s)
    end,
)

"""
    SPECIFIC_HEAT_FROM_ENERGY

Self-validation of `c_v = -β² ∂ε/∂β` via `ForwardDiff.derivative` on
`Energy(:per_site)`.  Requires the model's `fetch` methods to accept a
`Real` `beta` kwarg (TFIM was relaxed in PR #115; new models inherit
this requirement).
"""
const SPECIFIC_HEAT_FROM_ENERGY = ThermoIdentity(
    "c_v = -β² ∂ε/∂β  (ForwardDiff)",
    Type[Energy{:per_site}, SpecificHeat],
    function (model, bc, params)
        β = params.β
        dε_dβ = ForwardDiff.derivative(b -> fetch(model, Energy(:per_site), bc; beta=b), β)
        c_v = fetch(model, SpecificHeat(), bc; beta=β)
        return -β^2 * Float64(dε_dβ), Float64(c_v)
    end,
)

"""
    DEFAULT_IDENTITIES

The identity set evaluated by `verify_thermodynamic_identities` when
the caller does not pass an explicit `identities` kwarg.
"""
const DEFAULT_IDENTITIES = ThermoIdentity[GIBBS_RELATION, SPECIFIC_HEAT_FROM_ENERGY]

# ──────────────────────────────────────────────────────────────────────
# Dispatch-existence helper
# ──────────────────────────────────────────────────────────────────────

# Capture the catch-all `fetch(::AbstractQAtlasModel, ::AbstractQuantity,
# ::BoundaryCondition; ...)` once at file load; any (model, quantity, bc)
# triple whose `which(fetch, ...)` returns this same Method object is
# *not* implemented (the catch-all just throws an informative error).
const _CATCH_ALL_FETCH_METHOD = which(
    fetch, Tuple{AbstractQAtlasModel,AbstractQuantity,BoundaryCondition}
)

function _has_dispatch(model, ::Type{Q}, bc) where {Q}
    return which(fetch, Tuple{typeof(model),Q,typeof(bc)}) !== _CATCH_ALL_FETCH_METHOD
end

function _can_run(identity::ThermoIdentity, model, bc)
    all(Q -> _has_dispatch(model, Q, bc), identity.requires)
end

# ──────────────────────────────────────────────────────────────────────
# Harness
# ──────────────────────────────────────────────────────────────────────

"""
    verify_thermodynamic_identities(model, bc;
                                     βs,
                                     identities=DEFAULT_IDENTITIES,
                                     rtol=1e-8, atol=1e-10)
        -> Vector{IdentityCheckResult}

For every `(identity, β)` pair, evaluate `identity.check(model, bc, (;β))`
and record whether `lhs ≈ rhs` within `(rtol, atol)`.  Identities that
require a quantity not dispatchable on `(model, bc)` are recorded as
`:skipped` (numeric fields `NaN`).

The skip semantics let the same call be made on any model — including
ones missing some of the registry — without spurious test failures.
This is what makes the harness genuinely model-agnostic.

Use from `@testset`s with `@test all(r.status === :pass for r in results)`
when every required quantity is implemented (TFIM today), or
`@test all(r.status !== :fail for r in results)` for partially-populated
models.
"""
function verify_thermodynamic_identities(
    model::AbstractQAtlasModel,
    bc::BoundaryCondition;
    βs::AbstractVector{<:Real},
    identities::AbstractVector{ThermoIdentity}=DEFAULT_IDENTITIES,
    rtol::Real=1e-8,
    atol::Real=1e-10,
)
    results = IdentityCheckResult[]
    for identity in identities
        runnable = _can_run(identity, model, bc)
        for β in βs
            params = (; β=β)
            if !runnable
                push!(
                    results,
                    IdentityCheckResult(
                        model, bc, identity.name, params, NaN, NaN, NaN, NaN, :skipped
                    ),
                )
                continue
            end
            lhs, rhs = identity.check(model, bc, params)
            abs_err = abs(lhs - rhs)
            rel_err = abs_err / max(abs(lhs), abs(rhs), eps())
            status = (abs_err ≤ atol || rel_err ≤ rtol) ? :pass : :fail
            push!(
                results,
                IdentityCheckResult(
                    model, bc, identity.name, params, lhs, rhs, abs_err, rel_err, status
                ),
            )
        end
    end
    return results
end

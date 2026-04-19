# deprecate/legacy_fetch.jl — Symbol-dispatch shim retained for
# backward compatibility with pre-v0.13 call sites.
#
# This file is loaded LAST by src/QAtlas.jl so it can route a legacy
# call into any of the already-registered concrete-struct `fetch`
# methods.  Scheduled for removal in v1.0 (see src/deprecate/README.md).

# The canonical Symbol-dispatch shim previously lived in
# src/core/type.jl.  It is redefined here so the legacy surface is
# centralised in one place and emits a one-shot deprecation hint.

@doc """
    fetch(m::Symbol, q::Symbol, bc::BoundaryCondition = Infinite(); kwargs...)

**Deprecated in v0.13.** Prefer the concrete-struct form:

    fetch(TFIM(; J=1.0, h=1.0), Energy(), OBC(N=24); beta=5.0)

This symbol-based signature is kept for drop-in backward compatibility
— the call is routed through `Model(m; kwargs...)` + `Quantity(q)` and
forwarded to the concrete `fetch(::AbstractQAtlasModel, ::AbstractQuantity, ::BoundaryCondition)`
method registered by each model.  A one-shot informational log is
emitted the first time the shim is reached from a given site.
""" fetch(::Symbol, ::Symbol, ::BoundaryCondition)

function fetch(m::Symbol, q::Symbol, bc::BoundaryCondition=Infinite(); kwargs...)
    @info """QAtlas symbol-dispatch `fetch(:$(m), :$(q), …)` will be \
removed in v1.0; prefer the concrete-struct API
    (e.g. `fetch(TFIM(; J, h), Energy(), OBC(N=…); kwargs...)`).""" maxlog = 1 _id = (
        :qatlas_legacy_fetch_symbol, m, q
    )
    return fetch(Model(m; kwargs...), Quantity(q), bc; kwargs...)
end

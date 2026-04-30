# core/registry.jl — declarative implementation registry.
#
# Each (model, quantity, bc) triple that has a directly-implemented
# `fetch` method also gets a one-liner `@register` declaration next to
# it, capturing the human metadata that `methods(fetch)` cannot:
#   * `method`      — algorithm tag (`:bdg`, `:dense_ed`, `:analytic`,
#                      `:transfer_matrix`, `:bethe_ansatz`, `:tba`, `:pfaffian`,
#                      `:not_implemented`, …)
#   * `reliability` — `:high` (closed-form + literature-tested),
#                      `:medium` (ED only / cross-check),
#                      `:low` (heuristic, not validated),
#                      `:not_implemented`.  Aligned with the
#                      (a)/(b)/(c) test categories in #118.
#   * `tested_in`   — relative path to the test file that validates
#                      this triple (or `nothing` if no dedicated test).
#   * `references`  — short literature pointers (author + year).
#   * `notes`       — caller-facing caveats (granularity, kwargs, …).
#
# `implementation_status()` returns the registry as
# `Vector{NamedTuple}`, which is `Tables.jl`-compatible without us
# taking a Tables dependency — downstream users wanting a `DataFrame`
# can call `DataFrame(implementation_status())` themselves.
#
# Conversion fallbacks (e.g. `Energy(:per_site)` at OBC routed through
# the `Energy(:total)` native + `÷ N`) are *not* registered separately:
# the registry reflects native implementations, and the routing is
# automatic by design.

"""
    Implementation

A single `(model, quantity, bc)` row of the QAtlas implementation
registry.  See [`@register`](@ref) for how rows are added and
[`implementation_status`](@ref) for how to query them.
"""
struct Implementation
    model::Type
    quantity::Type
    bc::Type
    method::Symbol
    reliability::Symbol
    tested_in::Union{String,Nothing}
    references::Vector{String}
    notes::String
end

"""
    REGISTRY :: Vector{Implementation}

Module-level mutable vector populated at include-time by `@register`
calls scattered across `src/models/.../<Model>_registry.jl` files.
Public read API: [`implementation_status`](@ref).
"""
const REGISTRY = Implementation[]

"""
    register!(model_T, quantity_T, bc_T;
              method=:unknown, reliability=:unknown,
              tested_in=nothing, references=String[], notes="")

Push a new [`Implementation`](@ref) row into [`REGISTRY`](@ref).
Usually called via the [`@register`](@ref) macro for ergonomics.
"""
function register!(
    model_T::Type, quantity_T::Type, bc_T::Type;
    method::Symbol=:unknown,
    reliability::Symbol=:unknown,
    tested_in::Union{String,Nothing}=nothing,
    references::AbstractVector{<:AbstractString}=String[],
    notes::AbstractString="",
)
    push!(
        REGISTRY,
        Implementation(
            model_T, quantity_T, bc_T,
            method, reliability,
            tested_in, String[r for r in references], String(notes),
        ),
    )
    return nothing
end

"""
    @register Model Quantity BC method=… reliability=… tested_in=… references=… notes=…

Thin macro around [`register!`](@ref).  Lets each model file register
its native fetch methods declaratively, e.g.

```julia
@register TFIM Energy{:total} OBC method=:bdg reliability=:high \\
    tested_in="test/models/test_TFIM_thermal.jl" \\
    references=["Pfeuty 1970"]
```

The three positional arguments are spliced as types; the remaining
`key=value` pairs are forwarded as keyword arguments to
[`register!`](@ref).
"""
macro register(model_T, quantity_T, bc_T, kwargs...)
    kw_exprs = map(kwargs) do kw
        kw isa Expr && kw.head === :(=) ||
            error("@register: expected key=value, got $kw")
        return Expr(:kw, kw.args[1], esc(kw.args[2]))
    end
    return Expr(
        :call, register!,
        Expr(:parameters, kw_exprs...),
        esc(model_T), esc(quantity_T), esc(bc_T),
    )
end

# ──────────────────────────────────────────────────────────────────────
# Smoke-check helper: distinguish "registered + fetch method exists" from
# "registered but only the catch-all error method matches".  The
# catch-all lives in core/type.jl; we capture it here at registry load
# time so future fetch additions can't accidentally re-shadow it.
# ──────────────────────────────────────────────────────────────────────

const _CATCH_ALL_FETCH_METHOD = which(
    fetch, Tuple{AbstractQAtlasModel,AbstractQuantity,BoundaryCondition}
)

"""
    has_native_fetch(impl::Implementation) -> Bool

`true` iff `which(fetch, (impl.model, impl.quantity, impl.bc))` resolves
to a method *more specific* than the catch-all in `core/type.jl`.
Conversion fallbacks (e.g. the generic `Energy{:per_site}` ↔
`Energy{:total}` router) count as "native" because they are a real
dispatchable implementation — they just live above the model layer.

Used by `test/core/test_registry.jl` to detect registry rows that
silently lost their backing fetch method.
"""
function has_native_fetch(impl::Implementation)
    m = which(fetch, Tuple{impl.model,impl.quantity,impl.bc})
    return m !== _CATCH_ALL_FETCH_METHOD
end

# ──────────────────────────────────────────────────────────────────────
# Query API
# ──────────────────────────────────────────────────────────────────────

_to_nt(e::Implementation) = (
    model=e.model, quantity=e.quantity, bc=e.bc,
    method=e.method, reliability=e.reliability,
    tested_in=e.tested_in, references=e.references, notes=e.notes,
)

"""
    implementation_status() -> Vector{NamedTuple}
    implementation_status(model::AbstractQAtlasModel)
    implementation_status(::Type{<:AbstractQAtlasModel})
    implementation_status(quantity::AbstractQuantity)
    implementation_status(::Type{<:AbstractQuantity})
    implementation_status(queue::AbstractVector)

Return registry rows as `NamedTuple`s (`Tables.jl`-compatible without a
Tables dependency).

- No-arg: every registered triple.
- `model` / `quantity` (instance or type): rows whose corresponding type
  field matches exactly (no subtype walking — model parameters are part
  of the identity here).
- `queue`: a vector of `(model, quantity, bc)` triples (each component
  may be either an instance or a type).  Returns one row per queue
  entry that is registered, dropping entries that are not.

Use this to plan downstream work — e.g. before writing tests for a new
ThermalMPS workload, query the queue you intend to validate against.
"""
implementation_status() = [_to_nt(e) for e in REGISTRY]

implementation_status(::Type{M}) where {M<:AbstractQAtlasModel} =
    [_to_nt(e) for e in REGISTRY if e.model === M]
implementation_status(model::AbstractQAtlasModel) = implementation_status(typeof(model))

implementation_status(::Type{Q}) where {Q<:AbstractQuantity} =
    [_to_nt(e) for e in REGISTRY if e.quantity === Q]
implementation_status(quantity::AbstractQuantity) = implementation_status(typeof(quantity))

function implementation_status(queue::AbstractVector)
    out = NamedTuple[]
    for q in queue
        length(q) == 3 || error(
            "implementation_status(queue): each element must be a " *
            "(model, quantity, bc) triple; got length $(length(q))",
        )
        m_T = q[1] isa Type ? q[1] : typeof(q[1])
        q_T = q[2] isa Type ? q[2] : typeof(q[2])
        bc_T = q[3] isa Type ? q[3] : typeof(q[3])
        for e in REGISTRY
            if e.model === m_T && e.quantity === q_T && e.bc === bc_T
                push!(out, _to_nt(e))
                break
            end
        end
    end
    return out
end

# ──────────────────────────────────────────────────────────────────────
# Markdown rendering
# ──────────────────────────────────────────────────────────────────────

# Strip the leading `QAtlas.` (and any submodule prefix) from a `Type`
# so the rendered table reads `TFIM`, `Energy{:total}`, `OBC` rather
# than `QAtlas.TFIM` etc.
function _short_type(T::Type)
    s = string(T)
    return replace(s, r"^QAtlas\.|Main\.QAtlas\." => "")
end

"""
    implementation_status_markdown([io::IO=stdout], entries=implementation_status())

Render `entries` (any iterable of `NamedTuple` rows from
[`implementation_status`](@ref)) as a GitHub-flavoured Markdown table
to `io`.
"""
function implementation_status_markdown(
    io::IO=stdout, entries=implementation_status()
)
    println(io, "| Model | Quantity | BC | Method | Reliability | Tested in | References |")
    println(io, "|---|---|---|---|---|---|---|")
    for e in entries
        println(
            io,
            "| ", _short_type(e.model),
            " | ", _short_type(e.quantity),
            " | ", _short_type(e.bc),
            " | `", e.method, "`",
            " | `", e.reliability, "`",
            " | ", something(e.tested_in, "—"),
            " | ", isempty(e.references) ? "—" : join(e.references, "; "),
            " |",
        )
    end
    return nothing
end

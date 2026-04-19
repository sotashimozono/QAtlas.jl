# src/deprecate/

Shims keeping the pre-v0.13 API alive.  Everything in this directory is
scheduled for removal in **v1.0**; a `git rm -r src/deprecate/` at that
time should cleanly retire the legacy surface with no collateral damage.

## Contents

| File | Purpose |
|---|---|
| `legacy_fetch.jl` | `fetch(m::Symbol, q::Symbol, bc; kwargs...)` shim + one-shot deprecation warning. |

## Include order

`src/QAtlas.jl` loads `deprecate/*` **last**, after every model struct
and every `fetch` method has been registered.  This lets the shim route
a legacy symbol call into the new concrete-struct dispatch.

## Policy

- New files here must be self-contained and import-only (they reference
  public types from `core/` + model files, nothing else).
- No new public API lives here — only back-compat wrappers for names
  that existed in v0.12 and earlier.

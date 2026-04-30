# models/quantum/TFIM/TFIM_registry.jl — declarative implementation map.
#
# One `@register` line per natively-implemented (model, quantity, bc)
# triple.  Conversion fallbacks (e.g. `Energy(:per_site)` at OBC routed
# through the `Energy(:total)` native + `÷ N`) are *not* listed here:
# the registry tracks native implementations and the routing is
# automatic.  See `src/core/registry.jl` for the metadata schema.
#
# When you add a new fetch method to TFIM, add a sibling `@register`
# line here.  `test/core/test_registry.jl` will fail loudly if the
# registry and the dispatch table drift apart.

# ── Energy (granularity-aware) ─────────────────────────────────────────
@register(
    TFIM, Energy{:total}, OBC,
    method=:bdg, reliability=:high,
    tested_in="test/models/test_TFIM_thermal.jl",
    references=["Pfeuty 1970"],
    notes="Total ⟨H⟩(β) via the BdG spectrum; ground state when no β kwarg.",
)
@register(
    TFIM, Energy{:per_site}, Infinite,
    method=:bdg, reliability=:high,
    tested_in="test/models/test_TFIM_thermal.jl",
    references=["Pfeuty 1970"],
    notes="Per-site ε(β) by QuadGK over the PBC dispersion Λ(k).",
)

# ── Spectrum / criticality ────────────────────────────────────────────
@register(
    TFIM, MassGap, Infinite,
    method=:analytic, reliability=:high,
    tested_in="test/models/test_TFIM_massgap.jl",
    references=["Pfeuty 1970"],
    notes="Δ_∞(J,h) = 2|h - J| — closed-form Ising gap.",
)
@register(
    TFIM, MassGap, OBC,
    method=:bdg, reliability=:high,
    tested_in="test/models/test_TFIM_massgap.jl",
    references=["Pfeuty 1970"],
    notes="Smallest positive BdG eigenvalue of the OBC chain.",
)
@register(
    TFIM, CentralCharge, Infinite,
    method=:analytic, reliability=:high,
    tested_in="test/models/test_TFIM_central_charge.jl",
    references=["Belavin-Polyakov-Zamolodchikov 1984"],
    notes="c = 1/2 at the critical point (h = J), 0 otherwise.",
)

# ── Free-fermion thermal (per-site) — meta-defined in TFIM_thermal.jl ──
@register(TFIM, FreeEnergy,       OBC,      method=:bdg, reliability=:high, tested_in="test/models/test_TFIM_thermal.jl")
@register(TFIM, FreeEnergy,       Infinite, method=:bdg, reliability=:high, tested_in="test/models/test_TFIM_thermal.jl")
@register(TFIM, ThermalEntropy,   OBC,      method=:bdg, reliability=:high, tested_in="test/models/test_TFIM_thermal.jl")
@register(TFIM, ThermalEntropy,   Infinite, method=:bdg, reliability=:high, tested_in="test/models/test_TFIM_thermal.jl")
@register(TFIM, SpecificHeat,     OBC,      method=:bdg, reliability=:high, tested_in="test/models/test_TFIM_thermal.jl")
@register(TFIM, SpecificHeat,     Infinite, method=:bdg, reliability=:high, tested_in="test/models/test_TFIM_thermal.jl")
@register(TFIM, MagnetizationX,   OBC,      method=:bdg, reliability=:high, tested_in="test/models/test_TFIM_thermal.jl")
@register(TFIM, MagnetizationX,   Infinite, method=:bdg, reliability=:high, tested_in="test/models/test_TFIM_thermal.jl")
@register(TFIM, SusceptibilityXX, OBC,      method=:bdg, reliability=:high, tested_in="test/models/test_TFIM_thermal.jl")
@register(TFIM, SusceptibilityXX, Infinite, method=:bdg, reliability=:high, tested_in="test/models/test_TFIM_thermal.jl")

# ── Local one-site observables (per-site index, not bulk-averaged) ────
@register(TFIM, EnergyLocal,         OBC, method=:bdg, reliability=:high, tested_in="test/models/test_TFIM_local.jl")
@register(TFIM, MagnetizationXLocal, OBC, method=:bdg, reliability=:high, tested_in="test/models/test_TFIM_local.jl")
@register(TFIM, MagnetizationZLocal, OBC, method=:bdg, reliability=:high, tested_in="test/models/test_TFIM_local.jl")

# ── Z-axis dynamics + correlations (BdG time evolution) ───────────────
@register(TFIM, SusceptibilityZZ,    OBC, method=:bdg, reliability=:high, tested_in="test/models/test_TFIM_dynamics.jl")
@register(TFIM, ZZStructureFactor,   OBC, method=:bdg, reliability=:high, tested_in="test/models/test_TFIM_dynamics.jl")
@register(TFIM, ZZCorrelation{:static}, OBC, method=:bdg, reliability=:high, tested_in="test/models/test_TFIM_dynamics.jl")
@register(
    TFIM, ZZCorrelation{:dynamic}, OBC,
    method=:bdg, reliability=:high,
    tested_in="test/models/test_TFIM_dynamics.jl",
    notes="Single (i, j, t) point; for sweeps loop the kwargs.",
)
@register(
    TFIM, ZZCorrelation{:lightcone}, OBC,
    method=:bdg, reliability=:high,
    tested_in="test/models/test_TFIM_dynamics.jl",
    notes="C(r,t) lightcone slice for fixed center; takes a `times` vector.",
)
@register(TFIM, XXCorrelation{:dynamic}, OBC, method=:bdg, reliability=:high, tested_in="test/models/test_TFIM_dynamics.jl")

# ── Entanglement (T = 0; β kwarg defaults to Inf) ─────────────────────
@register(
    TFIM, VonNeumannEntropy, OBC,
    method=:bdg, reliability=:high,
    tested_in="test/models/test_TFIM_entanglement.jl",
    references=["Peschel 2003"],
    notes="Free-fermion correlation-matrix method; pass subsystem length ℓ.",
)

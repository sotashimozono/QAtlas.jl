# Contributing to QAtlas.jl

Thank you for your interest in contributing to QAtlas! This guide will help you understand how the project is organized and what we look for in contributions.

## What is QAtlas?

QAtlas is a **dictionary of rigorous results in quantum and statistical physics**. It stores analytically known exact values in `src/` and verifies them against independent numerical calculations in `test/`. The core value proposition is that every stored result is cross-validated from at least two independent theoretical or computational sources.

## Design Principles

### `src/` is a leaf — no lattice-package dependencies

The source code in `src/` does not depend on Lattice2D, QuasiCrystal, or any lattice construction package. It contains pure functions that map `(model, quantity, bc) → value`. Lattice-package dependencies live exclusively in `test/` via `[extras]`.

### Accumulate results first, refactor later

Adding a new rigorous result is always more valuable than perfecting the code structure. If a new result doesn't fit cleanly into the existing directory layout, **add it anyway and verify it** — the structure can be refactored later without losing the verified result.

### Physical correctness is paramount

A value in `src/` is only considered rigorous once it has been **independently verified** in `test/`. Internal consistency checks (e.g., scaling relations) are necessary but not sufficient. We require **cross-verification from independent theoretical sources** whenever possible.

### Every numerical value traces to a derivation

For each new rigorous result, the accompanying `docs/src/calc/` derivation must be **complete and step-by-step** (see [md/docs-conventions.md](md/docs-conventions.md) for the depth standard). Forbidden phrases such as "it can be shown" / "we omit details" / "standard calculation gives" must not appear.

## Repository Structure

```text
src/
├── core/                              # Type hierarchy, aliases, BC / Quantity structs
│   ├── type.jl                        #   AbstractQAtlasModel, OBC(N), PBC(N), Infinite, Quantity{S}
│   ├── alias.jl                       #   Symbol -> canonical-name dispatch table
│   └── quantities.jl                  #   MagnetizationX, SusceptibilityZZ, ZZCorrelation{M}, ...
├── deprecate/                         # Pre-v0.13 Symbol-dispatch shims (removable at v1.0)
├── universalities/                    # Universality{C}(...) parametric type
│   ├── Universality.jl                #   base type + exponent table
│   ├── Ising2D.jl                     #   2D Ising (exact + 3D numerical bootstrap)
│   ├── KPZ.jl, Percolation.jl, ...
│   └── E8.jl                          #   E8 mass ratios
├── universalities/ONModel.jl          #   XY / Heisenberg / O(N) nonlinear σ
└── models/                            # Layout: `<class>/<Model>/<Model>.jl`.
    │                                  # Per-model directory so each Model has
    │                                  # its own dir even if currently single
    │                                  # file (room to grow into thermal /
    │                                  # dynamics / local siblings as needed).
    ├── classical/
    │   └── IsingSquare/
    │       └── IsingSquare.jl         #   IsingSquare(; J, Lx, Ly)
    └── quantum/
        ├── TFIM/
        │   ├── TFIM.jl                #   TFIM(; J, h) + Energy + MassGap + CentralCharge
        │   ├── TFIM_thermal.jl        #   FreeEnergy, ThermalEntropy, SpecificHeat, M_x, χ_xx
        │   ├── TFIM_dynamics.jl       #   Majorana covariance, {XX,YY,ZZ}Correlation{M}, structure factors
        │   └── TFIM_local.jl          #   MagnetizationXLocal, MagnetizationZLocal, EnergyLocal
        ├── Heisenberg/
        │   └── Heisenberg.jl          #   Heisenberg1D() + ExactSpectrum + GroundStateEnergyDensity
        ├── XXZ/
        │   └── XXZ.jl                 #   XXZ1D(; J, Δ) + Luttinger K, u, central charge, Energy
        └── tightbinding/              # family dir (multiple tight-binding lattices)
            └── regular/               # Bloch-diagonalisable periodic lattices
                ├── Honeycomb.jl       #   Honeycomb(; t, Lx, Ly)
                ├── Kagome.jl
                ├── Lieb.jl
                └── Triangular.jl
            # Future: tightbinding/{quasicrystalline,fractal,disordered}/
            # for Fibonacci, Penrose, Sierpinski, Anderson, ...

test/
├── util/                              # Reusable verification helpers (dense + sparse)
│   ├── classical_partition.jl         #   Brute-force partition function
│   ├── tight_binding.jl               #   Real-space TB Hamiltonian
│   ├── spinhalf_ed.jl                 #   Dense spin-1/2 many-body ED
│   ├── sparse_ed.jl                   #   Sparse ED + KrylovKit Lanczos ground state
│   └── bloch.jl                       #   Generic Bloch Hamiltonian builder
├── standalone/                        # Lattice-independent tests
└── verification/                      # Cross-validation against Lattice2D

docs/src/                              # Documenter + GitHub Pages
├── calc/                              # Step-by-step derivations (Zettelkasten)
├── models/, universalities/, methods/ # API + results + narrative (links back to calc/)
└── verification/                      # Cross-check narrative

md/                                    # Dev memos (Japanese OK, not published)
├── api.md, roadmap.md
├── docs-conventions.md                # Docs depth standard — READ THIS BEFORE WRITING calc/
└── models/, universalities/           # Per-topic implementation notes
```

## The v0.13 API

Every model is a **concrete struct** that carries its physical parameters as typed fields. Quantities are likewise concrete structs (not `Symbol`s). Boundary conditions carry their own size.

```julia
# Canonical 3-dispatch form
fetch(TFIM(; J=1.0, h=1.0), Energy(), OBC(24); beta=5.0)
fetch(XXZ1D(; J=1.0, Δ=0.5), LuttingerParameter(), Infinite())
fetch(IsingSquare(; J=1.0, Lx=4, Ly=4), PartitionFunction(); β=0.44)
```

Legacy Symbol calls `fetch(:TFIM, :energy, OBC(); N=24, J=1.0, h=1.0, beta=5.0)` still work through the thin shim layer in `src/deprecate/` but emit a deprecation log.

See [md/api.md](md/api.md) for the full v0.13 design rationale (type hierarchy, quantity naming, deprecation path).

## How to Contribute

### Adding a new rigorous result

1. **Define a concrete model struct** (if the model does not yet exist). Physical parameters go in typed fields; boundary conditions and lattice size use the `OBC(N)` / `PBC(N)` / `Infinite` types — never a loose kwarg.

   ```julia
   struct MyModel <: QAtlas.AbstractQAtlasModel
       J::Float64
       h::Float64
   end
   MyModel(; J::Real=1.0, h::Real=0.0) = MyModel(Float64(J), Float64(h))

   function QAtlas.fetch(m::MyModel, ::Energy, bc::OBC; kwargs...)
       N = QAtlas._bc_size(bc, kwargs)
       ...
   end
   ```

2. **Cite the source precisely** in the docstring — not just "Author (Year)" but the specific equation, table, or theorem number. The `md/docs-conventions.md` policy applies equally to inline docstrings: no "it can be shown".

   ```julia
   # Good: traceable to a single line in the literature
   """
       fetch(model::IsingSquare, ::SpontaneousMagnetization) -> Float64

   Yang (1952) Phys. Rev. 85, 808 — from the spontaneous magnetization
   formula M(T) = (1 − sinh⁻⁴(2βJ))^{1/8}.
   """
   ```

3. **Write a derivation note under `docs/src/calc/`** covering the exact computation the new result is based on. Follow the Main result / Setup / Calculation / References / Used by skeleton of [md/docs-conventions.md](md/docs-conventions.md). Non-trivial algebra must be shown line by line; integrals must be evaluated (contour, poles, residues), not cited as "standard calculation gives".

4. **Link the derivation** in the corresponding `docs/src/models/…` or `docs/src/universalities/…` page via `See full derivation: [...](../../calc/…md)`.

5. **Write tests that independently verify the result.** The project uses a three-tier structure:

   | Tier               | Location             | Purpose                                                       |
   | ------------------ | -------------------- | ------------------------------------------------------------- |
   | Standalone         | `test/standalone/`   | Lattice-independent checks: special values, scaling relations |
   | Verification       | `test/verification/` | Cross-check against Lattice2D-based calculations              |
   | Cross-verification | `test/verification/` | Connect universality exponents to model-specific results      |

6. **Include at least one cross-verification** that links your result to an existing QAtlas result derived from a **different theoretical source**. This is what makes QAtlas rigorous.

### Adding a universality class

Use the `Universality{C}` parametric type with dimension `d` as a keyword:

```julia
fetch(Universality(:MyClass), CriticalExponents(); d=2)
```

- Use `Rational{Int}` for exact values (e.g. `β = 1//8`).
- For numerical estimates, include `_err` fields for uncertainty (e.g. `β = 0.32642, β_err = 0.00001`).
- Verify scaling relations in tests: `α + 2β + γ = 2` (Rushbrooke), `γ = β(δ−1)` (Widom), etc.

### Adding a tight-binding lattice

The generic Bloch builder in `test/util/bloch.jl` works for **any** Lattice2D topology. For lattices already supported by Lattice2D, a test file is usually all that is needed:

```julia
λ_bloch = bloch_tb_spectrum(MyTopology, Lx, Ly, t)
H = build_tight_binding(lat, t)
λ_real = sort(eigvals(Symmetric(H)))
@test λ_bloch ≈ λ_real atol=1e-10
```

If you want a hardcoded Bloch formula in `src/`, place it in `models/quantum/tightbinding/` as a concrete struct:

```julia
struct MyLattice <: QAtlas.AbstractQAtlasModel
    t::Float64
    Lx::Int
    Ly::Int
end
```

**Topology-name collisions with Lattice2D.** The identifiers `Honeycomb`, `Kagome`, `Lieb`, `Triangular` exist in both QAtlas and Lattice2D. To avoid `UndefVarError` in code that `using`s both packages, QAtlas **does not export** any of these names. Always qualify them as `QAtlas.Honeycomb()`, `QAtlas.Kagome()`, etc. The backward-compat alias `Graphene = Honeycomb` **is** exported because the name does not collide.

## Documentation

Writing a new `docs/src/calc/` note? First read [md/docs-conventions.md](md/docs-conventions.md). The depth standard is enforced by:

- grep check: `docs/src/calc/<new-note>.md` must have zero matches for `it can be shown`, `we omit`, `standard calculation`, `as in [Author] [Journal]`, `one can verify`, `it is easy to see`, `it follows immediately`.
- Structure check: `## Main result`, `## Setup`, `## Calculation`, `## References`, `## Used by` sections, in that order.
- Documenter build: `julia --project=docs -e 'include("docs/make.jl")'` must complete with no `Error:` lines.

The exemplars for depth are [docs/src/calc/jw-tfim-bdg.md](docs/src/calc/jw-tfim-bdg.md) and [docs/src/calc/bethe-ansatz-heisenberg-e0.md](docs/src/calc/bethe-ansatz-heisenberg-e0.md).

Public docs (`docs/src/`) are **English only**. Dev memos under `md/` may be in Japanese.

## Things to Watch Out For

- **`QAtlas.fetch` vs `Base.fetch`**: always qualify as `QAtlas.fetch(...)` in test code to avoid the name collision with Julia's built-in `Base.fetch`.
- **Bond counting in small PBC systems**: Lattice2D's `bonds(lat)` double-counts bonds when `Lx = 2` or `Ly = 2` with periodic boundaries. Both the transfer-matrix and brute-force paths use the same convention, so they agree — but be aware of this when writing new models.
- **OBC vs PBC for gap analysis**: in the ordered phase of the TFIM (h ≪ J), the lowest ED "gap" is actually the Z₂ tunneling splitting (exponentially small in N), not the physical excitation gap. This is correct physics but can be surprising in tests.
- **ForwardDiff compatibility**: if you want a `fetch` function to support automatic differentiation, relax type constraints from `Float64` to `Real` and avoid LAPACK-only operations like `eigvals(Symmetric(T))` — use `tr(T^n)` instead, which works with dual numbers.
- **Entanglement tests and RAM**: the ED-based entanglement-entropy tests (`test_entanglement_central_charge.jl`) scale as $O(2^N)$ memory at $N \ge 14$. The default PR-CI profile runs $N \le 14$ via sparse + KrylovKit Lanczos. `QATLAS_TEST_FULL=1` enables $N = 16$ (~24 GB peak).

## Before Submitting a PR

1. Run the full test suite with full threading:

   ```bash
   julia -e 'using Pkg; Pkg.activate("."); Pkg.test(; julia_args=`-t auto --heap-size-hint=96G`)'
   ```

   On CI (GitHub Actions) this is automatic via `julia-args: "--threads=auto"` in `.github/workflows/CI.yml`.

2. If your PR adds new features, update the version in `Project.toml` (patch bump per PR; minor for larger additions) — see [md/roadmap.md](md/roadmap.md) for the current version-bump policy.

3. If you added a new `docs/src/calc/` derivation, build the docs locally and confirm:

   ```bash
   julia --project=docs -e 'include("docs/make.jl")'
   ```

   completes with exit 0 and no `Error:` lines in stderr.

4. If you find a discrepancy between QAtlas values and the literature, please open an issue with the full reference (journal, volume, page, equation number) so we can investigate.

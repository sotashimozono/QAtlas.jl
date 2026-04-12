# Contributing to QAtlas.jl

Thank you for your interest in contributing to QAtlas! This guide will help you understand how the project is organized and what we look for in contributions.

## What is QAtlas?

QAtlas is a **dictionary of rigorous results in quantum and statistical physics**. It stores analytically known exact values in `src/` and verifies them against independent numerical calculations in `test/`. The core value proposition is that every stored result is cross-validated from at least two independent theoretical or computational sources.

## Design Principles

### `src/` is a leaf — no lattice-package dependencies

The source code in `src/` does not depend on Lattice2D, QuasiCrystal, or any lattice construction package. It contains pure functions that map `(model, quantity) → value`. Lattice-package dependencies live exclusively in `test/` via `[extras]`.

### Accumulate results first, refactor later

Adding a new rigorous result is always more valuable than perfecting the code structure. If a new result doesn't fit cleanly into the existing directory layout, **add it anyway and verify it** — the structure can be refactored later without losing the verified result.

### Physical correctness is paramount

A value in `src/` is only considered rigorous once it has been **independently verified** in `test/`. Internal consistency checks (e.g., scaling relations) are necessary but not sufficient. We require **cross-verification from independent theoretical sources** whenever possible.

## Repository Structure

```text
src/
├── core/                          # Type definitions, aliases, utilities
├── universalities/                # Universality class exponents
│   ├── Universality.jl            # Universality{C} parametric type
│   ├── Ising2D.jl                 # 2D Ising (exact + 3D numerical)
│   ├── KPZ.jl, Percolation.jl     # Other universality classes
│   └── E8.jl                      # E8 mass ratios
└── models/
    ├── classical/                 # IsingSquare (Z, T_c, M(T))
    └── quantum/
        ├── tightbinding/          # Graphene, Kagome, Lieb, Triangular
        └── Heisenberg.jl          # Dimer + 4-site PBC

test/
├── util/                          # Reusable verification helpers
│   ├── classical_partition.jl     # Brute-force partition function
│   ├── tight_binding.jl           # Real-space TB Hamiltonian
│   ├── spinhalf_ed.jl             # Spin-1/2 many-body ED
│   └── bloch.jl                   # Generic Bloch Hamiltonian builder
├── standalone/                    # Lattice-independent tests
└── verification/                  # Cross-validation against Lattice2D
```

## How to Contribute

### Adding a new rigorous result

1. **Write the result in `src/`** using a dispatch tag and `fetch`:

   ```julia
   struct MyModel end
   struct MyQuantity end
   fetch(::MyModel, ::MyQuantity; params...) = ...
   ```

2. **Cite your sources precisely** — not just "Author (Year)" but the specific equation, table, or theorem number:

   ```julia
   # Good: traceable to a single line in the literature
   # β = 1/8: Yang (1952) Phys. Rev. 85, 808 — from the spontaneous
   #          magnetization formula M(T) = (1−sinh⁻⁴(2βJ))^{1/8}.

   # Insufficient: ambiguous origin
   # References: Onsager (1944).
   ```

3. **Write tests that independently verify the result.** We use a three-tier test structure:

   | Tier               | Location               | Purpose                                                       |
   | ------------------ | ---------------------- | ------------------------------------------------------------- |
   | Standalone         | `test/standalone/`     | Lattice-independent checks: special values, scaling relations |
   | Verification       | `test/verification/`   | Cross-check against Lattice2D-based calculations              |
   | Cross-verification | `test/verification/`   | Connect universality exponents to model-specific results      |

4. **Include at least one cross-verification** that links your result to an existing QAtlas result derived from a *different theoretical source*. This is what makes QAtlas rigorous.

### Adding a universality class

Use the `Universality{C}` parametric type with dimension `d` as a keyword:

```julia
fetch(Universality(:MyClass), CriticalExponents(); d=2)
```

- Use `Rational{Int}` for exact values (e.g., `β = 1//8`).
- For numerical estimates, include `_err` fields for uncertainty (e.g., `β = 0.32642, β_err = 0.00001`).
- Verify scaling relations in tests: `α + 2β + γ = 2` (Rushbrooke), `γ = β(δ−1)` (Widom), etc.

### Adding a tight-binding lattice

The generic Bloch builder in `test/util/bloch.jl` works for **any** Lattice2D topology. For lattices already supported by Lattice2D, you only need a test file:

```julia
λ_bloch = bloch_tb_spectrum(MyTopology, Lx, Ly, t)
H = build_tight_binding(lat, t)
λ_real = sort(eigvals(Symmetric(H)))
@test λ_bloch ≈ λ_real atol=1e-10
```

If you want to add a hardcoded Bloch formula to `src/`, place it in `models/quantum/tightbinding/`. Note that topology names that collide with Lattice2D exports (e.g., `Kagome`, `Lieb`) should **not** be exported from QAtlas — access them as `QAtlas.Kagome()`.

## Things to Watch Out For

- **`QAtlas.fetch` vs `Base.fetch`**: Always qualify as `QAtlas.fetch(...)` in test code to avoid the name collision with Julia's built-in `Base.fetch`.

- **Bond counting in small PBC systems**: Lattice2D's `bonds(lat)` double-counts bonds when `Lx = 2` or `Ly = 2` with periodic boundaries. Both the transfer-matrix and brute-force paths use the same convention, so they agree — but be aware of this when writing new models.

- **OBC vs PBC for gap analysis**: In the ordered phase of the TFIM (h ≪ J), the lowest ED "gap" is actually the Z₂ tunneling splitting (exponentially small), not the physical excitation gap. This is correct physics but can be surprising in tests.

- **ForwardDiff compatibility**: If you want a `fetch` function to support automatic differentiation, relax type constraints from `Float64` to `Real` and avoid LAPACK-only operations like `eigvals(Symmetric(T))` — use `tr(T^n)` instead, which works with dual numbers.

## Before Submitting a PR

1. Run the full test suite and confirm it passes:
   ```bash
   julia --project=. -e 'using Pkg; Pkg.test()'
   ```

2. If your PR adds new features, update the version in `Project.toml`.

3. If you find a discrepancy between QAtlas values and the literature, please open an issue with the full reference (journal, volume, page, equation number).

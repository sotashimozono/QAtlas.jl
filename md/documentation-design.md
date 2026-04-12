# Documentation Design — Blueprint

## Goal

Build a **scientific reference documentation** for QAtlas where every
stored result is:

1. Precisely stated (formula + conditions)
2. Traceable to a specific publication (author, year, journal, equation)
3. Justified by a derivation sketch (enough to verify, not a textbook)
4. Connected to QAtlas's implementation (`fetch` call)
5. Cross-validated by independent tests (linked to test files)

This is not typical API documentation — it is a **curated scientific
reference** that happens to be implemented in Julia.

## Workflow

1. **This file** (`md/documentation-design.md`): agree on structure + templates
2. **Skeleton** in `docs/src/`: create empty pages with headers
3. **Content**: fill in each page, starting from the most established results
4. **Deploy**: Documenter.jl builds + GitHub Pages

## Tree Structure

```
docs/src/
├── index.md                              # Project overview, philosophy
│
├── models/
│   ├── index.md                          # Model catalog (table)
│   │
│   ├── classical/
│   │   ├── index.md                      # Classical models overview
│   │   └── ising-square.md              # Z, T_c, M(T)
│   │
│   └── quantum/
│       ├── index.md                      # Quantum models overview
│       ├── tfim.md                       # BdG, energy, gap, thermal
│       ├── heisenberg.md               # Dimer, 4-site, Bethe ansatz
│       └── tightbinding/
│           ├── index.md                  # TB overview + generic Bloch
│           ├── graphene.md              # Honeycomb
│           ├── kagome.md               # Flat band +2t
│           ├── lieb.md                 # Flat band E=0
│           └── triangular.md           # Frustration
│
├── universalities/
│   ├── index.md                          # What are universality classes?
│   ├── ising.md                         # d=2 exact, d=3 bootstrap, MF
│   ├── percolation.md                   # d=2 exact, d=3 MC
│   ├── potts.md                         # q=3, q=4
│   ├── kpz.md                          # Growth exponents
│   ├── on-models.md                    # XY, Heisenberg
│   ├── mean-field.md                   # Landau baseline
│   └── e8.md                          # Mass ratios
│
├── verification/
│   ├── index.md                          # Verification philosophy
│   ├── cross-checks.md                 # Universality ↔ model connections
│   ├── entanglement.md                 # c from S(l), Calabrese-Cardy
│   └── disordered.md                  # IRFP, random singlet
│
├── methods/
│   ├── index.md                          # Computational methods overview
│   ├── transfer-matrix.md              # Partition function via T^Lx
│   ├── bloch-hamiltonian.md           # Generic Bloch builder
│   ├── exact-diagonalization.md       # Spin-1/2 ED util
│   ├── automatic-differentiation.md   # ForwardDiff on Z
│   └── entanglement-entropy.md        # Reduced density matrix + vN
│
└── api/
    └── reference.md                      # Auto-generated API docs
```

## Page Templates

### Model Result Page (e.g., `ising-square.md`)

Each model page contains one or more **result cards**. A page can have
multiple results (e.g., IsingSquare has Z, T_c, M(T)).

```markdown
# Classical 2D Ising Model on the Square Lattice

## Overview

[1-paragraph physical description: what model, what regime, why important]

---

## Partition Function

### Statement

$$Z(L_x, L_y, \beta, J) = \mathrm{Tr}(T^{L_x})$$

where $T$ is the $2^{L_y} \times 2^{L_y}$ symmetric transfer matrix...

### Physical Context

- Valid for: finite $L_x \times L_y$ square lattice with PBC
- Regime: any $\beta \geq 0$, $J > 0$ (ferromagnetic)
- Special cases: $\beta = 0 \Rightarrow Z = 2^N$, ...

### Derivation Sketch

The transfer matrix encodes the Boltzmann weight of adding one row...
[Key steps, not a full derivation. Enough for a reader to verify.]

### References

- L. Onsager, "Crystal Statistics. I.", Phys. Rev. **65**, 117 (1944).
- B. M. McCoy and T. T. Wu, *The Two-Dimensional Ising Model*,
  Harvard University Press (1973), Ch. 2.

### QAtlas API

```julia
Z = QAtlas.fetch(IsingSquare(), PartitionFunction();
                 Lx=4, Ly=4, β=0.44, J=1.0)
```

### Verification

| Test file | Method | What is checked |
|-----------|--------|-----------------|
| `test_ising_2x2_classical.jl` | Brute-force 2^N | Z_TM ≈ Z_BF |
| `test_ising_ad_thermodynamics.jl` | ForwardDiff | ⟨E⟩, C_v, ⟨Σσσ⟩ |

### Connections

- [Critical Temperature](#critical-temperature): $T_c$ from the same model
- [Ising Universality](../universalities/ising.md): $\alpha = 0$ confirmed
  via specific heat from $Z$

---

## Critical Temperature

### Statement

$$T_c = \frac{2J}{\ln(1 + \sqrt{2})} \approx 2.269\,J$$

...

[Same structure repeats for each result card]
```

### Universality Class Page (e.g., `ising.md`)

```markdown
# Ising Universality Class

## Overview

The Ising universality class describes Z₂-symmetric second-order
phase transitions. In $d = 2$ it is exactly solved (Virasoro minimal
model M(3,4), central charge $c = 1/2$).

## Critical Exponents

### $d = 2$ (Exact)

| Exponent | Value | Derivation | Reference |
|----------|-------|------------|-----------|
| $\beta$  | $1/8$ | Yang magnetization | Yang (1952) Phys. Rev. **85**, 808 |
| $\nu$    | $1$   | Coulomb gas | den Nijs (1979) J. Phys. A **12**, 1857 |
| ...      | ...   | ... | ... |

**Scaling relations** (verified algebraically in QAtlas via Rational arithmetic):
- Rushbrooke: $\alpha + 2\beta + \gamma = 2$ ✓
- Widom: $\gamma = \beta(\delta - 1)$ ✓
- ...

### $d = 3$ (Conformal Bootstrap)

| Exponent | Value | Uncertainty | Reference |
|----------|-------|-------------|-----------|
| $\beta$ | 0.32642 | ±0.00001 | Kos+ (2016) JHEP **08**, 036, Table 2 |
| ...

### $d \geq 4$ (Mean-Field)

See [Mean-Field](mean-field.md).

## Cross-Verification

| Exponent | Extracted from | Method | Agreement |
|----------|---------------|--------|-----------|
| $\beta = 1/8$ | Yang $M(T)$ near $T_c$ | log-log slope | 1% |
| $\nu z = 1$ | TFIM gap $\Delta(N)$ | log-log regression | 5% |
| $c = 1/2$ | TFIM $S(l)$ | Calabrese-Cardy OBC | 10% (N=14) |

See [Cross-Checks](../verification/cross-checks.md) for details.

## QAtlas API

```julia
# Exact (Rational)
e = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=2)

# Numerical (Float64 + _err)
e = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=3)
```
```

### Methods Page (e.g., `transfer-matrix.md`)

```markdown
# Transfer Matrix Method

## What QAtlas Computes

`IsingSquare.PartitionFunction` uses the transfer matrix to compute
the exact partition function $Z = \mathrm{Tr}(T^{L_x})$ for finite
$L_x \times L_y$ PBC square lattices.

## Physical Origin

The transfer matrix formulation arises from rewriting the partition
function as a matrix product over rows...

[Derivation sketch with key equations]

## Implementation

$T$ is a $2^{L_y} \times 2^{L_y}$ real symmetric matrix.
$Z = \mathrm{tr}(T^{L_x})$ is computed via matrix exponentiation
(not eigendecomposition, for ForwardDiff compatibility).

### Why `tr(T^Lx)` instead of `eigvals`

The original implementation used `eigvals(Symmetric(T))` followed by
$Z = \sum_i \lambda_i^{L_x}$. This was changed to `tr(T^{L_x})` to
support `ForwardDiff.Dual` numbers for automatic differentiation,
since LAPACK's eigensolve does not accept dual inputs.

## References

- R. J. Baxter, *Exactly Solved Models in Statistical Mechanics*
  (1982), Ch. 7 (Transfer matrix formalism).

## Connections

- [IsingSquare](../models/classical/ising-square.md): uses this method
- [AD Verification](../verification/cross-checks.md): ForwardDiff on Z
```

## Content Priority

Phase 1 — most established results (already well-tested):
1. `models/classical/ising-square.md` (Z, T_c, M(T))
2. `universalities/ising.md` (d=2 exact + cross-checks)
3. `verification/index.md` (philosophy)
4. `methods/transfer-matrix.md`

Phase 2 — quantum models:
5. `models/quantum/tfim.md`
6. `models/quantum/heisenberg.md`
7. `methods/exact-diagonalization.md`

Phase 3 — tight-binding + remaining universalities:
8. `models/quantum/tightbinding/*.md`
9. `universalities/{percolation,potts,kpz,...}.md`
10. `methods/bloch-hamiltonian.md`

Phase 4 — advanced topics:
11. `verification/entanglement.md`
12. `verification/disordered.md`
13. `methods/entanglement-entropy.md`

## Open Questions

1. **methods/ depth**: How much derivation to include? Enough to
   verify (a few key equations + reference) or textbook-level?
   → Proposed: "enough to verify" + pointer to textbook for full derivation.

2. **Numerical vs exact distinction**: Should we visually distinguish
   exact results (Rational) from numerical estimates (Float64 + err)?
   → Proposed: use different formatting or a badge system.

3. **Versioning of references**: Some numerical values (3D Ising) are
   updated as better bootstrap bounds are published. How to track?
   → Proposed: each numerical value includes the paper year; we update
   when a significantly improved estimate is published.

4. **Connection to Documenter.jl**: The `api/reference.md` page can
   use `@autodocs` to pull docstrings. Other pages are hand-written.
   The `make.jl` needs to specify the page order matching our tree.

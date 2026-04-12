# Verification Philosophy

## Overview

QAtlas is a reference library: every stored value must be
**demonstrably correct**. Unlike typical numerical libraries where
testing checks that code runs without error, QAtlas's tests
constitute **independent physical verifications** that the stored
analytical results are consistent with each other and with
first-principles calculations.

The verification strategy has three tiers, each with a distinct
purpose and level of independence.

---

## Tier 1: Standalone Tests

**Location**: `test/standalone/`

**Principle**: Verify properties of QAtlas source values using
*only* QAtlas itself -- no Lattice2D, no ED, no external
computation. These tests check internal consistency and mathematical
identities.

**Examples**:

| Test file | What is checked |
|-----------|-----------------|
| `test_bethe_ansatz.jl` | $e_0 = J(1/4 - \ln 2)$ numerical value, $J$ scaling, consistency with finite-$N$ spectra |
| `test_ising_onsager_yang.jl` | $T_c = 2J/\ln(1+\sqrt{2})$, $\sinh(2\beta_c J) = 1$, $M(T_c) = 0$, $M(\beta \to \infty) = 1$ |
| `test_ising_square_pfaffian.jl` | $Z(\beta=0) = 2^N$, $L_x \leftrightarrow L_y$ symmetry |
| `test_universality_exponents.jl` | Scaling relations: $\alpha + 2\beta + \gamma = 2$, $\gamma = \nu(2 - \eta)$, hyperscaling |

**What this tier cannot do**: It cannot detect systematic errors
where the source formula itself is wrong, because there is no
independent computation to compare against.

---

## Tier 2: Verification Tests

**Location**: `test/verification/`

**Principle**: Cross-check QAtlas source values against an
**independent computation** from a different code path. Typically
this means comparing hardcoded analytical formulas
(`src/models/`) against numerical diagonalization of real-space
Hamiltonians built from `Lattice2D`.

The two paths must agree to numerical precision:

- **Path A** (src): Analytical formula, Bloch dispersion, or
  Bethe ansatz result -- implemented once, tested many times.
- **Path B** (test): Real-space Hamiltonian construction via
  `Lattice2D.bonds()` + exact diagonalization -- completely
  generic, no model-specific knowledge.

**Examples**:

| Test file | Path A (src) | Path B (test) |
|-----------|-------------|---------------|
| `test_heisenberg_dimer.jl` | Hardcoded $\{-3J/4, J/4, J/4, J/4\}$ | `build_spinhalf_heisenberg` + `eigvals` |
| `test_graphene_tight_binding.jl` | Closed-form Bloch $E_\pm = \pm t\sqrt{\ldots}$ | `build_tight_binding` + `eigvals` |
| `test_bloch_generic.jl` | Hardcoded Bloch formulas | Generic `bloch_tb_spectrum` builder |
| `test_tfim_gap_closure.jl` | BdG quasiparticle spectrum | Full $2^N$ ED via `build_tfim` |

**What this tier cannot do**: Both paths could share a systematic
error in the underlying physics (e.g., a wrong sign convention).
Tier 3 addresses this.

---

## Tier 3: Cross-Verification (Universality ↔ Models)

**Location**: `test/verification/test_universality_cross_check.jl`

**Principle**: Connect results from two **independently sourced
theoretical lines** -- typically a universality-class prediction
(sourced from CFT / Coulomb gas) and a model-specific exact result
(sourced from Onsager / Yang / Bethe). If both agree, the QAtlas
data is validated from two independent origins in the physics
literature.

This is the tier that makes QAtlas's data **physically verified**,
not just computationally consistent.

See [Cross-Verification: Universality ↔ Models](cross-checks.md)
for the complete table of all cross-checks.

---

## Design Principles

1. **Every `src/` value must have at least one Tier 2 test**.
   No analytical formula enters the source code without being
   verified against an independent numerical calculation.

2. **Cross-verification closes the loop**. Tier 3 tests connect
   different theoretical frameworks (e.g., Bethe ansatz vs ED,
   Yang's formula vs CFT exponents), ensuring that QAtlas's data
   is not just internally consistent but physically correct.

3. **Multiple computation paths are a feature, not redundancy**.
   For tight-binding models, three independent paths exist:
   hardcoded Bloch formula, generic Bloch builder, and real-space
   ED. For the TFIM, two paths exist: BdG quasiparticle spectrum
   and full $2^N$ ED. Agreement across all paths provides strong
   evidence of correctness.

4. **Test sizes are kept small enough for exact comparison**.
   All verification tests use system sizes where ED is feasible
   ($N \leq 16$ for spin models, $L_x, L_y \leq 6$ for
   tight-binding), enabling comparison to machine precision
   rather than statistical tolerance.

---

## Further Reading

- [Cross-Verification Table](cross-checks.md) -- all 8 universality
  ↔ model cross-checks
- [Entanglement Verification](entanglement.md) -- central charge
  extraction from $S(l)$
- [Disordered Systems](disordered.md) -- IRFP and random singlet
  verification

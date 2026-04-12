# O(n) Models: XY and Heisenberg

## Overview

The O($n$) model describes a classical spin $\mathbf{S}_i \in \mathbb{R}^n$
with $|\mathbf{S}_i| = 1$ on each lattice site, coupled by the
Hamiltonian

$$H = -J \sum_{\langle i,j \rangle} \mathbf{S}_i \cdot \mathbf{S}_j$$

The global symmetry is the orthogonal group O($n$). Special cases:

| $n$ | Name | Symmetry | Physical realisations |
|-----|------|----------|----------------------|
| 1 | [Ising](ising.md) | $\mathbb{Z}_2$ | Uniaxial magnets, lattice gas |
| 2 | XY | O(2) $\cong$ U(1) | Superfluid $^4$He, thin-film magnets |
| 3 | Heisenberg | O(3) $\cong$ SU(2) | Isotropic magnets (EuO, EuS) |

---

## XY Model ($n = 2$)

### $d = 2$ --- Berezinskii-Kosterlitz-Thouless (BKT) Transition

The 2D XY model has **no spontaneous symmetry breaking** at any
$T > 0$ (Mermin-Wagner theorem). Instead, it undergoes a BKT
transition at $T_{\mathrm{BKT}}$: below $T_{\mathrm{BKT}}$
correlations decay algebraically (quasi-long-range order), while
above they decay exponentially.

The BKT transition is **not** described by standard power-law
critical exponents. Its hallmarks are:

| Property | Value / behaviour | Note |
|----------|------------------|------|
| $\eta(T_{\mathrm{BKT}})$ | $1/4$ | Universal jump at the transition |
| Correlation length above $T_{\mathrm{BKT}}$ | $\xi \sim \exp(b/\sqrt{T - T_{\mathrm{BKT}}})$ | Essential singularity, not power law |
| Superfluid stiffness | Universal jump $\rho_s(T_{\mathrm{BKT}}^-) = 2T_{\mathrm{BKT}}/\pi$ | Nelson-Kosterlitz criterion |
| Vortex-antivortex pairs | Bound for $T < T_{\mathrm{BKT}}$, unbound above | Topological defect mechanism |

!!! warning "No standard critical exponents"
    The BKT transition has $\nu = \infty$ in the usual sense
    (essential singularity). Standard exponents $\beta, \gamma,
    \delta$ are not defined for this transition because there is
    no true order parameter. QAtlas stores the universal value
    $\eta = 1/4$ at $T_c$ and flags the BKT nature.

### $d = 3$ --- Conformal Bootstrap

In $d = 3$, the XY model has a conventional second-order transition.
The most precise exponents come from the O(2) conformal bootstrap.

| Exponent | Value | Uncertainty | Reference |
|----------|-------|-------------|-----------|
| $\alpha$ | $-0.0146(8)$ | — | Chester, Landry, Liu, Poland, Simmons-Duffin, Su, Vichi (2020) |
| $\beta$ | $0.3485(2)$ | — | JHEP **2020**, 142 |
| $\gamma$ | $1.3177(5)$ | — | " |
| $\nu$ | $0.6717(1)$ | — | " |
| $\eta$ | $0.0381(2)$ | — | " |

---

## Heisenberg Model ($n = 3$)

### $d \leq 2$ --- Mermin-Wagner Theorem

The Mermin-Wagner theorem forbids spontaneous breaking of a
continuous symmetry in $d \leq 2$ at $T > 0$ for short-range
interactions. Consequently:

- **$d = 1$**: The 1D Heisenberg chain is disordered at all $T > 0$.
  At $T = 0$ the quantum spin-1/2 chain is critical with $c = 1$
  (Luttinger liquid), described by the SU(2)$_1$ WZW model.
- **$d = 2$**: The 2D classical Heisenberg model has no finite-$T$
  phase transition. (In contrast, the 2D XY model has the BKT
  transition, which does not break the continuous symmetry.)

### $d = 3$ --- Conformal Bootstrap

In $d = 3$, the Heisenberg model has a conventional second-order
transition. Best estimates from the O(3) conformal bootstrap:

| Exponent | Value | Uncertainty | Reference |
|----------|-------|-------------|-----------|
| $\alpha$ | $-0.1336(15)$ | — | Chester, Landry, Liu, Poland, Simmons-Duffin, Su, Vichi (2020) |
| $\beta$ | $0.3689(3)$ | — | JHEP **2020**, 142 |
| $\gamma$ | $1.3960(9)$ | — | " |
| $\nu$ | $0.7112(5)$ | — | " |
| $\eta$ | $0.0378(3)$ | — | " |

---

## $d \geq 4$ --- Mean-Field

The upper critical dimension for O($n$) models (all $n \geq 1$) is
$d_c = 4$. For $d \geq 4$ the exponents take [mean-field](mean-field.md)
values: $\beta = 1/2$, $\nu = 1/2$, $\gamma = 1$, $\eta = 0$.

---

## QAtlas API

```julia
using QAtlas

# XY d = 3: numerical (Float64 + _err)
e_xy = QAtlas.fetch(Universality(:XY), CriticalExponents(); d=3)
# (β = 0.3485, β_err = 2e-4, ν = 0.6717, ν_err = 1e-4, ...)

# Heisenberg d = 3
e_heis = QAtlas.fetch(Universality(:Heisenberg), CriticalExponents(); d=3)
# (β = 0.3689, β_err = 3e-4, ν = 0.7112, ν_err = 5e-4, ...)

# d ≥ 4: mean-field
e_mf = QAtlas.fetch(Universality(:XY), CriticalExponents(); d=4)
# → same as fetch(MeanField(), CriticalExponents())
```

---

## References

- V. L. Berezinskii, Sov. Phys. JETP **32**, 493 (1971) --- BKT
  transition (part I).
- J. M. Kosterlitz, D. J. Thouless, J. Phys. C **6**, 1181 (1973)
  --- BKT transition.
- N. D. Mermin, H. Wagner, Phys. Rev. Lett. **17**, 1133 (1966) ---
  absence of long-range order in $d \leq 2$.
- S. M. Chester, W. Landry, J. Liu, D. Poland, D. Simmons-Duffin,
  N. Su, A. Vichi, "Carving out OPE space and precise O(2) model
  critical exponents", JHEP **2020**, 142 --- O(2) and O(3)
  bootstrap.
- D. R. Nelson, J. M. Kosterlitz, Phys. Rev. Lett. **39**, 1201
  (1977) --- universal superfluid-stiffness jump.

---

## Connections

- **Ising**: the $n = 1$ case; see [Ising](ising.md).
- **Heisenberg chain ($T = 0$)**: quantum $S = 1/2$ chain with
  $c = 1$; see [Heisenberg model](../models/quantum/heisenberg.md).
- **Mean-Field**: baseline for $d \geq 4$; see [mean-field](mean-field.md).

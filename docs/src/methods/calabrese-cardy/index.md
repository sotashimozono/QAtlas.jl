# Calabrese-Cardy Formula

## Overview

The Calabrese-Cardy formula gives the bipartite entanglement entropy of
a 1D critical system as a function of the subsystem size, the total
system size, and the **central charge** $c$ of the underlying conformal
field theory.

This formula is the primary tool for extracting the central charge from
a ground-state wavefunction, and is used throughout QAtlas's
[entanglement verification](../../verification/entanglement.md).

---

## Statement

### Periodic boundary conditions (PBC)

For a system of $N$ sites on a ring, with subsystem $A$ consisting of
$l$ consecutive sites:

$$S(l) = \frac{c}{3}\ln\!\left[\frac{N}{\pi a}\sin\!\left(\frac{\pi l}{N}\right)\right] + s_1$$

The prefactor $c/3$ arises because PBC creates **two** entanglement
cuts (between sites $l, l+1$ and sites $N, 1$).

### Open boundary conditions (OBC)

For a system of $N$ sites on an open chain, with subsystem $A$ being
sites $1, \ldots, l$:

$$S(l) = \frac{c}{6}\ln\!\left[\frac{2N}{\pi a}\sin\!\left(\frac{\pi l}{N}\right)\right] + s_1'$$

The prefactor is $c/6$ because OBC creates only **one** entanglement
cut (between sites $l$ and $l+1$).

!!! warning "OBC vs PBC prefactor"
    Using $c/3$ (PBC) for an OBC system halves the extracted central
    charge. This is a common mistake — always match the prefactor to
    the boundary conditions of the system being studied.

### Parameters

- $c$: central charge of the CFT (universal)
- $a$: UV cutoff (lattice spacing, typically $a = 1$)
- $s_1, s_1'$: non-universal additive constants

---

## Physical Origin

In a 1D critical system, the ground state is described by a conformal
field theory. The entanglement entropy of a subsystem of size $l$ in
a total system of size $N$ is determined by the two-point function of
**twist operators** in the orbifold CFT $\mathcal{C}^n / \mathbb{Z}_n$
(replica trick):

$$S = -\lim_{n \to 1} \frac{\partial}{\partial n} \mathrm{Tr}(\rho_A^n)$$

The $\sin(\pi l / N)$ factor is the **chord length** in the conformal
mapping from the infinite plane to the cylinder (PBC) or strip (OBC).

---

## Central Charge Extraction

To extract $c$ from numerical entropy data $\{S(l)\}$:

1. Compute $\xi(l) = \ln[(2N/\pi)\sin(\pi l/N)]$ (OBC conformal coordinate)
2. Perform linear regression: $S = (\text{slope}) \cdot \xi + \text{const}$
3. $c = 6 \times \text{slope}$ (OBC) or $c = 3 \times \text{slope}$ (PBC)

**Practical tips**:
- Skip boundary-adjacent points ($l = 1, 2$ and $l = N-1, N-2$) to
  reduce lattice artifacts.
- For the Heisenberg chain, use only **even** $l$ values to suppress
  the SU(2) alternating correction $(-1)^l f(l)$.

---

## Known Central Charges

| System | $c$ | QAtlas verification |
|--------|-----|---------------------|
| [TFIM](../../models/quantum/tfim.md) at $h = J$ | $1/2$ | [→](../../verification/entanglement.md) |
| [Heisenberg chain](../../models/quantum/heisenberg.md) | $1$ | [→](../../verification/entanglement.md) |
| Free boson (Luttinger liquid) | $1$ | — |
| Free Dirac fermion | $1/2$ | — |

---

## References

- P. Calabrese, J. Cardy, "Entanglement entropy and quantum field
  theory", J. Stat. Mech. **0406**, P06002 (2004) — original derivation.
  Eq. (7) for PBC, Eq. (19) for OBC.
- P. Calabrese, J. Cardy, "Entanglement entropy and conformal field
  theory", J. Phys. A **42**, 504005 (2009) — review with extensions.
- C. Holzhey, F. Larsen, F. Wilczek, "Geometric and renormalized
  entropy in conformal field theory", Nucl. Phys. B **424**, 443 (1994)
  — earlier result for the $c/3$ formula.

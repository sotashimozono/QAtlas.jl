# Entanglement Entropy and Central Charge

## Overview

The entanglement entropy of a 1D critical ground state encodes the
**central charge** $c$ of the underlying conformal field theory via
the [Calabrese-Cardy formula](../methods/calabrese-cardy/index.md).
QAtlas extracts $c$ from exact-diagonalization ground states and
compares against the known CFT values, providing a verification that
connects quantum information (entanglement) to conformal field theory
(central charge).

---

## Method

### Calabrese-Cardy Formula (OBC)

For an open chain of $N$ sites with subsystem $A = \{1, \ldots, l\}$:

$$S(l) = \frac{c}{6}\ln\!\left[\frac{2N}{\pi}\sin\!\left(\frac{\pi l}{N}\right)\right] + s_1'$$

The prefactor is $c/6$ (not $c/3$) because OBC creates only **one**
entanglement cut. See
[calc/calabrese-cardy-obc-vs-pbc.md](../calc/calabrese-cardy-obc-vs-pbc.md)
for the derivation and physical reason for the factor-of-2 difference.

### Extraction Procedure

1. Compute the ground-state wavefunction $|\psi_0\rangle$ via ED.
2. For each bipartition $l = 1, \ldots, N-1$, compute the reduced
   density matrix $\rho_A = \mathrm{tr}_B |\psi_0\rangle\langle\psi_0|$
   and the von Neumann entropy $S(l) = -\mathrm{tr}(\rho_A \ln \rho_A)$.
3. Define the conformal coordinate
   $\xi(l) = \ln[(2N/\pi)\sin(\pi l/N)]$.
4. Perform linear regression: $S = \text{slope} \cdot \xi + \text{const}$.
5. Extract $c = 6 \times \text{slope}$.

---

## Results

### TFIM at $h = J$ (Ising CFT, $c = 1/2$)

| Parameter | Value |
|-----------|-------|
| System size | $N = 14$, OBC |
| Boundary conditions | Open |
| Expected $c$ | $1/2$ |
| Extracted $c$ | $\approx 0.5$ (within 10%) |
| Dominant error source | Finite-size corrections at small $N$ |

At the TFIM critical point $h = J$, the system is in the Ising CFT
universality class with $c = 1/2$. The entanglement entropy $S(l)$
shows the characteristic logarithmic growth, and the extracted
central charge agrees with $1/2$ to approximately 10% at $N = 14$.

### Heisenberg Chain (SU(2)$_1$ WZW, $c = 1$)

| Parameter | Value |
|-----------|-------|
| System size | $N = 12$, OBC |
| Boundary conditions | Open |
| Expected $c$ | $1$ |
| Extracted $c$ | $\approx 1.0$ (within 20%, even-$l$ fit) |
| Dominant error source | SU(2) alternating correction $(-1)^l f(l)$ |

The spin-1/2 Heisenberg chain is a $c = 1$ Luttinger liquid
(SU(2)$_1$ WZW model). A crucial practical detail: the SU(2)
symmetry produces an **alternating correction** $(-1)^l$ to $S(l)$
that contaminates the fit if all $l$ values are used. Restricting
the fit to **even $l$ only** suppresses this correction and yields
$c \approx 1$ within 20%.

### Gapped Phases: Area Law

Away from criticality, the entanglement entropy saturates to a
constant as $l$ increases (area law). QAtlas verifies this for:

- **TFIM $h \gg J$**: paramagnetic phase, product-state ground state,
  $S \to 0$.
- **TFIM $h \ll J$**: ferromagnetic phase, $S$ saturates to a small
  constant (boundary contribution).

The absence of logarithmic growth in the gapped phase confirms that
the Calabrese-Cardy formula applies specifically at criticality.

---

## Limitations

### Peschel Method and Kramers-Wannier Duality

For the TFIM, an alternative route to entanglement entropy uses the
**Peschel method**: compute the correlation matrix of Jordan-Wigner
fermions and extract the entropy from its eigenvalues. This avoids
full ED and scales as $O(N^3)$.

However, the Peschel method requires the ground state to be
expressible as a free-fermion Slater determinant. The Kramers-Wannier
(KW) duality maps the TFIM to a dual model where the boundary
conditions may differ. Under OBC, the KW duality introduces a
subtlety at the boundary that can shift the effective system size
by one site. QAtlas uses the full ED approach to avoid this issue.

### Finite-Size Accuracy

At the system sizes accessible to ED ($N \leq 16$), the extracted
$c$ has $\sim 10$--$20\%$ error. This is expected: the Calabrese-Cardy
formula is asymptotically exact as $N \to \infty$, and subleading
corrections are significant at small $N$. The test asserts
**qualitative** agreement (correct order of magnitude and sign)
rather than high precision.

---

## QAtlas Test

```julia
# test/verification/test_entanglement_central_charge.jl

# TFIM: c ≈ 1/2
c_tfim = extract_central_charge(tfim_ground_state, N=14, bc=:OBC)
@test abs(c_tfim - 0.5) / 0.5 < 0.10

# Heisenberg: c ≈ 1 (even-l fit)
c_heis = extract_central_charge(heis_ground_state, N=12, bc=:OBC, even_only=true)
@test abs(c_heis - 1.0) / 1.0 < 0.20
```

---

## Connections

- **Calabrese-Cardy formula**: [method page](../methods/calabrese-cardy/index.md),
  [OBC vs PBC derivation](../calc/calabrese-cardy-obc-vs-pbc.md)
- **Ising universality**: $c = 1/2$ is the central charge of the
  minimal model $\mathcal{M}(3,4)$; see [Ising](../universalities/ising.md)
- **Cross-check #3**: TFIM $c = 1/2$ appears in the
  [cross-verification table](cross-checks.md)
- **Disordered systems**: entanglement in random chains;
  see [disordered](disordered.md)

---

## References

- P. Calabrese, J. Cardy, J. Stat. Mech. **0406**, P06002 (2004) ---
  original Calabrese-Cardy formula.
- I. Peschel, J. Phys. A **36**, L205 (2003) --- free-fermion
  entanglement from correlation matrix.
- G. Vidal, J. I. Latorre, E. Rico, A. Kitaev, Phys. Rev. Lett.
  **90**, 227902 (2003) --- entanglement entropy in spin chains.

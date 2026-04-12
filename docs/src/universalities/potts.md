# Potts Universality Classes

## Overview

The $q$-state Potts model generalises the Ising model ($q = 2$) to
a spin variable $\sigma_i \in \{1, 2, \ldots, q\}$ with Hamiltonian

$$H = -J \sum_{\langle i,j \rangle} \delta_{\sigma_i, \sigma_j}$$

where $\delta$ is the Kronecker delta. The global symmetry group is
$S_q$ (permutation of the $q$ colours).

In $d = 2$ the phase transition is second-order for $q \leq 4$ and
first-order for $q > 4$. For $q = 2, 3, 4$ the critical exponents
are known **exactly** via the Coulomb gas mapping.

---

## $q = 2$ --- Ising

The $q = 2$ Potts model is equivalent to the Ising model with a
rescaled coupling. See [Ising universality class](ising.md).

---

## $q = 3$ --- Three-State Potts ($d = 2$)

**Symmetry**: $S_3$ (permutation of 3 colours).

| Exponent | Value | Physical meaning | Reference |
|----------|-------|-----------------|-----------|
| $\alpha$ | $1/3$ | Specific heat: $C \sim \lvert T - T_c\rvert^{-1/3}$ | Baxter (1973); den Nijs (1979) |
| $\beta$ | $1/9$ | Order parameter | Alexander (1975); Nienhuis (1984) |
| $\gamma$ | $13/9$ | Susceptibility | Scaling relation |
| $\nu$ | $5/6$ | Correlation length | Coulomb gas |
| $\eta$ | $4/15$ | Anomalous dimension | Coulomb gas |
| $\delta$ | $14$ | Critical isotherm | $\delta = 1 + \gamma/\beta$ |
| $c$ | $4/5$ | Central charge (CFT) | Minimal model $\mathcal{M}(5,6)$ |

### Scaling Relations

$$\alpha + 2\beta + \gamma = \tfrac{1}{3} + \tfrac{2}{9} + \tfrac{13}{9} = 2 \quad\checkmark$$
$$\gamma = \nu(2 - \eta) = \tfrac{5}{6}\bigl(2 - \tfrac{4}{15}\bigr) = \tfrac{13}{9} \quad\checkmark$$
$$2 - \alpha = d\nu \implies \tfrac{5}{3} = 2 \cdot \tfrac{5}{6} \quad\checkmark$$

---

## $q = 4$ --- Four-State Potts ($d = 2$)

**Symmetry**: $S_4$.

The $q = 4$ case is **marginal**: the transition is second-order but
sits at the boundary of the first-order regime. Logarithmic
corrections to scaling appear, making numerical extraction of
exponents notoriously difficult.

| Exponent | Value | Note | Reference |
|----------|-------|------|-----------|
| $\alpha$ | $2/3$ | Strong specific-heat divergence | Baxter (1973) |
| $\beta$ | $1/12$ | Very small — order parameter onset is slow | den Nijs (1979) |
| $\gamma$ | $7/6$ | Susceptibility | Scaling relation |
| $\nu$ | $2/3$ | Correlation length | Coulomb gas |
| $\eta$ | $1/4$ | Same anomalous dimension as Ising $d = 2$ | — |
| $\delta$ | $15$ | Same as Ising $d = 2$ | $\delta = 1 + \gamma/\beta$ |
| $c$ | $1$ | Central charge (CFT) | Orbifold $c = 1$ theory |

### Logarithmic Corrections

At $q = 4$ the RG $\beta$-function has a **double zero** at the
fixed point, leading to multiplicative logarithmic corrections of
the form

$$M \sim (T_c - T)^{1/12} \lvert\ln(T_c - T)\rvert^{1/8}$$

These corrections make finite-size scaling analysis substantially
harder than for $q = 2$ or $q = 3$.

### Scaling Relations

$$\alpha + 2\beta + \gamma = \tfrac{2}{3} + \tfrac{1}{6} + \tfrac{7}{6} = 2 \quad\checkmark$$
$$2 - \alpha = d\nu \implies \tfrac{4}{3} = 2 \cdot \tfrac{2}{3} \quad\checkmark$$

---

## QAtlas API

```julia
using QAtlas

# q = 3, d = 2: exact (Rational{Int})
e3 = QAtlas.fetch(Universality(:Potts3), CriticalExponents(); d=2)
# (β = 1//9, ν = 5//6, γ = 13//9, η = 4//15, ...)

# q = 4, d = 2: exact (Rational{Int})
e4 = QAtlas.fetch(Universality(:Potts4), CriticalExponents(); d=2)
# (β = 1//12, ν = 2//3, γ = 7//6, η = 1//4, ...)
```

---

## References

- R. J. Baxter, "Potts model at the critical temperature",
  J. Phys. C **6**, L445 (1973) --- exact critical temperature.
- R. J. Baxter, *Exactly Solved Models in Statistical Mechanics*
  (Academic Press, 1982), Ch. 12.
- B. Nienhuis, "Critical behavior of two-dimensional spin models and
  charge asymmetry in the Coulomb gas", J. Stat. Phys. **34**, 731
  (1984) --- Coulomb gas derivation of exact exponents.
- M. P. M. den Nijs, J. Phys. A **12**, 1857 (1979) --- exponent
  relations via Coulomb gas.
- J. Salas, A. D. Sokal, "Logarithmic corrections and finite-size
  scaling in the two-dimensional 4-state Potts model",
  J. Stat. Phys. **88**, 567 (1997) --- log corrections.

---

## Connections

- **Ising**: the $q = 2$ case is the [Ising universality class](ising.md).
- **Percolation**: the $q \to 1$ limit gives [percolation](percolation.md)
  via the Kasteleyn-Fortuin mapping.
- **Scaling relations**: verified using the same algebraic framework
  as [calc/ising-scaling-relations.md](../calc/ising-scaling-relations.md).

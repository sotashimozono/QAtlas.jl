# Percolation Universality Class

## Overview

Percolation describes a **geometric** phase transition: sites (or
bonds) on a lattice are occupied independently with probability $p$.
At $p = p_c$ an infinite connected cluster first appears. Near $p_c$
the cluster statistics exhibit power-law behaviour governed by
critical exponents that depend only on the spatial dimension $d$.

Unlike thermal phase transitions, percolation has no Hamiltonian and
no temperature. The control parameter is the occupation probability
$p$ and the "order parameter" is the probability $P_\infty$ that a
randomly chosen site belongs to the infinite cluster.

**Symmetry**: Permutation symmetry of the $q \to 1$ limit of the
$q$-state Potts model (Kasteleyn-Fortuin mapping).

---

## $d = 2$ --- Exact Critical Exponents

In $d = 2$, the percolation exponents are known exactly via the
Coulomb gas mapping and SLE (Schramm-Loewner evolution).

| Exponent | Value | Physical meaning | Reference |
|----------|-------|-----------------|-----------|
| $\beta$ | $5/36$ | Infinite-cluster density: $P_\infty \sim (p - p_c)^\beta$ | Nienhuis (1982) J. Stat. Phys. **34**, 731 |
| $\nu$ | $4/3$ | Correlation length: $\xi \sim \lvert p - p_c\rvert^{-\nu}$ | den Nijs (1979); Nienhuis (1982) |
| $\gamma$ | $43/18$ | Mean cluster size: $\langle s \rangle \sim \lvert p - p_c\rvert^{-\gamma}$ | Scaling relation $\gamma = \nu(2 - \eta)$ |
| $\eta$ | $5/24$ | Pair connectivity: $g(r) \sim r^{-(d-2+\eta)}$ at $p_c$ | Coulomb gas |
| $\alpha$ | $-2/3$ | Cluster-number density | $\alpha = 2 - d\nu = 2 - 8/3$ |
| $\delta$ | $91/5$ | $P_\infty \sim h^{1/\delta}$ at $p_c$ with ghost field | Scaling relation |
| $\sigma$ | $36/91$ | Cluster-size distribution: $n_s \sim s^{-\tau}$ with $\tau = 1 + 1/(\sigma\beta\delta)$ | — |
| $\tau$ | $187/91$ | Fisher exponent: $n_s \sim s^{-\tau} f(s/s_\xi)$ | — |

### Scaling Relations

The exponents satisfy all standard scaling relations, verified
algebraically in QAtlas (see [calculation note](../calc/ising-scaling-relations.md)
for the general framework):

$$\alpha + 2\beta + \gamma = -\tfrac{2}{3} + \tfrac{5}{18} + \tfrac{43}{18} = 2 \quad\checkmark$$
$$\gamma = \nu(2 - \eta) = \tfrac{4}{3}\left(2 - \tfrac{5}{24}\right) = \tfrac{43}{18} \quad\checkmark$$

---

## $d = 3$ --- Numerical

In $d = 3$, no exact solution is known. The best estimates come from
large-scale Monte Carlo simulations.

| Exponent | Value | Reference |
|----------|-------|-----------|
| $\beta$ | $0.4181(8)$ | Wang, Zhou, et al. (2013) Phys. Rev. E **87**, 052107 |
| $\nu$ | $0.8765(12)$ | " |
| $\gamma$ | $1.793(3)$ | " |
| $\eta$ | $-0.046(8)$ | " |

---

## $d \geq 6$ --- Mean-Field

The upper critical dimension for percolation is $d_c = 6$. For
$d \geq 6$ the exponents take [mean-field](mean-field.md) values on
a Bethe lattice:

$$\beta = 1, \quad \nu = 1/2, \quad \gamma = 1, \quad \eta = 0, \quad \alpha = -1$$

Note that $\alpha < 0$ (no specific-heat divergence) and
$\beta = 1$ (linear onset of $P_\infty$) differ from the thermal
mean-field values, reflecting the geometric nature of the transition.

---

## QAtlas API

```julia
using QAtlas

# d = 2: exact (Rational{Int})
e = QAtlas.fetch(Universality(:Percolation), CriticalExponents(); d=2)
# (β = 5//36, ν = 4//3, γ = 43//18, η = 5//24, ...)

# d = 3: numerical (Float64 + _err)
e3 = QAtlas.fetch(Universality(:Percolation), CriticalExponents(); d=3)

# d ≥ 6: mean-field
e6 = QAtlas.fetch(Universality(:Percolation), CriticalExponents(); d=6)
```

---

## References

- B. Nienhuis, "Exact critical point and critical exponents of O(n)
  models in two dimensions", Phys. Rev. Lett. **49**, 1062 (1982).
- M. P. M. den Nijs, "A relation between the temperature exponents
  of the eight-vertex and q-state Potts model", J. Phys. A **12**,
  1857 (1979).
- D. Stauffer, A. Aharony, *Introduction to Percolation Theory*
  (Taylor & Francis, 1994).
- G. Grimmett, *Percolation* (Springer, 1999).

---

## Connections

- **Potts $q \to 1$**: The Kasteleyn-Fortuin mapping relates bond
  percolation to the $q \to 1$ limit of the [Potts model](potts.md).
- **Scaling relations**: verified in [calc/ising-scaling-relations.md](../calc/ising-scaling-relations.md)
  (same algebraic framework, different exponent values).

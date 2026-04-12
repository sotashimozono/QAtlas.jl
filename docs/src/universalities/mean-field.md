# Mean-Field (Landau)

## Overview

The Landau mean-field theory provides the **baseline reference** for
critical exponents. It describes the behaviour of any system above
its upper critical dimension $d_c$, where fluctuations are
irrelevant and the saddle-point approximation to the partition
function becomes exact.

The mean-field exponents are independent of spatial dimension (for
$d \geq d_c$), symmetry group, and microscopic details. They
represent the simplest possible universality class and serve as the
starting point for $\epsilon$-expansion calculations
($d = d_c - \epsilon$).

---

## Critical Exponents

| Exponent | Value | Physical meaning |
|----------|-------|-----------------|
| $\alpha$ | $0$ | Specific heat: finite jump (discontinuity), no divergence |
| $\beta$ | $1/2$ | Order parameter: $M \sim (T_c - T)^{1/2}$ |
| $\gamma$ | $1$ | Susceptibility: $\chi \sim \lvert T - T_c\rvert^{-1}$ |
| $\delta$ | $3$ | Critical isotherm: $M \sim h^{1/3}$ |
| $\nu$ | $1/2$ | Correlation length: $\xi \sim \lvert T - T_c\rvert^{-1/2}$ |
| $\eta$ | $0$ | No anomalous dimension: $G(r) \sim r^{-(d-2)}$ (Ornstein-Zernike) |

### Derivation from Landau Free Energy

The Landau free energy functional is

$$F[m] = \int d^d x \left[\frac{r}{2}m^2 + \frac{u}{4}m^4 + \frac{c}{2}(\nabla m)^2 - hm\right]$$

where $r \propto (T - T_c)$ changes sign at the transition. Minimising
$F$ with respect to $m$ (saddle-point):

- $h = 0$, $r < 0$: $m = \pm\sqrt{-r/u} \sim (T_c - T)^{1/2}$ gives $\beta = 1/2$.
- $r = 0$: $um^3 = h$ gives $\delta = 3$.
- $h = 0$: $\chi = \partial m/\partial h = 1/|r| \sim |T - T_c|^{-1}$ gives $\gamma = 1$.
- The Gaussian propagator $G(k) = 1/(r + ck^2)$ gives $\xi^2 = c/|r|$, so $\nu = 1/2$ and $\eta = 0$.

### Scaling Relations

All four standard scaling relations are satisfied exactly:

$$\alpha + 2\beta + \gamma = 0 + 1 + 1 = 2 \quad\checkmark$$
$$\gamma = \beta(\delta - 1) = \tfrac{1}{2}\cdot 2 = 1 \quad\checkmark$$
$$\gamma = \nu(2 - \eta) = \tfrac{1}{2}\cdot 2 = 1 \quad\checkmark$$

The **Josephson (hyperscaling) relation** $2 - \alpha = d\nu$ gives
$d = 4$. This relation is violated for $d > 4$ because the Gaussian
fixed point acquires dangerous irrelevant variables. The mean-field
exponents remain valid for $d > 4$, but hyperscaling does not hold.

---

## Upper Critical Dimensions

| Model | Symmetry | $d_c$ | Reference |
|-------|----------|-------|-----------|
| [Ising](ising.md) | $\mathbb{Z}_2$ | 4 | Ginzburg criterion |
| [XY / Heisenberg](on-models.md) | O($n$), $n \geq 2$ | 4 | Wilson-Fisher (1972) |
| [Percolation](percolation.md) | $S_q$, $q \to 1$ | 6 | Toulouse (1974) |
| Directed percolation | — | $4 + 1$ | Janssen (1981) |

---

## QAtlas API

```julia
using QAtlas

# Mean-field exponents (no dimension needed)
e = QAtlas.fetch(MeanField(), CriticalExponents())
# (β = 1//2, ν = 1//2, γ = 1//1, η = 0//1, δ = 3//1, α = 0//1)

# Equivalently, any universality class at d ≥ d_c returns these
e_ising4 = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=4)
# → same values as MeanField()
```

All mean-field exponents are stored as `Rational{Int}`, enabling
exact algebraic verification of scaling relations.

---

## References

- L. D. Landau, "On the theory of phase transitions",
  Zh. Eksp. Teor. Fiz. **7**, 19 (1937) --- original Landau theory.
- V. L. Ginzburg, "Some remarks on phase transitions of the second
  kind and the microscopic theory of ferroelectric materials",
  Sov. Phys. Solid State **2**, 1824 (1960) --- Ginzburg criterion
  for the validity of mean-field theory.
- K. G. Wilson, M. E. Fisher, Phys. Rev. Lett. **28**, 240 (1972) ---
  $\epsilon$-expansion showing mean-field is exact for $d \geq 4$.

---

## Connections

- Every universality class in QAtlas reduces to mean-field at or
  above its upper critical dimension.
- The $\epsilon$-expansion computes corrections to these exponents
  perturbatively in $\epsilon = d_c - d$.
- **Scaling relations**: verified in [calc/ising-scaling-relations.md](../calc/ising-scaling-relations.md).

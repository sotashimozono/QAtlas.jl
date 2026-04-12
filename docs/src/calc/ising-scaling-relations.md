# Ising Critical Exponents: Scaling Relations

## Setup

Near a second-order phase transition at $T = T_c$, thermodynamic
quantities diverge (or vanish) as power laws characterised by
critical exponents. For a system with reduced temperature
$t = (T - T_c)/T_c$ and external field $h$:

| Exponent | Definition | Quantity |
|---|---|---|
| $\alpha$ | $C \sim |t|^{-\alpha}$ | Specific heat |
| $\beta$ | $M \sim (-t)^{\beta}$ ($t < 0$, $h = 0$) | Order parameter |
| $\gamma$ | $\chi \sim |t|^{-\gamma}$ | Susceptibility |
| $\delta$ | $M \sim h^{1/\delta}$ ($t = 0$) | Critical isotherm |
| $\nu$ | $\xi \sim |t|^{-\nu}$ | Correlation length |
| $\eta$ | $G(r) \sim r^{-(d-2+\eta)}$ ($t = 0$) | Anomalous dimension |

For the 2D Ising universality class, the exact values are:

$$\alpha = 0,\quad \beta = \tfrac{1}{8},\quad \gamma = \tfrac{7}{4},\quad
  \delta = 15,\quad \nu = 1,\quad \eta = \tfrac{1}{4}$$

## Calculation

### Rushbrooke inequality (equality)

$$\alpha + 2\beta + \gamma = 2$$

**Physical content:** the specific heat ($\alpha$), order parameter
($\beta$), and susceptibility ($\gamma$) are related through the
free-energy scaling form $f(t, h) = |t|^{2-\alpha}\Phi(h/|t|^{\beta\delta})$.
Differentiating twice with respect to $h$ gives $\chi$, connecting
$\alpha$, $\beta$, and $\gamma$.

**2D Ising check:**

$$0 + 2 \cdot \tfrac{1}{8} + \tfrac{7}{4} = 0 + \tfrac{1}{4} + \tfrac{7}{4} = 2 \quad\checkmark$$

### Widom relation

$$\gamma = \beta(\delta - 1)$$

**Physical content:** the susceptibility divergence ($\gamma$) is
related to the order parameter ($\beta$) and critical isotherm
($\delta$) through the equation of state $h \sim M|M|^{\delta-1}$
near $T_c$.

**2D Ising check:**

$$\tfrac{1}{8}\cdot(15 - 1) = \tfrac{1}{8}\cdot 14 = \tfrac{7}{4} = \gamma \quad\checkmark$$

### Fisher relation

$$\gamma = \nu(2 - \eta)$$

**Physical content:** the susceptibility is the integral of the
correlation function, $\chi = \int G(r)\,d^d r$. Combining the
power-law decay $G(r) \sim r^{-(d-2+\eta)}$ with the correlation
length $\xi \sim |t|^{-\nu}$ as the cutoff gives this relation.

**2D Ising check:**

$$1 \cdot (2 - \tfrac{1}{4}) = \tfrac{7}{4} = \gamma \quad\checkmark$$

### Josephson (hyperscaling) relation

$$2 - \alpha = d\nu$$

**Physical content:** the singular part of the free energy density
scales as $f_s \sim \xi^{-d} \sim |t|^{d\nu}$, while $f_s \sim |t|^{2-\alpha}$
by definition. This relation involves the spatial dimension $d$
explicitly and is valid only below the upper critical dimension
$d_c = 4$.

**2D Ising check ($d = 2$):**

$$2 - 0 = 2 \cdot 1 = 2 \quad\checkmark$$

### Algebraic verification in QAtlas

In the QAtlas test suite, all four scaling relations are verified
using exact `Rational{Int}` arithmetic, avoiding any floating-point
round-off:

```julia
alpha, beta, gamma, delta, nu, eta = 0//1, 1//8, 7//4, 15//1, 1//1, 1//4

@assert alpha + 2*beta + gamma == 2       # Rushbrooke
@assert gamma == beta * (delta - 1)       # Widom
@assert gamma == nu * (2 - eta)           # Fisher
@assert 2 - alpha == 2 * nu               # Josephson (d=2)
```

All assertions pass with zero error.

## Result

$$\boxed{\begin{aligned}
  \alpha + 2\beta + \gamma &= 2 & &\text{(Rushbrooke)} \\
  \gamma &= \beta(\delta - 1) & &\text{(Widom)} \\
  \gamma &= \nu(2 - \eta) & &\text{(Fisher)} \\
  2 - \alpha &= d\nu & &\text{(Josephson)}
\end{aligned}}$$

All four relations are satisfied exactly by the 2D Ising exponents
$\{\alpha = 0,\, \beta = 1/8,\, \gamma = 7/4,\, \delta = 15,\,
\nu = 1,\, \eta = 1/4\}$.

## References

- L. P. Kadanoff, W. Gotze, D. Hamblen, R. Hecht, E. A. S. Lewis, V. V. Palciauskas, M. Rayl, J. Swift, D. Aspnes, J. Kane, Rev. Mod. Phys. **39**, 395 (1967) — scaling relations review.
- B. Widom, J. Chem. Phys. **43**, 3898 (1965) — Widom scaling.
- M. E. Fisher, Phys. Rev. **180**, 594 (1969) — Fisher relation.
- B. D. Josephson, Proc. Phys. Soc. **92**, 269 (1967) — hyperscaling.

## Used by

- [Ising Universality Class](../universalities/ising.md)
- `test/standalone/test_universality_exponents.jl`

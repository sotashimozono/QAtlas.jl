# KPZ Universality Class

## Overview

The Kardar-Parisi-Zhang (KPZ) universality class describes
**non-equilibrium stochastic growth** of interfaces. It is
fundamentally different from equilibrium critical phenomena: there is
no partition function, no free energy, and the relevant exponents
characterise dynamic scaling of a growing surface rather than
static thermodynamic singularities.

The KPZ equation in $d$ spatial dimensions is

$$\frac{\partial h}{\partial t} = \nu_0 \nabla^2 h + \frac{\lambda}{2}(\nabla h)^2 + \eta(\mathbf{x}, t)$$

where $h(\mathbf{x}, t)$ is the interface height, $\nu_0$ is a
smoothing coefficient, $\lambda$ is the non-linear growth coupling,
and $\eta$ is Gaussian white noise with
$\langle\eta(\mathbf{x},t)\eta(\mathbf{x}',t')\rangle = 2D\,\delta^d(\mathbf{x}-\mathbf{x}')\delta(t-t')$.

**Systems in this class**: ballistic deposition, Eden growth,
polynuclear growth, directed polymers in random media, TASEP
(totally asymmetric simple exclusion process).

---

## $1+1$ Dimensions --- Exact Exponents

In $1+1$D (one spatial + one temporal dimension), the KPZ exponents
are known exactly.

### Growth Exponents

| Exponent | Value | Definition | Reference |
|----------|-------|-----------|-----------|
| $\beta$ | $1/3$ | Growth exponent: $W(t) \sim t^\beta$ at early times | KPZ (1986) |
| $\alpha$ | $1/2$ | Roughness exponent: $W_{\mathrm{sat}} \sim L^\alpha$ | KPZ (1986) |
| $z$ | $3/2$ | Dynamic exponent: $t_{\times} \sim L^z$ | KPZ (1986) |

Here $W(t) = \sqrt{\langle(h - \langle h\rangle)^2\rangle}$ is the
interface width (roughness).

!!! note "Not equilibrium critical exponents"
    These are **not** the standard $\alpha, \beta$ of thermal
    phase transitions. The KPZ $\beta$ is the growth exponent
    (width vs time), and $\alpha$ is the roughness exponent
    (saturation width vs system size). Do not confuse with
    order-parameter or specific-heat exponents.

### Galilean Invariance Constraint

The non-linear term $(\nabla h)^2$ endows the KPZ equation with
**Galilean invariance** under tilted-frame transformations. This
symmetry enforces the exact relation

$$\alpha + z = 2$$

which, combined with the scaling relation $z = \alpha / \beta$,
fixes all three exponents from a single one:

$$\alpha = 1/2, \quad z = 3/2, \quad \beta = \alpha/z = 1/3$$

### Exact Distribution

Beyond the exponents, the full probability distribution of the
height fluctuations is known exactly in $1+1$D:

- **Flat initial condition**: $\chi \sim t^{1/3} \xi_{\mathrm{GOE}}$
  (Tracy-Widom GOE distribution)
- **Curved initial condition**: $\chi \sim t^{1/3} \xi_{\mathrm{GUE}}$
  (Tracy-Widom GUE distribution)

This was proven rigorously via the connection to the TASEP and
random matrix theory (Sasamoto-Spohn 2010, Amir-Corwin-Quastel
2011).

---

## Higher Dimensions

For $d \geq 2$ spatial dimensions, no exact solution is known. The
upper critical dimension of KPZ (if it exists) remains an open
problem. Numerical estimates for $2+1$D give $\alpha \approx 0.393$,
$\beta \approx 0.240$, $z \approx 1.607$.

---

## QAtlas API

```julia
using QAtlas

# 1+1D KPZ: exact growth exponents
g = QAtlas.fetch(Universality(:KPZ), GrowthExponents(); d=1)
# (β = 1//3, α = 1//2, z = 3//2)
```

Note the use of `GrowthExponents()` rather than `CriticalExponents()`
to reflect the non-equilibrium nature of the KPZ class.

---

## References

- M. Kardar, G. Parisi, Y.-C. Zhang, "Dynamic scaling of growing
  interfaces", Phys. Rev. Lett. **56**, 889 (1986) --- original KPZ
  equation and exponent prediction.
- T. Sasamoto, H. Spohn, "One-dimensional Kardar-Parisi-Zhang
  equation: an exact solution and its universality",
  Phys. Rev. Lett. **104**, 230602 (2010) --- exact height
  distribution.
- G. Amir, I. Corwin, J. Quastel, "Probability distribution of the
  free energy of the continuum directed random polymer in 1 + 1
  dimensions", Comm. Pure Appl. Math. **64**, 466 (2011).
- I. Corwin, "The Kardar-Parisi-Zhang equation and universality
  class", Random Matrices Theory Appl. **1**, 1130001 (2012) ---
  review.

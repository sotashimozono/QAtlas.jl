# Ising Universality Class

## Overview

The Ising universality class describes second-order phase transitions
with $\mathbb{Z}_2$ symmetry breaking. It is the most fundamental
universality class in statistical physics, exactly solved in $d = 2$
and determined to extraordinary precision in $d = 3$ via the conformal
bootstrap.

**Symmetry**: $\mathbb{Z}_2$ (spin-flip $\sigma \to -\sigma$).

**Models in this class**: 2D classical Ising, 1+1D TFIM at $h = J$,
liquid-gas critical point, binary alloy order-disorder, uniaxial
ferromagnets.

**CFT ($d = 2$)**: Virasoro minimal model $\mathcal{M}(3,4)$, central
charge $c = 1/2$.

---

## $d = 2$ — Exact Critical Exponents

All critical exponents of the 2D Ising universality class are known
**exactly** as rational numbers.

| Exponent | Value | Physical meaning | Derivation | Reference |
|----------|-------|-----------------|------------|-----------|
| $\alpha$ | $0$ | Specific heat: $C \sim \ln\lvert T - T_c\rvert$ | Onsager exact solution | Onsager (1944) Phys. Rev. **65**, 117 |
| $\beta$ | $1/8$ | Order parameter: $M \sim (T_c - T)^\beta$ | Yang magnetization formula | Yang (1952) Phys. Rev. **85**, 808 |
| $\gamma$ | $7/4$ | Susceptibility: $\chi \sim \lvert T - T_c\rvert^{-\gamma}$ | From scaling or exact calculation | Fisher (1964) J. Math. Phys. **5**, 944 |
| $\nu$ | $1$ | Correlation length: $\xi \sim \lvert T - T_c\rvert^{-\nu}$ | Coulomb gas mapping | den Nijs (1979) J. Phys. A **12**, 1857 |
| $\eta$ | $1/4$ | Anomalous dimension: $G(r) \sim r^{-(d-2+\eta)}$ | Two-point function at $T_c$ | Kadanoff (1966) Physics **2**, 263 |
| $\delta$ | $15$ | Critical isotherm: $M \sim h^{1/\delta}$ at $T_c$ | Scaling relation $\delta = 1 + \gamma/\beta$ | — |
| $c$ | $1/2$ | Central charge (CFT) | Minimal model $\mathcal{M}(3,4)$ | BPZ (1984) Nucl. Phys. B **241**, 333 |

### Scaling Relations (Algebraically Exact)

QAtlas stores these values as `Rational{Int}`, so the following
scaling relations hold with **zero floating-point error**:

$$\alpha + 2\beta + \gamma = 2 \qquad \text{(Rushbrooke)}$$
$$\gamma = \beta(\delta - 1) \qquad \text{(Widom)}$$
$$\gamma = \nu(2 - \eta) \qquad \text{(Fisher)}$$
$$2 - \alpha = d\nu \qquad \text{(Josephson, } d = 2\text{)}$$

### CFT Primary Operators

| Operator | Conformal dimension $h$ | Physical correspondence |
|----------|------------------------|------------------------|
| $\mathbb{1}$ (identity) | $0$ | — |
| $\sigma$ (spin field) | $1/16$ | Continuum limit of $\sigma^z$ |
| $\varepsilon$ (energy density) | $1/2$ | Continuum limit of $\sigma^z\sigma^z$ |

The scaling dimensions are $\Delta = h + \bar{h} = 2h$ (diagonal in the
unitary minimal model). The exponents are related to the scaling
dimensions via $\eta = 2\Delta_\sigma - (d - 2) = 2 \times 1/8 = 1/4$
and $\nu = 1/(d - \Delta_\varepsilon) = 1/(2 - 1) = 1$.

---

## $d = 3$ — Conformal Bootstrap

In $d = 3$, the exponents are not known exactly but have been
determined to extraordinary precision via the conformal bootstrap
program.

| Exponent | Value | Uncertainty | Reference |
|----------|-------|-------------|-----------|
| $\alpha$ | $0.11009$ | $\pm 0.00001$ | Kos, Poland, Simmons-Duffin, Vichi (2016) |
| $\beta$ | $0.32642$ | $\pm 0.00001$ | JHEP **08**, 036, Table 2 |
| $\gamma$ | $1.23708$ | $\pm 0.00001$ | " |
| $\delta$ | $4.78984$ | $\pm 0.00001$ | " |
| $\nu$ | $0.62997$ | $\pm 0.00001$ | " |
| $\eta$ | $0.03630$ | $\pm 0.00005$ | " |

These values satisfy the scaling relations approximately (within
the stated uncertainties).

---

## $d \geq 4$ — Mean-Field

Above the upper critical dimension $d_c = 4$, fluctuations are
irrelevant and the exponents take the [mean-field](mean-field.md)
values: $\beta = 1/2$, $\nu = 1/2$, $\gamma = 1$, $\eta = 0$, $\delta = 3$.

---

## Cross-Verification in QAtlas

The 2D Ising exponents are **not just stored values** — they are
cross-checked against independent QAtlas results:

| Exponent | Extracted from | Method | Test file |
|----------|---------------|--------|-----------|
| $\beta = 1/8$ | [Yang $M(T)$](../models/classical/ising-square.md#spontaneous-magnetization-yang) near $T_c$ | log-log slope | `test_universality_cross_check.jl` |
| $\nu z = 1$ | [TFIM gap](../models/quantum/tfim.md#energy-gap-and-quantum-phase-transition) $\Delta(N)$ at $h = J$ | log-log regression | `test_universality_cross_check.jl` |
| $c = 1/2$ | [TFIM entanglement](../models/quantum/tfim.md#central-charge-from-entanglement-entropy) $S(l)$ | Calabrese-Cardy OBC | `test_entanglement_central_charge.jl` |
| $\alpha = 0$ | IsingSquare specific heat near $T_c$ | ForwardDiff on $Z$ | `test_universality_cross_check.jl` |

Each row connects a **universality-level claim** (Source A: CFT) to a
**model-level computation** (Source B: Onsager / Yang / BdG), providing
independent physical validation.

---

## QAtlas API

```julia
# d = 2: exact (Rational{Int})
e = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=2)
# (β = 1//8, ν = 1//1, γ = 7//4, η = 1//4, δ = 15//1, α = 0//1, c = 1//2)

# d = 3: numerical (Float64 + _err)
e3 = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=3)
# (α = 0.11009, α_err = 1.0e-5, β = 0.32642, β_err = 1.0e-5, ...)

# d ≥ 4: mean-field
e4 = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=4)
# → same as fetch(MeanField(), CriticalExponents())

# Backward-compatible alias
e_old = QAtlas.fetch(Ising2D(), CriticalExponents())  # same as d=2
```

---

## Connections

- **Models**: [IsingSquare](../models/classical/ising-square.md),
  [TFIM](../models/quantum/tfim.md)
- **E8 extension**: [E8 mass spectrum](e8.md) — longitudinal field at
  $h = J$ produces the $E_8$ integrable field theory with 8 stable
  particles (Zamolodchikov, 1989).
- **Potts generalization**: [Potts](potts.md) — the Ising model is the
  $q = 2$ case of the $q$-state Potts model.

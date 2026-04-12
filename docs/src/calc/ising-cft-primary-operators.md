# Ising CFT: Primary Operators and Scaling Dimensions

## Setup

The 2D Ising model at the critical point is described by the simplest
non-trivial conformal field theory: the Virasoro minimal model
$\mathcal{M}(3,4)$ with central charge $c = 1/2$.

This CFT has exactly **three primary operators**, fewer than any other
non-trivial unitary CFT. All critical exponents of the
[Ising universality class](../universalities/ising.md) follow from
their scaling dimensions.

## The Three Primary Operators

| Operator | Symbol | Conformal dimension $h = \bar{h}$ | Scaling dimension $\Delta = 2h$ | Physical correspondence |
| -------- | ------ | --------------------------------- | ------------------------------- | ----------------------- |
| Identity | $\mathbb{1}$ | $0$ | $0$ | Trivial |
| Spin field | $\sigma$ | $1/16$ | $1/8$ | Continuum limit of $\sigma^z_i$ |
| Energy density | $\varepsilon$ | $1/2$ | $1$ | Continuum limit of $\sigma^z_i \sigma^z_{i+1}$ |

### Origin from the Kac table

The minimal model $\mathcal{M}(p, p') = \mathcal{M}(3, 4)$ has primary
operators labeled by integers $(r, s)$ with $1 \leq r \leq p-1$,
$1 \leq s \leq p'-1$, with conformal dimensions:

$$h_{r,s} = \frac{(p' r - p s)^2 - (p' - p)^2}{4 p p'}$$

For $\mathcal{M}(3,4)$:

| $(r,s)$ | $h_{r,s}$ | Operator |
| ------- | --------- | -------- |
| $(1,1)$ | $0$ | $\mathbb{1}$ |
| $(1,2)$ | $1/2$ | $\varepsilon$ |
| $(2,1)$ | $1/16$ | $\sigma$ |

(The remaining entries $(2,2)$, $(1,3)$, $(2,3)$ are identified with
these three by the Kac table symmetry $h_{r,s} = h_{p-r, p'-s}$.)

## Connection to Critical Exponents

The critical exponents are determined by the scaling dimensions via
standard CFT relations:

$$\eta = 2\Delta_\sigma - (d - 2) = 2 \times \frac{1}{8} - 0 = \frac{1}{4}$$

$$\nu = \frac{1}{d - \Delta_\varepsilon} = \frac{1}{2 - 1} = 1$$

$$\beta = \frac{\nu \Delta_\sigma}{1} = \frac{1}{8}$$

These reproduce the exact exponents in the
[Ising universality table](../universalities/ising.md).

## Relevance to E8

When the Ising CFT is perturbed by the **spin field** $\sigma$
(corresponding to a longitudinal magnetic field in the TFIM), the
resulting massive field theory is integrable and has the
[E8 mass spectrum](e8-mass-spectrum-derivation.md). See
[Ising CFT magnetic perturbation](ising-cft-magnetic-perturbation.md)
for the mechanism.

When perturbed by the **energy density** $\varepsilon$ (corresponding
to moving away from $T_c$ or tuning $h \neq J$ in the TFIM), the
result is a free massive fermion — no E8 structure.

## References

- A. A. Belavin, A. M. Polyakov, A. B. Zamolodchikov, "Infinite
  conformal symmetry in two-dimensional quantum field theory",
  Nucl. Phys. B **241**, 333 (1984) — BPZ: minimal models.
- P. Di Francesco, P. Mathieu, D. Sénéchal, *Conformal Field Theory*,
  Springer (1997), Ch. 7 — pedagogical treatment of minimal models.

## Used by

- [Ising universality class](../universalities/ising.md) — exponents from CFT
- [E8 universality](../universalities/e8.md) — σ perturbation
- [Ising CFT magnetic perturbation](ising-cft-magnetic-perturbation.md) — next step

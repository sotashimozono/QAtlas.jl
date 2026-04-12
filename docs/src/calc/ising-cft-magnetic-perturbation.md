# Ising CFT + Magnetic Perturbation → Integrable Massive Field Theory

## Setup

Consider the 1D TFIM at the quantum critical point $h = J$, perturbed
by a small longitudinal magnetic field $\lambda$:

$$H = -J\sum_i \sigma^z_i\sigma^z_{i+1} - J\sum_i \sigma^x_i - \lambda\sum_i \sigma^z_i$$

In the continuum (field theory) limit, this becomes the Ising CFT
perturbed by the spin operator $\sigma$:

$$\mathcal{A} = \mathcal{A}_{\text{Ising CFT}} + \lambda \int d^2x\, \sigma(x)$$

The question is: **what is the spectrum of this perturbed theory?**

## Why This Perturbation Is Special

### Relevant vs irrelevant perturbations

The spin field $\sigma$ has scaling dimension $\Delta_\sigma = 1/8 < 2$
(the spacetime dimension), so this perturbation is **relevant**: it
drives the system away from criticality and opens a mass gap.

### Integrability

Zamolodchikov (1989) made the remarkable discovery that this particular
perturbation **preserves integrability**. The key argument:

1. The Ising CFT has an infinite set of conserved charges (from the
   Virasoro algebra).
2. Under the $\sigma$ perturbation, **eight** of these charges survive
   as exact conservation laws of the massive theory.
3. The existence of these conservation laws constrains the S-matrix to
   be purely elastic (no particle production in scattering).
4. The bootstrap equations for the elastic S-matrix have a unique
   solution: the $E_8$ affine Toda field theory.

### Contrast with $\varepsilon$ perturbation

If instead we perturb by the energy density $\varepsilon$
($\Delta_\varepsilon = 1$), the resulting theory is a **free massive
Majorana fermion**. This corresponds to moving $h$ away from $J$ (the
thermal perturbation). No E8 structure appears — the spectrum is a
single particle with mass $m \propto |h - J|$.

## The Two Perturbation Directions

Starting from the Ising critical point:

| Direction | Operator | $\Delta$ | Physical parameter | Result |
| --------- | -------- | -------- | ------------------ | ------ |
| Thermal | $\varepsilon$ | $1$ | $h - J$ (transverse field) | Free massive fermion |
| Magnetic | $\sigma$ | $1/8$ | $\lambda$ (longitudinal field) | **E8 integrable theory** |

The thermal direction is "boring" (free fermion), while the magnetic
direction is "extraordinary" (E8). This asymmetry reflects the
different fusion rules of $\varepsilon$ and $\sigma$ in the Ising CFT.

## Mass Gap Scaling

The mass of the lightest particle scales with the perturbation strength:

$$m_1 \propto |\lambda|^{8/15}$$

The exponent $8/15$ follows from dimensional analysis:
$[\lambda] = [\text{energy}]^{2 - \Delta_\sigma} = [\text{energy}]^{15/8}$,
so $[m] = [\lambda]^{1/(2-\Delta_\sigma)} = [\lambda]^{8/15}$.

The ratios $m_n / m_1$ are universal constants determined by the E8
algebra — see [E8 mass spectrum derivation](e8-mass-spectrum-derivation.md).

## Connection to Experiment

The compound CoNb₂O₆ is a quasi-1D Ising ferromagnet. Near its quantum
critical point (transverse field $B_c \approx 5.5\,\text{T}$), the
inter-chain coupling acts as a weak longitudinal field $\lambda$,
realizing the E8 perturbation. Coldea et al. (2010) observed the first
two particles with $m_2/m_1 = 1.618 \pm 0.015$, consistent with the
golden ratio $\varphi = (1+\sqrt{5})/2 \approx 1.618$.

## References

- A. B. Zamolodchikov, "Integrals of motion and S-matrix of the
  (scaled) $T = T_c$ Ising model with magnetic field", Int. J. Mod.
  Phys. A **4**, 4235 (1989) — discovery of E8 integrability.
- G. Delfino, "Integrable field theory and critical phenomena: the
  Ising model in a magnetic field", J. Phys. A **37**, R45 (2004) —
  comprehensive review.
- R. Coldea et al., "Quantum criticality in an Ising chain:
  experimental evidence for emergent E8 symmetry", Science **327**,
  177 (2010) — experimental confirmation in CoNb₂O₆.

## Used by

- [E8 universality](../universalities/e8.md)
- [E8 mass spectrum derivation](e8-mass-spectrum-derivation.md) — mass ratios
- [TFIM](../models/quantum/tfim.md) — as the unperturbed theory
- [Ising CFT primary operators](ising-cft-primary-operators.md) — operator content

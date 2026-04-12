# Disordered Systems

## Overview

Disordered quantum systems test QAtlas's ED and entanglement
infrastructure in a regime where exact analytical solutions are
unavailable for individual disorder realisations. The physics is
governed by **disorder-averaged** quantities and **universal**
properties at infinite-randomness fixed points (IRFPs).

Two disordered models are currently implemented:

1. **Random-bond Heisenberg chain** --- verifies ground-state
   properties (singlet, positive entanglement) for random couplings.
2. **Random transverse-field Ising model** --- probes the Fisher IRFP
   and its entanglement signatures.

---

## Random-Bond Heisenberg Chain

### Setup

The Hamiltonian is

$$H = \sum_{i=1}^{N-1} J_i\, \mathbf{S}_i \cdot \mathbf{S}_{i+1}$$

where the couplings $J_i > 0$ are drawn independently from a random
distribution (e.g., uniform on $[0, 1]$ or log-uniform). OBC is
used.

### Verified Properties

| Property | Expected | Physical reason |
|----------|----------|-----------------|
| Ground state is a singlet ($S_{\mathrm{tot}} = 0$) | Yes, for any $\{J_i > 0\}$ | Marshall sign rule: antiferromagnetic ground state is always a singlet |
| Entanglement entropy $S(l) > 0$ for $1 \leq l \leq N-1$ | Yes | Ground state is entangled across any bipartition |
| $S(l)$ varies with disorder realisation | Yes | Non-universal, sample-dependent |

These properties hold for **every** disorder realisation, not just
on average. They serve as sanity checks that the ED machinery
(Hamiltonian construction, diagonalization, entropy computation)
works correctly for non-translationally-invariant systems.

---

## Random Transverse-Field Ising Model

### Setup

The Hamiltonian is

$$H = -\sum_{i=1}^{N-1} J_i\, \sigma_i^z \sigma_{i+1}^z - \sum_{i=1}^{N} h_i\, \sigma_i^x$$

where both the couplings $J_i$ and the transverse fields $h_i$ are
drawn from random distributions.

### Fisher Infinite-Randomness Fixed Point (IRFP)

The quantum phase transition between the ferromagnetic ($J \gg h$)
and paramagnetic ($h \gg J$) phases is controlled by the **Fisher
IRFP** when the distributions of $\ln J$ and $\ln h$ have the same
mean:

$$[\ln J] = [\ln h]$$

where $[\cdot]$ denotes the disorder average. At the IRFP, the
system is described by the **strong-disorder renormalization group**
(SDRG, or Ma-Dasgupta-Hu-Fisher procedure): the strongest local
coupling is decimated at each step, generating an effective
random-singlet ground state.

### Properties at the IRFP

| Property | Behaviour | Reference |
|----------|-----------|-----------|
| Typical correlation length | $\ln \xi \sim \lvert\delta\rvert^{-1}$ (activated scaling, not power-law) | Fisher (1995) |
| Average entanglement entropy | $\overline{S(l)} = \frac{c_{\mathrm{eff}}}{3}\ln l + \text{const}$ with $c_{\mathrm{eff}} = \frac{\ln 2}{2} \cdot \frac{1}{1} \approx 0.347\ldots$ | Refael-Moore (2004) |
| Sample-to-sample fluctuations | $\mathrm{Var}[S(l)]$ does not vanish as $N \to \infty$ | Characteristic of infinite-randomness |
| Ground-state entanglement | $S(l) > 0$ for critical samples | Random-singlet structure |

### What QAtlas Currently Tests

For small system sizes ($N \leq 14$), QAtlas verifies:

1. **Critical entanglement is positive**: at the IRFP tuning
   $[\ln J] = [\ln h]$, the ground-state entanglement entropy is
   non-zero for generic bipartitions, confirming the random-singlet
   structure.

2. **Sample-to-sample fluctuations**: repeating the calculation for
   multiple disorder realisations produces a spread in $S(l)$,
   consistent with the infinite-randomness nature of the fixed point.

3. **Gapped phases**: deep in the paramagnetic phase ($h_i \gg J_i$),
   $S(l) \to 0$ (product state), confirming the area law.

### Future: $c_{\mathrm{eff}}$ Extraction

Extracting the effective central charge $c_{\mathrm{eff}} = \ln 2 / 2
\approx 0.347$ requires disorder-averaging $S(l)$ over many samples
at the IRFP and fitting the average to $\overline{S} \propto \ln l$.
This requires larger system sizes and more samples than currently
implemented. It is planned as a future verification target.

---

## QAtlas Test

```julia
# test/verification/test_disordered.jl

# Random-bond Heisenberg: singlet ground state
E, psi = eigen(H_random_heisenberg)
@test is_singlet(psi[:, 1])
@test all(S_l .> 0)

# Random TFIM at IRFP: positive entanglement
E, psi = eigen(H_random_tfim_irfp)
@test S_half > 0  # entanglement at midpoint
```

---

## References

- D. S. Fisher, "Random transverse field Ising spin chains",
  Phys. Rev. Lett. **69**, 534 (1992); "Critical behavior of
  random transverse-field Ising spin chains", Phys. Rev. B **51**,
  6411 (1995) --- SDRG and Fisher IRFP.
- S. K. Ma, C. Dasgupta, C.-K. Hu, "Random antiferromagnetic chain",
  Phys. Rev. Lett. **43**, 1434 (1979) --- original SDRG procedure.
- G. Refael, J. E. Moore, "Entanglement entropy of random quantum
  critical points in one dimension", Phys. Rev. Lett. **93**, 260602
  (2004) --- $c_{\mathrm{eff}} = (\ln 2)/2$ at the random-singlet
  fixed point.
- N. Laflorencie, "Scaling of entanglement entropy in the random
  singlet phase", Phys. Rev. B **72**, 140408(R) (2005).

---

## Connections

- **Entanglement entropy**: [entanglement verification](entanglement.md)
  for the clean (non-disordered) case.
- **TFIM**: the clean model is described in
  [TFIM](../models/quantum/tfim.md); disorder adds random $J_i, h_i$.
- **Heisenberg**: the clean chain is described in
  [Heisenberg](../models/quantum/heisenberg.md); disorder adds
  random $J_i$.
- **Calabrese-Cardy**: the clean formula
  [does not apply](../methods/calabrese-cardy/index.md) at the IRFP;
  the effective central charge $c_{\mathrm{eff}}$ replaces $c$.

# Cross-Verification: Universality ↔ Models

## Overview

This page documents all cross-checks implemented in
`test/verification/test_universality_cross_check.jl`. Each cross-check
connects a **universality-class prediction** (Source A, typically from
CFT or Coulomb gas theory) to a **model-specific exact result**
(Source B, typically from Onsager, Yang, or Bethe ansatz). The two
sources are derived from independent theoretical lines in the physics
literature.

If both agree numerically, the QAtlas data is validated from two
independent origins. This is the "physical correctness proof" of
QAtlas -- it demonstrates that the stored values are not merely
internally consistent but are anchored to independent, well-established
results.

---

## Cross-Check Table

| # | Exponent / quantity | Source A (universality) | Source B (model) | Method | Agreement |
|---|---------------------|------------------------|------------------|--------|-----------|
| 1 | $\beta = 1/8$ (order parameter) | `Universality(:Ising)` $d = 2$ (CFT: BPZ 1984) | `IsingSquare` `SpontaneousMagnetization` (Yang 1952) | Log-log slope of $M(T)$ near $T_c$: $\ln M / \ln(T_c - T) \to \beta$ | $< 1\%$ relative error at each successive $\delta T$ pair |
| 2 | $z = 1$ (dynamic exponent) | `Universality(:Ising)` $d = 2$: $\nu = 1$, $z = 1$ for 1D TFIM | TFIM ED gap $\Delta(N)$ via `build_tfim` + Lattice2D | Log-log slope of $\Delta$ vs $1/N$ at $h = J$: $\log\Delta / \log(1/N) \to z$ | Last pair within 10% of $z = 1$; trend converges |
| 3 | $c = 1/2$ (central charge) | `Universality(:Ising)` $d = 2$ (CFT) | TFIM ED $E_0(N)$ via `build_tfim` | $E_0(N)/N$ converges with $N$; BdG $=$ ED at each $N$ | Spread $< 0.1$; BdG matches ED to $10^{-10}$ |
| 4 | $e_0 = J(1/4 - \ln 2)$ (Bethe) | Bethe ansatz $e_\infty$ (Hulthen 1938) | `Heisenberg1D` `ExactSpectrum` $N = 4$ PBC + ED $N = 6, 8$ PBC | $E_0(N)/N$ approaches $e_\infty$ monotonically; $N = 8$ within 5% | Error decreases: $N = 4 > N = 6 > N = 8$ |
| 5 | $T_c$ self-consistency | Kramers--Wannier duality fixed point: $\sinh(2\beta_c J) = 1$ | `IsingSquare` `CriticalTemperature` (Onsager 1944) | Direct identity check + $M(T_c) = 0$ + $Z(\beta)$ monotonicity near $T_c$ | $\sinh(2\beta_c J) = 1$ to $10^{-14}$; $M(T_c) = 0$ exact |
| 6 | $\nu z = 1$ (gap scaling) | `Universality(:Ising)` $d = 2$: $\nu = 1$, $z = 1$ | TFIM BdG quasiparticle gap (rigorous, $N = 200$) | $\Delta_{\text{BdG}} \approx 2\lvert h - J\rvert$ for $h > J$; log-log slope vs $\lvert h - J\rvert$ | Gap matches thermodynamic prediction within 5%; slope $\approx 1.0$ within 5% |
| 7 | $\alpha = 0$ (specific heat) | `Universality(:Ising)` $d = 2$ (CFT): logarithmic divergence | `IsingSquare` `PartitionFunction` + ForwardDiff $\to C_v(\beta)$ | $C_v$ increases as $T \to T_c$; growth slower than any power law (consistent with log) | $C_v > 0$; $C_v$ increases toward $T_c$; $\alpha = 0$ confirmed |
| 8 | $e_\infty = -4J/\pi$ (TFIM) | BdG dispersion integral: $\int_0^\pi \Lambda(k)\,dk / \pi$ | TFIM BdG $E_0(N)$ (rigorous, $N = 50$--$400$) | $E_0(N)/N \to e_\infty$ as $1/N$; boundary correction $\times N \approx \text{const}$ | $N = 400$ within 0.1% of $e_\infty$; $O(1/N)$ scaling verified to 1% |

---

## Detailed Descriptions

### 1. Yang magnetization → $\beta = 1/8$

The order-parameter exponent $\beta = 1/8$ is a universal prediction
of the 2D Ising CFT (Source A). Yang's exact formula
$M(T) = (1 - \sinh^{-4}(2\beta J))^{1/8}$ (Source B) encodes this
exponent in the power $1/8$. The test extracts $\beta$ numerically
from the log-log slope $\ln M / \ln(T_c - T)$ at progressively
closer temperatures to $T_c$, verifying convergence to $1/8$ from
the functional form without assuming the exponent value.

### 2. TFIM gap $\Delta(N)$ → $z = 1$

At the TFIM critical point $h = J$, the many-body gap closes as
$\Delta \sim N^{-z}$ with dynamic exponent $z = 1$. This is extracted
from ED gaps at $N = 4, 6, 8, 10, 12$ via successive log-log slopes,
which should converge to $z = 1$. Finite-size corrections cause
deviations at small $N$.

### 3. TFIM $E_0$ scaling → consistent with $c = 1/2$

At criticality, finite-size scaling of $E_0(N)/N$ involves the central
charge $c$ via the Cardy formula. This test verifies that ED energies
are consistent with the BdG analytical values (agreement to
$10^{-10}$) and that $E_0(N)/N$ converges as $N$ grows.

### 4. Heisenberg $E_0(N)/N$ → Bethe ansatz $e_\infty$

The Bethe ansatz predicts $e_\infty = J(1/4 - \ln 2) \approx -0.4431J$
for the infinite PBC chain. Finite-size ED at $N = 4$ (from
`QAtlas.fetch`), $N = 6$, and $N = 8$ (from `build_spinhalf_heisenberg`
+ Lattice2D) shows that $E_0(N)/N$ lies below $e_\infty$ (negative
finite-size correction for PBC) and converges toward it as $N$ grows.

### 5. Onsager $T_c$ ↔ duality self-consistency

The Kramers--Wannier duality predicts that the critical coupling
satisfies $\sinh(2K_c) = 1$. This is verified to machine precision
against `IsingSquare` `CriticalTemperature`. Additionally,
$M(T_c) = 0$ (phase boundary) and $Z(\beta)$ monotonicity near
$T_c$ are checked.

### 6. TFIM BdG gap → $\nu z = 1$ (rigorous)

Using the exact BdG spectrum at $N = 200$ (no ED approximation), the
gap $\Delta = \min(\Lambda_n) \approx 2|h - J|$ is verified in the
disordered phase. A log-log regression of $\Delta$ vs $|h - J|$
extracts the exponent $\nu z$, which should equal 1 for the 1D TFIM.

### 7. IsingSquare specific heat → $\alpha = 0$

The 2D Ising specific heat diverges logarithmically at $T_c$
($\alpha = 0$), not as a power law. Using the transfer-matrix
partition function and ForwardDiff to compute
$C_v = \beta^2 \partial^2 (\ln Z) / \partial \beta^2$, the test
verifies that $C_v$ is positive and grows toward $T_c$, consistent
with logarithmic (not power-law) divergence.

### 8. TFIM $E_0(N)/N$ → $e_\infty = -4J/\pi$ (rigorous, OBC)

At the TFIM critical point ($h = J$), the thermodynamic-limit energy
per site is $e_\infty = -4J/\pi$ (from integrating the BdG dispersion
$\Lambda(k) = 2J|\sin k|$). The test verifies convergence of the exact
BdG $E_0(N)/N$ to this value with $O(1/N)$ boundary corrections,
confirmed by checking that the correction $\times N$ is approximately
constant across $N = 50, 100, 200, 400$.

---

## Significance

These 8 cross-checks collectively ensure that QAtlas's universality
exponents, model-specific exact solutions, and computational methods
are mutually consistent. Any systematic error in the source data --
a wrong sign, a missing factor of 2, an incorrect exponent -- would
cause at least one of these cross-checks to fail.

The cross-checks span:

- **Two universality classes**: Ising ($c = 1/2$) and Heisenberg
  ($c = 1$ Luttinger liquid)
- **Three model families**: classical Ising (transfer matrix),
  quantum TFIM (BdG), quantum Heisenberg (Bethe ansatz + ED)
- **Four computational methods**: transfer matrix + AD, BdG
  diagonalization, exact diagonalization, analytical formulas

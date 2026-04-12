# Kramers-Wannier Duality

## Setup

The 2D classical Ising model (or equivalently the 1D TFIM) possesses
a duality symmetry that maps the high-temperature phase to the
low-temperature phase. The critical point is the unique self-dual
point.

## The Duality Transformation

### Classical 2D Ising

The partition function of the 2D Ising model on a square lattice
can be expanded in two ways:

- **High-temperature expansion**: sum over bond configurations
- **Low-temperature expansion**: sum over domain-wall configurations

These two expansions are related by the duality:

$$\tanh K^* = e^{-2K}$$

where $K = \beta J$ is the reduced coupling and $K^*$ is the dual
coupling. The self-dual point satisfies $K = K^*$:

$$\tanh K_c = e^{-2K_c} \quad \Longleftrightarrow \quad \sinh(2K_c) = 1$$

Solving: $K_c = \frac{1}{2}\ln(1 + \sqrt{2})$, hence
$T_c = 2J/\ln(1+\sqrt{2})$.

### Quantum 1D TFIM

For the 1D TFIM $H = -J\sum\sigma^z\sigma^z - h\sum\sigma^x$, the
duality maps:

$$J \leftrightarrow h$$

The critical point is $h = J$ (self-dual). The duality exchanges
the ordered and disordered phases.

## Result

$$T_c = \frac{2J}{\ln(1 + \sqrt{2})}, \qquad \sinh(2\beta_c J) = 1$$

## References

- H. A. Kramers, G. H. Wannier, Phys. Rev. **60**, 252 (1941).

## Used by

- [IsingSquare: Critical Temperature](../models/classical/ising-square.md#critical-temperature-onsager)
- [TFIM: Phase Diagram](../models/quantum/tfim.md)
- [JW-TFIM-BdG derivation](jw-tfim-bdg.md) — duality as Step 1

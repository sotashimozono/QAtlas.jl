# Heisenberg Chain

## Overview

The spin-1/2 antiferromagnetic Heisenberg chain is the paradigmatic
model of quantum magnetism. Unlike the [TFIM](tfim.md), it is not
mappable to free fermions; its exact solution requires the Bethe ansatz
(1931), one of the earliest and most celebrated applications of
integrability in condensed matter physics.

$$H = J \sum_{i} \mathbf{S}_i \cdot \mathbf{S}_{i+1}$$

where $\mathbf{S}_i = \tfrac{1}{2}\boldsymbol{\sigma}_i$ are spin-1/2
operators and $J > 0$ is the antiferromagnetic exchange coupling.

**Parameters**: Exchange coupling $J$ (default 1.0; $J > 0$ AFM,
$J < 0$ FM).

**Key physics**: The 1D Heisenberg chain is gapless with algebraically
decaying spin correlations. At low energies it is described by a
$c = 1$ conformal field theory (free boson / Luttinger liquid), in
contrast to the $c = 1/2$ Ising CFT governing the [TFIM](tfim.md).
The ground state is a spin singlet ($S_{\text{tot}} = 0$) for any
finite even $N$ with AFM coupling.

---

## Dimer Spectrum (N = 2, OBC)

### Statement

For two spin-1/2 sites coupled by $H = J\,\mathbf{S}_1 \cdot \mathbf{S}_2$,
the Hilbert space $\mathbb{C}^2 \otimes \mathbb{C}^2$ decomposes into
a singlet ($S_{\text{tot}} = 0$) and a triplet ($S_{\text{tot}} = 1$):

$$\text{Spectrum} = \left\{-\frac{3J}{4},\; \frac{J}{4},\; \frac{J}{4},\; \frac{J}{4}\right\}$$

The singlet--triplet gap is $\Delta = J$.

### Derivation

Rewriting the Hamiltonian via the total-spin identity

$$\mathbf{S}_1 \cdot \mathbf{S}_2 = \frac{1}{2}\!\left[S_{\text{tot}}(S_{\text{tot}}+1) - \frac{3}{2}\right]$$

immediately gives $E_s = -3J/4$ for $S_{\text{tot}} = 0$ and
$E_t = J/4$ for $S_{\text{tot}} = 1$. Full details are in
**[Heisenberg Dimer: Singlet-Triplet](../../calc/heisenberg-dimer-singlet-triplet.md)**.

### References

- A. Auerbach, *Interacting Electrons and Quantum Magnetism*
  (Springer, 1994), Section 2 -- textbook derivation.
- Any quantum mechanics textbook (e.g., Sakurai, Ch. 4).

### QAtlas API

```julia
# Full 4-state spectrum of the dimer
λ = QAtlas.fetch(Heisenberg1D(), ExactSpectrum(); N=2, J=1.0, bc=:OBC)
# → [-0.75, 0.25, 0.25, 0.25]
```

### Verification

| Test file | Method | What is checked |
|-----------|--------|-----------------|
| `test_heisenberg_dimer.jl` | Lattice2D ED ($2 \times 1$ OBC chain) | $\lambda_{\text{ED}} = \lambda_{\text{exact}}$ for $J = 1.0, 0.5, 2.0, -1.0$ |
| `test_heisenberg_dimer.jl` | Analytical | Singlet--triplet gap $\Delta = J$ |
| `test_heisenberg_dimer.jl` | Structural | $\text{tr}\,H = 0$, ground-state wavefunction is the antisymmetric singlet |

---

## 4-Site PBC Ring (N = 4)

### Statement

For $N = 4$ spins on a periodic ring with Hamiltonian
$H = J \sum_{i=1}^{4} \mathbf{S}_i \cdot \mathbf{S}_{i+1}$
($\mathbf{S}_5 \equiv \mathbf{S}_1$), the $2^4 = 16$ eigenvalues are:

| Energy | Degeneracy | Spin sector |
|--------|------------|-------------|
| $-2J$ | 1 | Singlet ($S = 0$) |
| $-J$ | 3 | Triplet ($S = 1$) |
| $0$ | 7 | Mixed ($1 \times S = 0 + 2 \times S = 1$) |
| $+J$ | 5 | Quintet ($S = 2$) |

The ground state $E_0 = -2J$ is a unique singlet. The ferromagnetic
sector ($|\!\uparrow\uparrow\uparrow\uparrow\rangle$ etc.) sits at
$E = +J$, since each of the 4 bonds contributes $S^z_i S^z_{i+1} = 1/4$.

### Derivation

The 4-site ring can be solved by direct diagonalization of the
$16 \times 16$ Hamiltonian or by exploiting the full SU(2) symmetry
and translational invariance. The result is also obtainable as a
special case of the Bethe ansatz for $N = 4$.

### References

- A. Auerbach, *Interacting Electrons and Quantum Magnetism*
  (Springer, 1994), Section 2.3.
- H. Bethe, Z. Physik **71**, 205 (1931) -- original Bethe ansatz
  (applicable to any $N$).

### QAtlas API

```julia
# Full 16-state spectrum of the 4-site PBC ring
λ = QAtlas.fetch(Heisenberg1D(), ExactSpectrum(); N=4, J=1.0, bc=:PBC)
# → [-2.0, -1.0, -1.0, -1.0, 0.0, ..., 1.0, 1.0, 1.0, 1.0, 1.0]
```

### Verification

| Test file | Method | What is checked |
|-----------|--------|-----------------|
| `test_heisenberg_4site_pbc.jl` | Lattice2D ED ($4 \times 1$ PBC chain) | $\lambda_{\text{ED}} = \lambda_{\text{exact}}$ for $J = 1.0, 0.5, 2.0, -1.0$ |
| `test_heisenberg_4site_pbc.jl` | Degeneracy counting | $\{1, 3, 7, 5\}$ at $\{-2J, -J, 0, +J\}$ |
| `test_heisenberg_4site_pbc.jl` | Structural | $\text{tr}\,H = 0$, quintet eigenvalue via $\langle\uparrow\uparrow\uparrow\uparrow|H|\uparrow\uparrow\uparrow\uparrow\rangle = J$ |
| `test_heisenberg_4site_pbc.jl` | Finite-size | $E_0/N = -0.5 < e_{\infty} \approx -0.443$ (overshoot from finite-size correction) |

---

## Bethe Ansatz Ground-State Energy Density

### Statement

In the thermodynamic limit ($N \to \infty$, PBC), the ground-state
energy per site of the spin-1/2 AFM Heisenberg chain is

$$e_0 = J\!\left(\frac{1}{4} - \ln 2\right) \approx -0.4431\,J$$

This exact result was first obtained by Hulthen (1938) from the Bethe
ansatz solution.

### Derivation

The Bethe ansatz parameterizes eigenstates in the $M$-down-spin sector
by a set of rapidities $\{\lambda_j\}$ satisfying the Bethe equations.
In the thermodynamic limit, the ground-state rapidity distribution
satisfies a linear integral equation whose solution yields the energy
through integration. The full calculation is given in
**[Bethe Ansatz: Heisenberg $e_0$](../../calc/bethe-ansatz-heisenberg-e0.md)**.

### Finite-size corrections

For a PBC chain of $N$ sites:

$$\frac{E_0(N)}{N} = e_0 + O\!\left(\frac{1}{N^2}\right)$$

with logarithmic corrections also present. Due to the negative sign
of the leading finite-size correction (PBC), $E_0(N)/N < e_0$ for
finite $N$.

### References

- H. Bethe, "Zur Theorie der Metalle. I. Eigenwerte und Eigenfunktionen
  der linearen Atomkette", Z. Physik **71**, 205 (1931) -- original
  Bethe ansatz.
- L. Hulthen, "Uber das Austauschproblem eines Kristalles",
  Ark. Mat. Astron. Fys. **26A**, No. 11, 1 (1938) -- first evaluation
  of $e_0 = 1/4 - \ln 2$.

### QAtlas API

```julia
# Exact thermodynamic-limit ground-state energy per site
e₀ = QAtlas.fetch(Heisenberg1D(), GroundStateEnergyDensity(); J=1.0)
# → -0.44314718055994530
```

### Verification

| Test file | Method | What is checked |
|-----------|--------|-----------------|
| `test_bethe_ansatz.jl` | Numerical value | $e_0 \approx J(1/4 - \ln 2)$ to machine precision |
| `test_bethe_ansatz.jl` | $J$ scaling | $e_0(J) = J \cdot e_0(1)$ |
| `test_bethe_ansatz.jl` | Consistency | $E_0(N=4)/4 < e_0$ (finite-size overshoot) |
| `test_universality_cross_check.jl` | ED at $N = 4, 6, 8$ PBC | $|E_0(N)/N - e_0|$ decreases with $N$; $N = 8$ within 5% |

---

## Connections

- **Universality**: The low-energy physics of the Heisenberg chain is
  described by a $c = 1$ Luttinger liquid (free boson CFT), **not** the
  $c = 1/2$ Ising CFT. This is a fundamentally different universality
  class from the [TFIM](tfim.md) / [IsingSquare](../classical/ising-square.md).
- **Entanglement verification**: The Heisenberg chain provides a
  $c = 1$ test case for the
  [Calabrese--Cardy formula](../../methods/calabrese-cardy/index.md),
  complementing the $c = 1/2$ test from the TFIM.
- **Cross-verification**: The Bethe ansatz $e_0$ is cross-checked against
  finite-size ED in
  [test_universality_cross_check.jl](../../verification/cross-checks.md),
  constituting an independent physical verification from two different
  theoretical lines (Bethe ansatz vs exact diagonalization).

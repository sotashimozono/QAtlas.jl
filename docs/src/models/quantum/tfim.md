# Transverse-Field Ising Model (TFIM)

## Overview

The one-dimensional transverse-field Ising model is one of the simplest
quantum many-body systems exhibiting a quantum phase transition. It is
exactly solvable via [Jordan-Wigner transformation](../../methods/jordan-wigner/index.md)
followed by Bogoliubov-de Gennes (BdG) diagonalization.

$$H = -J \sum_{i} \sigma^z_i \sigma^z_{i+1} - h \sum_{i} \sigma^x_i$$

**Parameters**: Ising coupling $J$ (default 1.0), transverse field $h$.

**Phase diagram**:

- $h/J < 1$: Ferromagnetic ordered phase ($\langle \sigma^z \rangle \neq 0$)
- $h/J = 1$: Quantum critical point ([Ising CFT, $c = 1/2$](../../universalities/ising.md))
- $h/J > 1$: Quantum paramagnetic phase ($\langle \sigma^z \rangle = 0$)

**Universality**: The critical point belongs to the
[2D Ising universality class](../../universalities/ising.md) via the
quantum-classical mapping (1+1D quantum ↔ 2D classical).

---

## Boundary Conditions

QAtlas supports three boundary conditions for the TFIM, each with
different physical content:

| BC       | `fetch` argument         | BdG size                   | Physical setting                        |
| -------- | ------------------------ | -------------------------- | --------------------------------------- |
| OBC      | `OBC()`                  | $2N \times 2N$ (numerical) | Open chain, $N$ sites, $N-1$ bonds      |
| PBC      | (not directly supported) | —                          | Ring of $N$ sites, $N$ bonds            |
| Infinite | `Infinite()`             | $k$-integral               | Thermodynamic limit, PBC $N \to \infty$ |

**OBC**: The BdG matrix is diagonalized numerically. Boundary effects
include the Z₂ tunneling splitting in the ordered phase and the
boundary energy correction $O(1/N)$ at criticality. See
[gap analysis](#energy-gap-and-quantum-phase-transition) for details.

**Infinite**: The quasiparticle dispersion
$\Lambda(k) = 2\sqrt{J^2 + h^2 - 2Jh\cos k}$ is integrated over
the Brillouin zone using Gauss-Kronrod quadrature (QuadGK.jl).

---

## Ground-State Energy

### Statement

The ground-state energy of the OBC TFIM with $N$ sites is

$$E_0 = -\sum_{n=1}^{N} \frac{\Lambda_n}{2}$$

where $\{\Lambda_n\}$ are the positive eigenvalues of the $2N \times 2N$
BdG matrix. At finite temperature $\beta = 1/(k_B T)$:

$$\langle H \rangle(\beta) = -\sum_{n=1}^{N} \frac{\Lambda_n}{2} \tanh\!\left(\frac{\beta \Lambda_n}{2}\right)$$

### Derivation

The TFIM is solved exactly via the
[Jordan-Wigner transformation](../../methods/jordan-wigner/index.md),
which maps the spin chain to free fermions after a
[Kramers-Wannier duality](../../calc/kramers-wannier-duality.md) step.
The full derivation — including why the duality is needed for the
$\sigma^z\sigma^z$ convention and the explicit construction of the
BdG matrix — is given in the calculation note
**[JW-TFIM-BdG](../../calc/jw-tfim-bdg.md)**.

The result is a $2N \times 2N$ real symmetric BdG matrix whose
eigenvalues come in $\pm\Lambda_n$ pairs. The positive eigenvalues
$\Lambda_n > 0$ are the quasiparticle energies, and the total energy
at inverse temperature $\beta$ is:

$$\langle H \rangle = -\sum_n \frac{\Lambda_n}{2} \tanh\!\left(\frac{\beta \Lambda_n}{2}\right)$$

!!! note "Thermodynamic limit"
For PBC in the $N \to \infty$ limit, the quasiparticle dispersion
is $\Lambda(k) = 2\sqrt{J^2 + h^2 - 2Jh\cos k}$, and the energy
per site becomes a $k$-integral evaluated by Gauss-Kronrod
quadrature (QuadGK.jl).

### References

- P. Pfeuty, "The one-dimensional Ising model with a transverse field",
  Ann. Phys. **57**, 79 (1970) — exact solution of the 1D TFIM.
- E. Lieb, T. Schultz, D. Mattis, "Two Soluble Models of an
  Antiferromagnetic Chain", Ann. Phys. **16**, 407 (1961) — JW
  transformation for spin chains.
- S. Sachdev, _Quantum Phase Transitions_, Cambridge University Press
  (2011), Ch. 5 — pedagogical treatment.

### QAtlas API

```julia
# Ground-state energy (β → ∞), OBC, N=16
E₀ = QAtlas.fetch(:TFIM, :energy, OBC(); N=16, J=1.0, h=0.5)

# Finite-temperature energy
E_β = QAtlas.fetch(:TFIM, :energy, OBC(); N=16, J=1.0, h=0.5, beta=2.0)

# Thermodynamic limit (PBC, N→∞)
ε = QAtlas.fetch(:TFIM, :energy, Infinite(); J=1.0, h=0.5, beta=2.0)
```

### Verification

| Test file                          | Method                    | What is checked                                        |
| ---------------------------------- | ------------------------- | ------------------------------------------------------ |
| `test_tfim_gap_closure.jl`         | Dense ED via `build_tfim` | $E_0^{\text{ED}} = E_0^{\text{BdG}}$ for $N = 4, 6, 8$ |
| `test_universality_cross_check.jl` | BdG at $N = 200$          | $E_0/N \to -4J/\pi$ at $h = J$                         |

---

## Finite-Temperature Observables

### Statement

At inverse temperature $\beta$ and for $N$ sites (OBC), the following
quantities are computed from the BdG spectrum $\{\Lambda_n\}$:

| Quantity      | Formula                                                                                  | Symbol           |
| ------------- | ---------------------------------------------------------------------------------------- | ---------------- |
| Free energy   | $F = -\frac{1}{\beta}\sum_n \ln\left[2\cosh\left(\frac{\beta\Lambda_n}{2}\right)\right]$ | `:free_energy`   |
| Entropy       | $S = \beta(\langle H \rangle - F)$                                                       | `:entropy`       |
| Specific heat | $C_v = -\beta^2 \frac{\partial \langle H \rangle}{\partial \beta}$                       | `:specific_heat` |

### Derivation

All quantities follow from the free-fermion partition function. For
independent modes with energies $\Lambda_n$:

$$\mathcal{Z} = \prod_n \left[2\cosh\!\left(\frac{\beta\Lambda_n}{2}\right)\right]$$

The free energy is $F = -\beta^{-1}\ln\mathcal{Z}$, and all other
thermodynamic quantities are obtained by appropriate $\beta$-derivatives.

### References

- S. Sachdev, _Quantum Phase Transitions_ (2011), Ch. 5.3.
- QAtlas: `src/models/quantum/TFIM_thermal.jl` — full implementation.

### QAtlas API

```julia
F = QAtlas.fetch(:TFIM, :free_energy, OBC(); N=16, J=1.0, h=0.5, beta=2.0)
S = QAtlas.fetch(:TFIM, :entropy, OBC(); N=16, J=1.0, h=0.5, beta=2.0)
Cv = QAtlas.fetch(:TFIM, :specific_heat, OBC(); N=16, J=1.0, h=0.5, beta=2.0)
```

### Verification

| Test file              | Method                 | What is checked                                       |
| ---------------------- | ---------------------- | ----------------------------------------------------- |
| `test_TFIM_thermal.jl` | Dense ED ($N \leq 10$) | Exact match of $F$, $S$, $C_v$ against independent ED |

---

## Energy Gap and Quantum Phase Transition

### Statement

The many-body energy gap $\Delta = E_1 - E_0$ equals the smallest BdG
quasiparticle energy $\Lambda_{\min}$. In the thermodynamic limit:

$$\Delta = 2|J - h|$$

At the critical point $h = J$, the gap closes as $\Delta \sim N^{-z}$
with dynamic exponent $z = 1$.

### Physical Context

- **Ordered phase** ($h < J$): for OBC with finite $N$, the "gap" seen
  by exact diagonalization is actually the Z₂ **tunneling splitting**
  between $|\!\uparrow\cdots\uparrow\rangle$ and
  $|\!\downarrow\cdots\downarrow\rangle$, which is exponentially small
  in $N$. This is distinct from the physical excitation gap
  $\Delta \approx 2(J - h)$.

- **Critical point** ($h = J$): $\Delta \sim \pi/(N)$ (finite-size
  gap for OBC).

- **Disordered phase** ($h > J$): $\Delta \approx 2(h - J)$, i.e.,
  the paramagnetic gap.

### References

- P. Pfeuty, Ann. Phys. **57**, 79 (1970), Eq. (3.6).
- S. Sachdev, _Quantum Phase Transitions_ (2011), §5.5.

### QAtlas API

```julia
# Infinite chain — closed form Δ = 2|h − J|
QAtlas.fetch(TFIM(; J=1.0, h=0.3), MassGap(), Infinite())
# → 1.4

QAtlas.fetch(TFIM(; J=1.0, h=1.0), MassGap(), Infinite())
# → 0.0   (critical: gap closes)

# OBC finite-N — smallest positive BdG eigenvalue
QAtlas.fetch(TFIM(; J=1.0, h=1.0), MassGap(), OBC(32))
# → π J / N ≈ 0.098   (finite-size CFT gap)
```

Symbol aliases recognised by the legacy layer include `:mass_gap`,
`:gap`, `:Δ`, `:excitation_gap`, and the capitalised `:MassGap`.

### Verification

| Test file                          | Method                  | What is checked                                           |
| ---------------------------------- | ----------------------- | --------------------------------------------------------- | --- | -------------------------------------- |
| `test_tfim_gap_closure.jl`         | Dense ED ($N = 4$–$12$) | Gap shrinks with $N$ at $h = J$                           |
| `test_tfim_gap_closure.jl`         | ED                      | Ordered-phase gap is Z₂ tunneling ($< 10^{-3}$ for $N=6$) |
| `test_universality_cross_check.jl` | BdG ($N = 200$)         | $\Delta \approx 2                                         | h-J | $; $\nu z = 1$ from log-log regression |

---

## Central Charge from Entanglement Entropy

### Statement

At the critical point $h = J$, the entanglement entropy of a
subsystem of $l$ consecutive sites in an $N$-site OBC chain obeys
the [Calabrese-Cardy formula](../../methods/calabrese-cardy/index.md):

$$S(l) = \frac{c}{6}\ln\!\left[\frac{2N}{\pi}\sin\!\left(\frac{\pi l}{N}\right)\right] + s_1$$

with central charge $c = 1/2$ (Ising CFT). See the
[Calabrese-Cardy method page](../../methods/calabrese-cardy/index.md)
for the general formula, OBC vs PBC prefactors, and extraction
procedure.

### Physical Context

!!! note "Peschel method limitation"
The TFIM is a free-fermion system, so in principle $S(l)$ can be
computed in $O(N^3)$ via the Peschel correlation-matrix method.
However, QAtlas's $\sigma^z\sigma^z$ convention maps to free fermions
only after [Kramers-Wannier duality](../../calc/kramers-wannier-duality.md),
and the BdG eigenvectors give **dual-fermion** correlators, not
spin-basis correlators. QAtlas therefore uses full
[exact diagonalization](../../methods/exact-diagonalization/index.md)
for entanglement entropy.

### References

- P. Calabrese, J. Cardy, J. Stat. Mech. **0406**, P06002 (2004), Eq. (19).
- I. Peschel, J. Phys. A **36**, L205 (2003).

### QAtlas API

Computed via ED (not `fetch` — test utility only):

```julia
include("test/util/spinhalf_ed.jl")
lat = build_lattice(Square, N, 1; boundary=OpenAxis())
H = build_tfim(lat, J, J)
ψ₀ = eigen(Symmetric(H)).vectors[:, 1]
S_l = entanglement_entropy(ψ₀, l, N)
```

### Verification

| Test file                             | Method             | What is checked                               |
| ------------------------------------- | ------------------ | --------------------------------------------- |
| `test_entanglement_central_charge.jl` | ED ($N = 10$–$16$) | $c_{\text{extracted}} \approx 0.5$ within 10% |
| `test_entanglement_central_charge.jl` | ED                 | $S(l)$ symmetric, maximal at $l = N/2$        |
| `test_entanglement_central_charge.jl` | ED ($h \gg J$)     | Area law: $S \approx 0$ in disordered phase   |

---

## Connections

- **Universality**: [Ising universality class](../../universalities/ising.md) —
  the TFIM critical point has $c = 1/2$, $\nu = 1$, $z = 1$.
- **Classical counterpart**: [IsingSquare](../classical/ising-square.md) —
  the 1+1D TFIM maps to the 2D classical Ising model via the
  quantum-classical correspondence ($\beta_{\text{classical}} \leftrightarrow$ imaginary time).
- **Disordered version**: [Random TFIM](../../verification/disordered.md) —
  the Fisher IRFP at $[\ln J]_{\text{avg}} = [\ln h]_{\text{avg}}$.
- **E8 spectrum**: [E8 universality](../../universalities/e8.md) —
  adding a longitudinal field $\lambda \sigma^z$ at $h = J$ breaks
  integrability and produces the E8 mass spectrum.

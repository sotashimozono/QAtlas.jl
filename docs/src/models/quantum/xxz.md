# XXZ Chain (1D, spin-1/2)

## Overview

The spin-1/2 XXZ chain generalises both the [Heisenberg
chain](heisenberg.md) (isotropic $\Delta = 1$) and the XX / free-fermion
point ($\Delta = 0$) through a single anisotropy parameter $\Delta$:

$$H = J \sum_{i} \bigl[\, S^x_i S^x_{i+1} + S^y_i S^y_{i+1}
                       + \Delta\, S^z_i S^z_{i+1} \,\bigr]$$

with $\mathbf{S}_i = \tfrac{1}{2}\boldsymbol{\sigma}_i$ and $J > 0$ the
antiferromagnetic exchange coupling.

**Parameters**: $J$ (exchange coupling, default `1.0`), $\Delta$
(anisotropy, default `0.0` = XX point).

**Phase diagram**:

| Regime | Phase | Description |
|--------|-------|-------------|
| $\Delta < -1$ | Gapped ferromagnet | Ising-like FM order |
| $\Delta = -1$ | Saturated ferromagnet | $e_0/J = -1/4$ |
| $-1 < \Delta \le 1$ | Luttinger liquid ($c = 1$) | Gapless, critical |
| $\Delta = 0$ | XX / free fermion | $e_0/J = -1/\pi$ |
| $\Delta = 1$ | Isotropic AF Heisenberg | $e_0/J = 1/4 - \ln 2$ |
| $\Delta > 1$ | Gapped Néel antiferromagnet | Ising-like AFM order |

---

## Ground-state energy density (infinite chain)

### Statement

QAtlas currently exposes the three canonical exact points:

$$
e_0/J =
\begin{cases}
-\,\tfrac{1}{4} & (\Delta = -1,\ \text{FM saturated}) \\[2pt]
-\,\tfrac{1}{\pi} & (\Delta = 0,\ \text{XX / free fermion}) \\[2pt]
\tfrac{1}{4} - \ln 2 & (\Delta = 1,\ \text{AF Heisenberg, Hulthén 1938})
\end{cases}
$$

For every other $\Delta$ the call returns `NaN` and emits a warning —
the general-$\Delta$ Yang–Yang integral has multiple inequivalent
normalisations in the literature and is tracked as a follow-up PR.

### References

- L. Hulthén, Ark. Mat. Astron. Fys. **26A**, No. 11 (1938) — evaluation
  at $\Delta = 1$.
- C. N. Yang, C. P. Yang, Phys. Rev. **150**, 321 (1966) — general
  Bethe-ansatz integral equation.
- M. Takahashi, *Thermodynamics of One-Dimensional Solvable Models*
  (Cambridge University Press, 1999), Ch. 4.

### QAtlas API

```julia
e_xx = QAtlas.fetch(XXZ1D(; J=1.0, Δ=0.0), Energy(), Infinite())
# → -0.3183098861837907          (= -1/π)

e_af = QAtlas.fetch(XXZ1D(; J=1.0, Δ=1.0), Energy(), Infinite())
# → -0.4431471805599453          (= 1/4 - ln 2)

e_fm = QAtlas.fetch(XXZ1D(; J=1.0, Δ=-1.0), Energy(), Infinite())
# → -0.25
```

The `GroundStateEnergyDensity()` quantity is an alias returning the
same value.

### Verification

| Test | Method | What is checked |
|------|--------|-----------------|
| `test_XXZ1D.jl` | Analytical | $e_0(\Delta = 0) = -1/\pi$ to $10^{-10}$ |
| `test_XXZ1D.jl` | Analytical | $e_0(\Delta = 1) = 1/4 - \ln 2$ to $10^{-10}$ |
| `test_XXZ1D.jl` | Analytical | $e_0(\Delta = -1) = -1/4$ to $10^{-14}$ |
| `test_XXZ1D.jl` | Warning behaviour | `NaN` + `general-Δ` warning for any other $\Delta$ |

---

## Central charge (critical regime only)

### Statement

For $-1 < \Delta < 1$ the chain flows to a $c = 1$ compactified-boson
(Luttinger-liquid) CFT in the IR:

$$c(\Delta) = 1, \qquad -1 < \Delta < 1.$$

Outside this window the chain is gapped; QAtlas returns `NaN` with a
warning.

### QAtlas API

```julia
QAtlas.fetch(XXZ1D(; Δ=0.3), CentralCharge(), Infinite())  # → 1.0
QAtlas.fetch(XXZ1D(; Δ=1.5), CentralCharge(), Infinite())  # → NaN (+ warn)
```

---

## Luttinger parameter $K$

### Statement

Across the full critical regime $-1 < \Delta \le 1$,

$$\boxed{\,K(\Delta) = \frac{\pi}{2\,(\pi - \gamma)},
         \qquad \gamma \equiv \arccos \Delta\,}$$

Canonical values:

| $\Delta$ | $\gamma$ | $K$ | Interpretation |
|----------|----------|-----|-----------------|
| $-1^{+}$ | $\pi^{-}$ | $\to \infty$ | FM boundary |
| $0$ | $\pi/2$ | $1$ | XX / free fermion |
| $1$ | $0$ | $1/2$ | AF Heisenberg |

Monotone decreasing in $\Delta$.

Full derivation: **[XXZ Luttinger parameters from Bethe ansatz
](../../calc/xxz-luttinger-parameters.md)**.

### References

- T. Giamarchi, *Quantum Physics in One Dimension* (Oxford, 2004),
  Ch. 6.
- F. D. M. Haldane, Phys. Rev. Lett. **45**, 1358 (1980);
  Phys. Rev. Lett. **47**, 1840 (1981) — bosonisation of the XXZ chain.

### QAtlas API

```julia
QAtlas.fetch(XXZ1D(; Δ=0.0), LuttingerParameter(), Infinite())  # → 1.0
QAtlas.fetch(XXZ1D(; Δ=1.0), LuttingerParameter(), Infinite())  # → 0.5
```

---

## Luttinger / spin-wave velocity $u$

### Statement

$$\boxed{\,u(\Delta) = J\cdot \frac{\pi}{2}\,\frac{\sin\gamma}{\gamma},
         \qquad \gamma \equiv \arccos \Delta\,}$$

Canonical values:

| $\Delta$ | $u / J$ | Identification |
|----------|---------|----------------|
| $0$ | $1$ | Free-fermion Fermi velocity $v_F$ |
| $1$ | $\pi/2$ | des Cloizeaux–Pearson spin-wave velocity |

`SpinWaveVelocity` is a **type-level alias** of `LuttingerVelocity`
(`const SpinWaveVelocity = LuttingerVelocity`). The two names denote
the same physical quantity for 1D critical spin chains; the alias
exists purely for readability in spin-chain contexts.

Full derivation: **[XXZ Luttinger parameters from Bethe ansatz
](../../calc/xxz-luttinger-parameters.md)**.

### QAtlas API

```julia
QAtlas.fetch(XXZ1D(; J=1.0, Δ=0.0), LuttingerVelocity(), Infinite())
# → 1.0

QAtlas.fetch(XXZ1D(; J=1.0, Δ=1.0), SpinWaveVelocity(), Infinite())
# → 1.5707963267948966  (= π/2)
```

---

## Legacy Symbol API

Symbol-dispatch calls are still routed through the v0.13 deprecation
layer (`src/deprecate/legacy_xxz.jl`):

```julia
QAtlas.fetch(:XXZ, :energy, Infinite(); J=1.0, Δ=0.0)              # → -1/π
QAtlas.fetch(:XXZ, :spin_wave_velocity, Infinite(); J=1.0, Δ=1.0)  # → π/2
QAtlas.fetch(:XXZ, :luttinger_parameter, Infinite(); J=1.0, Δ=0.0) # → 1.0
```

Recognised quantity aliases for the velocity family include
`:v_F`, `:v_LL`, `:fermi_velocity`, `:spin_wave_velocity`,
`:sound_velocity`, and the capitalised struct names.

---

## Connections

- **Heisenberg limit**: $\Delta = 1$ reproduces the Hulthén result
  cached in the [Heisenberg model page](heisenberg.md). The Heisenberg
  chain is the isotropic SU(2)-symmetric point of the XXZ family.
- **Free-fermion limit**: $\Delta = 0$ reduces to the XX chain, which
  Jordan–Wigner-maps to a 1D nearest-neighbour tight-binding model.
  The value $-1/\pi$ can equivalently be read off the free-fermion
  cosine band — see [JW → TFIM BdG](../../calc/jw-tfim-bdg.md) for the
  analogous mapping on the Ising side.
- **Universality**: The entire $-1 < \Delta < 1$ window sits in the
  $c = 1$ compactified-boson universality class (Luttinger liquid).
  The compactification radius varies continuously with $\Delta$ via
  $K(\Delta)$, which controls every long-distance exponent of the
  chain.

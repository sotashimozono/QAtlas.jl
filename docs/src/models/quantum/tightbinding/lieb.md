# Lieb Lattice

## Overview

The nearest-neighbor tight-binding model on the Lieb lattice
(line-centred square lattice) is the canonical example of a flat band
arising from bipartite sublattice imbalance. The Lieb lattice has
three sublattices per unit cell: A (corner site) and B, C (edge-centre
sites). Only A--B and A--C bonds exist; there is **no direct B--C
bond**. The lattice is therefore bipartite between $\{A\}$ and
$\{B, C\}$.

$$H = -t \sum_{\langle i,j \rangle} \bigl(c^\dagger_i c_j + c^\dagger_j c_i\bigr)$$

**Lattice properties**: 3 sublattices per unit cell (A corner, B right
edge, C top edge), bipartite, coordination number 4 (A) / 2 (B, C),
flat band at $E = 0$.

**Key physics**: Because the B and C sublattices have no mutual
hopping, a localized state on B and C sites can be constructed at
every momentum with zero energy, giving a completely flat band at
$E = 0$. This is a consequence of Lieb's theorem on bipartite
Hubbard models.

---

## Bloch Spectrum

### Statement

The $3 \times 3$ Bloch Hamiltonian in the sublattice basis is

$$H(\mathbf{k}) = -2t \begin{pmatrix} 0 & \cos(\theta_1/2) & \cos(\theta_2/2) \\ \cos(\theta_1/2) & 0 & 0 \\ \cos(\theta_2/2) & 0 & 0 \end{pmatrix}$$

where $\theta_1 = 2\pi m/L_x$ and $\theta_2 = 2\pi n/L_y$. The
zero B--C block yields an exact closed-form spectrum:

$$E_{mn} \in \left\{-E(\mathbf{k}),\; 0,\; +E(\mathbf{k})\right\}$$

$$E(\mathbf{k}) = 2t\,\sqrt{\cos^2\!\left(\frac{\pi m}{L_x}\right) + \cos^2\!\left(\frac{\pi n}{L_y}\right)}$$

The $E = 0$ eigenvalue appears at *every* momentum, forming the
dispersionless flat band. The spectrum is symmetric about zero
(bipartite chiral symmetry).

### Flat band and M-point band touching

For a generic $(L_x, L_y)$, the flat band contributes exactly
$L_x L_y$ zero eigenvalues. However, when both $L_x$ and $L_y$ are
**even**, the M-point $(\theta_1, \theta_2) = (\pi, \pi)$ is an
allowed momentum, and there
$E(\mathbf{k}) = 2t\sqrt{\cos^2(\pi/2) + \cos^2(\pi/2)} = 0$.
Both dispersive bands collapse to zero at this point, contributing
two *extra* zero modes. The total count of $E = 0$ eigenvalues is:

$$N_{E=0} = L_x L_y + \begin{cases} 2 & \text{if both } L_x, L_y \text{ even} \\ 0 & \text{otherwise} \end{cases}$$

### Band structure summary

| Band | Energy range | Notes |
|------|-------------|-------|
| Lower dispersive | $[-2\sqrt{2}\,t,\; 0)$ | Minimum at $\Gamma$: $E = -2t\sqrt{2}$ |
| Flat band | $0$ | $L_x L_y$ states |
| Upper dispersive | $(0,\; +2\sqrt{2}\,t]$ | Maximum at $\Gamma$: $E = +2t\sqrt{2}$ |

### Derivation

See **[Bloch Lieb Flat Band](../../../calc/bloch-lieb-flat-band.md)** for
the derivation of $H(\mathbf{k})$, the closed-form eigenvalues, and
the proof that $E = 0$ is an exact eigenvalue at every momentum.

### References

- E. H. Lieb, "Two Theorems on the Hubbard Model", Phys. Rev. Lett.
  **62**, 1201 (1989) -- Lieb's theorem on flat bands in bipartite
  lattices.
- H. Tasaki, "From Nagaoka's Ferromagnetism to Flat-Band Ferromagnetism
  and Beyond", Prog. Theor. Phys. **99**, 489 (1998) -- review of
  flat-band physics.

### QAtlas API

```julia
# Sorted single-particle spectrum, 3×3 Lieb PBC
λ = QAtlas.fetch(QAtlas.Lieb(), TightBindingSpectrum(); Lx=3, Ly=3, t=1.0)
# → 27 eigenvalues; exactly 9 zero modes (flat band, no M-point bonus)
```

### Verification

| Test file | Method | What is checked |
|-----------|--------|-----------------|
| `test_lieb_tight_binding.jl` | Real-space ED via Lattice2D | $\lambda_{\text{real}} = \lambda_{\text{Bloch}}$ for $2 \times 2$ through $4 \times 4$ |
| `test_lieb_tight_binding.jl` | Bipartite symmetry | $\sum \lambda_i = 0$ and $\{\lambda\} = \{-\lambda\}$ |
| `test_lieb_tight_binding.jl` | Flat band count | $N_{E=0} = L_x L_y$ (odd sizes) or $L_x L_y + 2$ (both even) |
| `test_lieb_tight_binding.jl` | Ground state | $E_{\min} = -2\sqrt{2}\,t$ (at $\Gamma$) |
| `test_lieb_tight_binding.jl` | $t$ scaling | $\lambda(t) = t \cdot \lambda(1)$ |
| `test_lieb_tight_binding.jl` | $3 \times 3$ M-point | No extra zero modes (both dimensions odd) |
| `test_bloch_generic.jl` | Generic Bloch builder | Closed form $=$ `bloch_tb_spectrum` |

---

## Connections

- **Lattice family**: Part of the [tight-binding model catalogue](index.md).
  The Lieb lattice is bipartite like [Honeycomb](honeycomb.md) but has
  a flat band due to its sublattice imbalance, unlike the
  equal-sublattice honeycomb.
- **Flat-band contrast**: Flat band at $E = 0$ (bipartite origin),
  compared to $E = +2t$ for the [Kagome](kagome.md) lattice
  (frustration origin).
- **Dice lattice**: The Dice ($T_3$) lattice shares the bipartite
  flat-band mechanism with the Lieb lattice (flat band at $E = 0$,
  hub vs rim sublattice imbalance).
- **Methods**: Computed via the
  [Bloch Hamiltonian](../../../methods/bloch-hamiltonian/index.md) method.

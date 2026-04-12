# Kagome Lattice

## Overview

The nearest-neighbor tight-binding model on the kagome lattice is the
archetypal example of a flat band arising from geometric frustration.
The kagome lattice has three sublattices (A, B, C) forming
corner-sharing triangles. It is *not* bipartite.

$$H = -t \sum_{\langle i,j \rangle} \bigl(c^\dagger_i c_j + c^\dagger_j c_i\bigr)$$

**Lattice properties**: 3 sublattices per unit cell, non-bipartite
(frustrated), coordination number 4, flat band at $E = +2t$.

**Key physics**: One of the three bands is completely dispersionless
at $E = +2t$. The flat band touches the upper dispersive band at the
$\Gamma$-point ($\mathbf{k} = 0$), where the spectrum degenerates to
$\{-4t, +2t, +2t\}$. The flat band arises from destructive interference
of hopping amplitudes around the corner-sharing triangles.

!!! note "Name clash"
    The dispatch tag `Kagome` is **not exported** from QAtlas to avoid
    a name clash with the `Lattice2D.Kagome` topology type. Use the
    fully qualified form `QAtlas.Kagome()` in user code.

---

## Bloch Hamiltonian

### Statement

In the three-sublattice basis, the Bloch Hamiltonian is the
$3 \times 3$ real symmetric matrix

$$H(\mathbf{k}) = -2t \begin{pmatrix} 0 & \cos(\theta_1/2) & \cos(\theta_2/2) \\ \cos(\theta_1/2) & 0 & \cos\bigl((\theta_2 - \theta_1)/2\bigr) \\ \cos(\theta_2/2) & \cos\bigl((\theta_2 - \theta_1)/2\bigr) & 0 \end{pmatrix}$$

where $\theta_1 = \mathbf{k} \cdot \mathbf{a}_1 = 2\pi m/L_x$ and
$\theta_2 = \mathbf{k} \cdot \mathbf{a}_2 = 2\pi n/L_y$, with
$\mathbf{a}_1, \mathbf{a}_2$ the primitive vectors of `Lattice2D`'s
kagome topology.

Each off-diagonal element collects both the intra-cell and
inter-cell hoppings of a sublattice pair (A--B, A--C, or B--C),
yielding the factor of 2 and the half-angle arguments.

### Flat band

Regardless of $\mathbf{k}$, one eigenvalue of $H(\mathbf{k})$ is
always $+2t$. This can be verified directly: the vector
$(1, -e^{i\theta_1/2}\cos(\theta_2/2)/\cos(\theta_1/2), \ldots)$
lies in the kernel of $H - 2tI$ at every momentum (appropriately
regularized at singular points).

For an $L_x \times L_y$ lattice, the flat band contributes $L_x L_y$
eigenvalues at $+2t$. The upper dispersive band also reaches $+2t$ at
the $\Gamma$-point, contributing one additional eigenvalue there, for a
total degeneracy of $L_x L_y + 1$.

### Band structure summary

| Band | Energy range | Degeneracy at $\Gamma$ |
|------|-------------|------------------------|
| Lower dispersive | $[-4t, +2t)$ | 1 (at $-4t$) |
| Upper dispersive | $(-4t, +2t]$ | 1 (at $+2t$, touching flat band) |
| Flat band | $+2t$ | $L_x L_y$ (all $\mathbf{k}$) |

### Derivation

See **[Bloch Kagome Flat Band](../../calc/bloch-kagome-flat-band.md)**
for the derivation of $H(\mathbf{k})$ and the proof that one
eigenvalue is identically $+2t$.

### References

- I. Syozi, "Statistics of Kagome Lattice", Prog. Theor. Phys. **6**,
  306 (1951) -- original lattice definition.
- D. L. Bergman, C. Wu, L. Balents, "Band touching from real-space
  topology in frustrated hopping models", Phys. Rev. B **78**, 125104
  (2008) -- flat band analysis.

### QAtlas API

```julia
# Sorted single-particle spectrum, 3×3 kagome PBC
λ = QAtlas.fetch(QAtlas.Kagome(), TightBindingSpectrum(); Lx=3, Ly=3, t=1.0)
# → 27 eigenvalues; 10 of them equal to +2.0 (9 flat band + 1 Γ-touch)
```

### Verification

| Test file | Method | What is checked |
|-----------|--------|-----------------|
| `test_kagome_tight_binding.jl` | Real-space ED via Lattice2D | $\lambda_{\text{real}} = \lambda_{\text{Bloch}}$ for $2 \times 2$ through $4 \times 4$ |
| `test_kagome_tight_binding.jl` | Flat band count | Exactly $L_x L_y + 1$ eigenvalues at $+2t$ |
| `test_kagome_tight_binding.jl` | Ground state | $E_{\min} = -4t$ (unique, at $\Gamma$) |
| `test_kagome_tight_binding.jl` | Structural | $\text{tr}\,H = 0$ |
| `test_kagome_tight_binding.jl` | $t$ scaling | $\lambda(t) = t \cdot \lambda(1)$ |
| `test_bloch_generic.jl` | Generic Bloch builder | Hardcoded $3 \times 3$ diagonalization $=$ `bloch_tb_spectrum` |

---

## Connections

- **Lattice family**: Part of the [tight-binding model catalogue](index.md).
  The kagome lattice is frustrated (non-bipartite) like
  [Triangular](triangular.md), but features a flat band that the
  triangular lattice lacks.
- **Flat-band contrast**: Flat band at $E = +2t$, compared to $E = 0$
  for the [Lieb](lieb.md) lattice. The kagome flat band arises from
  frustration, while the Lieb flat band arises from bipartite sublattice
  imbalance.
- **Methods**: Computed via the
  [Bloch Hamiltonian](../../methods/bloch-hamiltonian/index.md) method.

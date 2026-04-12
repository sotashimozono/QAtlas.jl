# Graphene (Honeycomb)

## Overview

The nearest-neighbor tight-binding model on the honeycomb lattice is
the standard model for the electronic band structure of graphene. It
features two sublattices (A, B), bipartite structure, and the
celebrated Dirac cones at the $K$ and $K'$ points of the Brillouin
zone.

$$H = -t \sum_{\langle i,j \rangle \in A\text{-}B} \bigl(c^\dagger_i c_j + c^\dagger_j c_i\bigr)$$

**Lattice properties**: 2 sublattices per unit cell, bipartite
(A vs B), coordination number 3, no flat band.

**Key physics**: The bipartite (chiral) symmetry guarantees that the
spectrum is symmetric about $E = 0$. The two bands touch linearly at
the $K$ and $K'$ points, forming massless Dirac cones with a linear
dispersion $E \sim \pm v_F |\mathbf{k} - \mathbf{K}|$ near the
Fermi level at half filling.

---

## Bloch Spectrum

### Statement

The $2 \times 2$ Bloch Hamiltonian in the sublattice basis is

$$H(\mathbf{k}) = \begin{pmatrix} 0 & f(\mathbf{k}) \\ f(\mathbf{k})^* & 0 \end{pmatrix}$$

with $f(\mathbf{k}) = -t\bigl[e^{i\mathbf{k}\cdot\boldsymbol{\delta}_1} + e^{i\mathbf{k}\cdot\boldsymbol{\delta}_2} + e^{i\mathbf{k}\cdot\boldsymbol{\delta}_3}\bigr]$,
where $\boldsymbol{\delta}_{1,2,3}$ are the three A-to-B nearest-neighbor
displacement vectors.

Using the unit-cell basis $(\mathbf{a}_1, \mathbf{a}_2)$ adopted by
`Lattice2D`'s `Honeycomb` topology and defining $\theta_1 = \mathbf{k} \cdot \mathbf{a}_1 = 2\pi m/L_x$,
$\theta_2 = \mathbf{k} \cdot \mathbf{a}_2 = 2\pi n/L_y$, the two
eigenvalues at each allowed momentum are

$$E_{mn,\pm} = \pm\, t\,\sqrt{3 + 2\cos\!\left(\frac{2\pi m}{L_x}\right) + 2\cos\!\left(\frac{2\pi n}{L_y}\right) + 2\cos\!\left(\frac{2\pi n}{L_y} - \frac{2\pi m}{L_x}\right)}$$

The full spectrum of $2 L_x L_y$ eigenvalues is obtained by
ranging over $m \in \{0, \ldots, L_x - 1\}$ and
$n \in \{0, \ldots, L_y - 1\}$.

### Derivation

See **[Bloch Honeycomb Dispersion](../../calc/bloch-honeycomb-dispersion.md)**
for the full derivation of $|f(\mathbf{k})|^2$ from the unit-cell
geometry.

### Dirac points

When $L_x$ and $L_y$ are both divisible by 3, the $K$ and $K'$ points
$(m/L_x, n/L_y) = (1/3, 2/3)$ and $(2/3, 1/3)$ are commensurate with
the finite Brillouin zone, and $|f(\mathbf{K})|^2 = 0$ exactly. Each
Dirac point contributes two zero-energy modes (one per band), giving 4
zero modes total. For example, at $3 \times 3$ the spectrum contains
exactly 4 zero eigenvalues.

### References

- P. R. Wallace, "The Band Theory of Graphite", Phys. Rev. **71**, 622
  (1947) -- original tight-binding calculation.
- A. H. Castro Neto et al., "The electronic properties of graphene",
  Rev. Mod. Phys. **81**, 109 (2009) -- comprehensive review.

### QAtlas API

```julia
# Sorted single-particle spectrum, 3×3 honeycomb PBC
λ = QAtlas.fetch(Graphene(), TightBindingSpectrum(); Lx=3, Ly=3, t=1.0)
# → 18 eigenvalues, including 4 zero modes at Dirac points
```

### Verification

| Test file | Method | What is checked |
|-----------|--------|-----------------|
| `test_graphene_tight_binding.jl` | Real-space ED via Lattice2D | $\lambda_{\text{real}} = \lambda_{\text{Bloch}}$ for $2 \times 2$ through $4 \times 4$ |
| `test_graphene_tight_binding.jl` | Chiral symmetry | $\sum \lambda_i = 0$ and $\{\lambda\} = \{-\lambda\}$ |
| `test_graphene_tight_binding.jl` | Dirac points ($3 \times 3$) | Exactly 4 zero modes |
| `test_graphene_tight_binding.jl` | $t$ scaling | $\lambda(t) = t \cdot \lambda(1)$ |
| `test_bloch_generic.jl` | Generic Bloch builder | Hardcoded formula $=$ generic `bloch_tb_spectrum` |

---

## Connections

- **Lattice family**: Part of the [tight-binding model catalogue](index.md).
  The honeycomb lattice is bipartite like
  [Lieb](lieb.md) but unlike [Kagome](kagome.md) and
  [Triangular](triangular.md).
- **Flat-band contrast**: Unlike the Kagome ($+2t$) and Lieb ($0$)
  lattices, the honeycomb has no flat band -- both bands are fully
  dispersive.
- **Methods**: Computed via the
  [Bloch Hamiltonian](../../methods/bloch-hamiltonian/index.md)
  method; the generic builder provides an independent cross-check.

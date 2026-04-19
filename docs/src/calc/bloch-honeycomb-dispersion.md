# Bloch Hamiltonian for the Honeycomb Lattice

## Setup

Nearest-neighbour tight-binding model on the honeycomb lattice with
two sublattices $A$ and $B$ and hopping amplitude $t$.

Lattice vectors:

$$\mathbf{a}_1 = (\sqrt{3},\,0), \qquad
  \mathbf{a}_2 = (\sqrt{3}/2,\,3/2)$$

(in units of the bond length). The three nearest-neighbour vectors
from an $A$ site to its three $B$ neighbours:

$$\boldsymbol{\delta}_1 = (0,\,1), \quad
  \boldsymbol{\delta}_2 = (\sqrt{3}/2,\,-1/2), \quad
  \boldsymbol{\delta}_3 = (-\sqrt{3}/2,\,-1/2)$$

## Calculation

### Bloch Hamiltonian

Fourier-transforming the tight-binding Hamiltonian
$H = -t \sum_{\langle i,j \rangle} c_i^\dagger c_j$ in the
two-site unit cell gives the $2 \times 2$ Bloch Hamiltonian:

$$H(\mathbf{k}) = \begin{pmatrix} 0 & f(\mathbf{k}) \\
  f^*(\mathbf{k}) & 0 \end{pmatrix}$$

where the off-diagonal element is

$$f(\mathbf{k}) = -t \sum_{i=1}^{3} e^{i\mathbf{k}\cdot\boldsymbol{\delta}_i}
  = -t\bigl(e^{ik_y} + e^{i(\sqrt{3}k_x/2 - k_y/2)}
    + e^{i(-\sqrt{3}k_x/2 - k_y/2)}\bigr)$$

### Dispersion relation

The eigenvalues are $E(\mathbf{k}) = \pm |f(\mathbf{k})|$. Computing
$|f(\mathbf{k})|^2$:

$$|f(\mathbf{k})|^2 = t^2\bigl[3 + 2\cos(\mathbf{k}\cdot\mathbf{a}_1)
  + 2\cos(\mathbf{k}\cdot\mathbf{a}_2)
  + 2\cos(\mathbf{k}\cdot(\mathbf{a}_2 - \mathbf{a}_1))\bigr]$$

Expanding the dot products with $\theta_1 = \mathbf{k}\cdot\mathbf{a}_1$,
$\theta_2 = \mathbf{k}\cdot\mathbf{a}_2$:

$$E(\mathbf{k}) = \pm t\sqrt{3 + 2\cos\theta_1 + 2\cos\theta_2
  + 2\cos(\theta_2 - \theta_1)}$$

### Finite-size spectrum

On a finite $L_x \times L_y$ lattice with periodic boundary conditions,
the allowed momenta are $\mathbf{k} = (2\pi m/L_x,\, 2\pi n/L_y)$,
giving eigenvalues:

$$E_{m,n} = \pm t\sqrt{3 + 2\cos(2\pi m/L_x)
  + 2\cos(2\pi n/L_y)
  + 2\cos(2\pi(n/L_y - m/L_x))}$$

### Dirac points

The bands touch at the Dirac points $K$ and $K'$ where $|f(\mathbf{k})| = 0$.
These are located at the corners of the Brillouin zone:

$$K = \frac{2\pi}{3}\left(\frac{1}{\sqrt{3}},\, 1\right), \qquad
  K' = \frac{2\pi}{3}\left(\frac{1}{\sqrt{3}},\, -1\right)$$

Near the Dirac points the dispersion is linear:

$$E(\mathbf{K} + \mathbf{q}) \approx \pm \frac{3t}{2}|\mathbf{q}|$$

with Fermi velocity $v_F = 3t/2$ (in units where the lattice constant
is unity). This is the origin of the massless Dirac fermion description
of graphene.

## Result

$$\boxed{E(\mathbf{k}) = \pm t\sqrt{3 + 2\cos(\mathbf{k}\cdot\mathbf{a}_1)
  + 2\cos(\mathbf{k}\cdot\mathbf{a}_2)
  + 2\cos\bigl(\mathbf{k}\cdot(\mathbf{a}_2 - \mathbf{a}_1)\bigr)}}$$

The spectrum has two bands touching at the Dirac points $K$, $K'$.
At half-filling the system is a semimetal with linear dispersion
$E \sim |\mathbf{q}|$ near the Dirac points.

## References

- P. R. Wallace, Phys. Rev. **71**, 622 (1947) — original band structure of graphite.
- A. H. Castro Neto, F. Guinea, N. M. R. Peres, K. S. Novoselov, A. K. Geim, Rev. Mod. Phys. **81**, 109 (2009) — comprehensive review.

## Used by

- [Honeycomb Tight-Binding Model](../models/quantum/tightbinding/honeycomb.md)

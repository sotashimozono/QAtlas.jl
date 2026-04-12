# Lieb Lattice Flat Band at E=0

## Setup

Nearest-neighbour tight-binding model on the Lieb lattice
(line-centred square lattice, also called the CuO$_2$ lattice).
The unit cell contains 3 sublattices:

- $A$: corner sites (e.g. Cu)
- $B$: edge-centre sites along $\mathbf{a}_1$ (e.g. O$_x$)
- $C$: edge-centre sites along $\mathbf{a}_2$ (e.g. O$_y$)

Lattice vectors:

$$\mathbf{a}_1 = (2,\, 0), \qquad \mathbf{a}_2 = (0,\, 2)$$

(in units of the nearest-neighbour distance). Crucially, there is
**no direct $B$-$C$ hopping** --- $B$ and $C$ sites connect only
through $A$ sites. We define $\theta_1 = \mathbf{k}\cdot\mathbf{a}_1$,
$\theta_2 = \mathbf{k}\cdot\mathbf{a}_2$.

## Calculation

### Bloch Hamiltonian

The $3 \times 3$ Bloch Hamiltonian in the basis $(A, B, C)$ is:

$$H(\mathbf{k}) = -2t
\begin{pmatrix}
  0 & \cos(\theta_1/2) & \cos(\theta_2/2) \\
  \cos(\theta_1/2) & 0 & 0 \\
  \cos(\theta_2/2) & 0 & 0
\end{pmatrix}$$

The zeros in the $B$-$C$ block reflect the absence of direct $B$-$C$
bonds. The Hamiltonian is real and symmetric for all $\mathbf{k}$.

### Characteristic polynomial

The secular equation $\det(H - \lambda I) = 0$ gives:

$$-\lambda^3 + 4t^2\bigl[\cos^2(\theta_1/2) + \cos^2(\theta_2/2)\bigr]\lambda = 0$$

Factoring:

$$\lambda\bigl[\lambda^2 - 4t^2\bigl(\cos^2(\theta_1/2) + \cos^2(\theta_2/2)\bigr)\bigr] = 0$$

### The flat band at $E = 0$

One root is $\lambda = 0$ for **all** $\mathbf{k}$. This is the
flat band. Its eigenvector is

$$|\psi_0(\mathbf{k})\rangle = \frac{1}{\mathcal{N}}
  \bigl(0,\; -\cos(\theta_2/2),\; \cos(\theta_1/2)\bigr)^T$$

which has zero weight on the $A$ sublattice --- the flat-band states
live entirely on the $B/C$ sublattice.

### Physical origin

The Lieb lattice is bipartite with sublattice imbalance: the $A$
sublattice has 1 site per unit cell while $B \cup C$ has 2 sites.
By Lieb's theorem (Lieb 1989), the difference in sublattice
counts guarantees a flat band at zero energy, and for the
half-filled Hubbard model, a ground-state spin
$S = |N_A - N_{B \cup C}|/2 \neq 0$ (flat-band ferromagnetism).

### Dispersive bands

The remaining two eigenvalues are:

$$\lambda_\pm(\mathbf{k}) = \pm E(\mathbf{k}), \qquad
  E(\mathbf{k}) = 2t\sqrt{\cos^2(\theta_1/2) + \cos^2(\theta_2/2)}$$

The bandwidth of each dispersive band is from
$E_{\min} = 2t$ (at zone corners) to $E_{\max} = 2t\sqrt{2}$ (at $\Gamma$),
so the dispersive bands touch the flat band only when both
$\cos(\theta_1/2) = 0$ and $\cos(\theta_2/2) = 0$, i.e. at the
$M = (\pi, \pi)$ point where $E(\mathbf{k}) = 0$. This is a Dirac-like
touching point.

## Result

$$\boxed{\text{Spectrum: } \{-E(\mathbf{k}),\; 0,\; +E(\mathbf{k})\}, \qquad
  E(\mathbf{k}) = 2t\sqrt{\cos^2(\theta_1/2) + \cos^2(\theta_2/2)}}$$

The flat band at $E = 0$ exists for all momenta $\mathbf{k}$ and
arises from the bipartite sublattice imbalance ($|A| < |B \cup C|$).

## References

- E. H. Lieb, Phys. Rev. Lett. **62**, 1201 (1989) — sublattice theorem and flat-band ferromagnetism.
- H. Tasaki, Prog. Theor. Phys. **99**, 489 (1998) — rigorous treatment of flat-band magnetism.

## Used by

- [Lieb Tight-Binding Model](../models/quantum/tightbinding/lieb.md)

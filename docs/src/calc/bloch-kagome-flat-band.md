# Kagome Bloch Hamiltonian and the Flat Band

## Setup

Nearest-neighbour tight-binding model on the kagome lattice with
hopping amplitude $t$. The kagome lattice has 3 sublattices per unit
cell, here labelled $A$ (hub), $B$, and $C$ (edge sites).

Lattice vectors:

$$\mathbf{a}_1 = (1,\, 0), \qquad
  \mathbf{a}_2 = (1/2,\, \sqrt{3}/2)$$

(in units of twice the nearest-neighbour distance). We define the
phase variables $\theta_1 = \mathbf{k}\cdot\mathbf{a}_1$ and
$\theta_2 = \mathbf{k}\cdot\mathbf{a}_2$.

## Calculation

### Bloch Hamiltonian

Fourier-transforming the tight-binding Hamiltonian in the three-site
unit cell yields the $3 \times 3$ Bloch matrix:

$$H(\mathbf{k}) = -2t
\begin{pmatrix}
  0 & \cos(\theta_1/2) & \cos(\theta_2/2) \\
  \cos(\theta_1/2) & 0 & \cos((\theta_2 - \theta_1)/2) \\
  \cos(\theta_2/2) & \cos((\theta_2 - \theta_1)/2) & 0
\end{pmatrix}$$

The cosine entries arise from $e^{i\theta/2} + e^{-i\theta/2} = 2\cos(\theta/2)$,
reflecting the two bonds per pair of sublattices within and across
unit cells.

### Characteristic polynomial

Setting $\lambda' = \lambda/(2t)$ and $c_1 = \cos(\theta_1/2)$,
$c_2 = \cos(\theta_2/2)$, $c_3 = \cos((\theta_2 - \theta_1)/2)$,
the secular equation is:

$$-\lambda'^3 + (c_1^2 + c_2^2 + c_3^2)\lambda' + 2c_1 c_2 c_3 = 0$$

Using the identity $c_1^2 + c_2^2 + c_3^2 + 2c_1 c_2 c_3
= 1 + 2\cos(\theta_1/2)\cos(\theta_2/2)\cos((\theta_2 - \theta_1)/2)$
and factoring, one root is always:

$$\lambda'_{\mathrm{flat}} = -1
  \qquad\Longleftrightarrow\qquad
  \lambda_{\mathrm{flat}} = +2t$$

### Proof of the flat band

Substituting $\lambda' = -1$ into the characteristic polynomial:

$$-(-1)^3 + (c_1^2 + c_2^2 + c_3^2)(-1) + 2c_1 c_2 c_3
  = 1 - c_1^2 - c_2^2 - c_3^2 + 2c_1 c_2 c_3$$

This expression vanishes identically by the trigonometric identity:

$$\cos^2\alpha + \cos^2\beta + \cos^2(\alpha + \beta)
  - 2\cos\alpha\cos\beta\cos(\alpha + \beta) = 1$$

with $\alpha = \theta_1/2$, $\beta = (\theta_2 - \theta_1)/2$,
$\alpha + \beta = \theta_2/2$. Therefore $\lambda = +2t$ is an
eigenvalue for all $\mathbf{k}$.

### Full spectrum

Factoring out the flat-band root, the remaining two bands are:

$$\lambda_\pm = -t \pm t\sqrt{4(c_1^2 + c_2^2 + c_3^2) - 3}$$

### $\Gamma$ point

At $\mathbf{k} = (0, 0)$: $c_1 = c_2 = c_3 = 1$, so
$c_1^2 + c_2^2 + c_3^2 = 3$ and:

$$\lambda_+ = -t + t\sqrt{9} = 2t, \qquad
  \lambda_- = -t - t\sqrt{9} = -4t$$

Eigenvalues: $\{-4t,\, 2t,\, 2t\}$. The flat band touches the
upper dispersive band at $\Gamma$, forming a quadratic band touching
point.

## Result

$$\boxed{\text{Eigenvalues: } \lambda_{\mathrm{flat}} = +2t\ (\text{all } \mathbf{k}), \quad
  \lambda_\pm(\mathbf{k}) = -t \pm t\sqrt{4(c_1^2 + c_2^2 + c_3^2) - 3}}$$

The kagome tight-binding model has an exactly flat band at
$E = +2t$ for all momenta. The flat band touches the upper dispersive
band quadratically at the $\Gamma$ point. The macroscopic degeneracy
of the flat band leads to localised eigenstates and is the origin of
flat-band ferromagnetism in the Hubbard model on the kagome lattice.

## References

- I. Syozi, Prog. Theor. Phys. **6**, 306 (1951) — original kagome lattice study.
- D. L. Bergman, C. Wu, L. Balents, Phys. Rev. B **78**, 125104 (2008) — flat band analysis.
- E. H. Lieb, Phys. Rev. Lett. **62**, 1201 (1989) — flat-band ferromagnetism.

## Used by

- [Kagome Tight-Binding Model](../models/quantum/tightbinding/kagome.md)

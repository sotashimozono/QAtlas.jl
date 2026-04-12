# Yang's Spontaneous Magnetization via Toeplitz Determinant

## Setup

The 2D Ising model on an infinite square lattice with nearest-neighbour
coupling $J$ at temperature $T < T_c$. The spontaneous magnetization
$M(T) = \langle \sigma_0 \rangle$ is computed as the limit of a
spin-spin correlation function at infinite separation:

$$M(T) = \lim_{|r|\to\infty} \langle \sigma_0 \sigma_r \rangle^{1/2}$$

## Calculation

### Correlation as a Toeplitz determinant

The correlation function $\langle \sigma_{0,0}\,\sigma_{n,0} \rangle$
along a lattice row can be written as the $n \times n$ Toeplitz
determinant

$$\langle \sigma_{0,0}\,\sigma_{n,0} \rangle = D_n(a)
  = \det\bigl[a_{j-k}\bigr]_{j,k=0}^{n-1}$$

where the symbol $a(\theta)$ is determined by the transfer-matrix
structure (Kaufman 1949; Montroll, Potts, Ward 1963):

$$a(\theta) = \frac{(1 - \alpha_1 e^{i\theta})(1 - \alpha_2 e^{-i\theta})}
  {|1 - \alpha_1 e^{i\theta}||1 - \alpha_2 e^{-i\theta}|}$$

with $\alpha_1 = e^{-2K_1^*}$, $\alpha_2 = e^{-2K_2}$, and
$K_i = \beta J_i$, $K_i^*$ the dual couplings.

### Application of Szego's theorem

The strong Szego limit theorem states that for a Toeplitz determinant
with symbol $a(\theta) = \exp[\sum_{k} (c_k e^{ik\theta} + c_{-k}e^{-ik\theta})]$,

$$\lim_{n\to\infty} D_n(a) = \exp\left[\sum_{k=1}^{\infty} k\, c_k c_{-k}\right] \cdot G^n$$

where $G = \exp[c_0]$ is the geometric mean of the symbol. For
$T < T_c$ the geometric mean is nonzero and yields a finite
magnetization in the limit $n\to\infty$.

### Evaluation for the isotropic case

For the isotropic lattice ($J_1 = J_2 = J$, so $K_1 = K_2 = K$),
Yang evaluated the Toeplitz determinant in the $n\to\infty$ limit
to obtain:

$$M(T)^2 = \lim_{n\to\infty} D_n(a)
  = \bigl(1 - \sinh^{-4}(2\beta J)\bigr)^{1/4}$$

Taking the square root:

$$M(T) = \bigl(1 - \sinh^{-4}(2\beta J)\bigr)^{1/8}$$

## Result

$$\boxed{M(T) = \left(1 - \frac{1}{\sinh^4(2\beta J)}\right)^{1/8}},
  \qquad T < T_c$$

The exponent $1/8$ in $M \sim (T_c - T)^{\beta}$ identifies the
order-parameter critical exponent

$$\beta = \tfrac{1}{8}$$

At $T = T_c$ (where $\sinh(2\beta_c J) = 1$) the magnetization
vanishes continuously, confirming the second-order nature of the
transition.

## References

- C. N. Yang, Phys. Rev. **85**, 808 (1952).
- B. Kaufman, Phys. Rev. **76**, 1232 (1949) — transfer-matrix diagonalisation.
- E. W. Montroll, R. B. Potts, J. C. Ward, J. Math. Phys. **4**, 308 (1963) — Toeplitz formulation.
- G. Szego, Math. Z. **6**, 167 (1920) — Toeplitz limit theorem.

## Used by

- [Ising Square: Spontaneous Magnetization](../models/classical/ising-square.md)
- [Ising Universality Class](../universalities/ising.md)

# Transfer Matrix: Symmetric Split Construction

## Setup

For the 2D Ising model on an $L_x \times L_y$ square lattice with PBC,
the partition function is $Z = \mathrm{Tr}(T^{L_x})$ where $T$ is the
transfer matrix acting on a row of $L_y$ spins.

**Goal**: Construct $T$ such that it is **symmetric** (for numerical
stability and AD compatibility).

## Calculation

### Bond decomposition

The total energy decomposes into horizontal bonds (within a row) and
vertical bonds (between adjacent rows):

$$H = -J\sum_{\text{horiz}} \sigma_i\sigma_j - J\sum_{\text{vert}} \sigma_i\sigma_j$$

### Asymmetric transfer matrix

A naive transfer matrix assigns all horizontal bonds to one row:

$$\tilde{T}_{\sigma,\sigma'} = e^{\beta J E_v(\sigma,\sigma')} \cdot e^{\beta J E_h(\sigma')}$$

This is **not symmetric**: $\tilde{T}_{\sigma,\sigma'} \neq \tilde{T}_{\sigma',\sigma}$.

### Symmetric split

Split the horizontal bonds equally between the two rows:

$$T_{\sigma,\sigma'} = e^{\beta J E_h(\sigma)/2} \cdot e^{\beta J E_v(\sigma,\sigma')} \cdot e^{\beta J E_h(\sigma')/2}$$

where:

$$E_h(\sigma) = \sum_{j=1}^{L_y} \sigma_j\sigma_{(j\bmod L_y)+1} \quad\text{(PBC in } y\text{)}$$

$$E_v(\sigma,\sigma') = \sum_{j=1}^{L_y} \sigma_j\sigma'_j$$

**Symmetry proof**: Since $E_v(\sigma,\sigma') = E_v(\sigma',\sigma)$,

$$T_{\sigma,\sigma'} = e^{\beta J E_h(\sigma)/2} \cdot e^{\beta J E_v(\sigma,\sigma')} \cdot e^{\beta J E_h(\sigma')/2} = T_{\sigma',\sigma}$$

### Why `tr(T^Lx)` instead of eigendecomposition

QAtlas computes $Z = \mathrm{tr}(T^{L_x})$ via matrix exponentiation
rather than $Z = \sum_i \lambda_i^{L_x}$ from eigenvalues.
This is because `eigvals(Symmetric(T))` uses LAPACK, which does not
support `ForwardDiff.Dual` inputs. Matrix multiplication and `tr` do
support duals, enabling [automatic differentiation](../methods/automatic-differentiation/index.md) of thermodynamic quantities.

## Result

$T$ is a $2^{L_y} \times 2^{L_y}$ **real symmetric** matrix.
$Z = \mathrm{tr}(T^{L_x})$ is the exact partition function.

## References

- R. J. Baxter, *Exactly Solved Models in Statistical Mechanics*
  (1982), Ch. 7.

## Used by

- [IsingSquare: Partition Function](../models/classical/ising-square.md#partition-function-transfer-matrix)
- [Transfer Matrix method](../methods/transfer-matrix/index.md)

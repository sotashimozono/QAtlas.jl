# Transfer Matrix Method

## What It Is

The transfer matrix method converts the partition function of a
finite 2D classical lattice model into a trace of a matrix power:

$$Z = \mathrm{tr}(T^{L_x})$$

where $T$ is the **transfer matrix** acting on a row of $L_y$ spins
(or more generally, on the degrees of freedom of a single row).
Each matrix element $T_{\sigma, \sigma'}$ encodes the Boltzmann
weight for two adjacent rows in configurations $\sigma$ and
$\sigma'$.

---

## When to Use It

| Condition | Satisfied? | Consequence |
|-----------|------------|-------------|
| 2D classical lattice model | Required | Transfer matrix is a row-to-row operator |
| Short-range interactions | Required | $T$ couples only adjacent rows |
| Finite $L_y$ | Required | $T$ is a $2^{L_y} \times 2^{L_y}$ matrix |
| PBC in $y$-direction | Typical | Ensures translational invariance within a row |

**Limitations**: The transfer matrix dimension grows **exponentially**
in $L_y$. For the Ising model, $T$ is $2^{L_y} \times 2^{L_y}$, so
$L_y \leq 20$ is the practical limit for dense storage.

---

## Construction: Symmetric Split

A naive transfer matrix that assigns all horizontal bonds to one row
is **not symmetric**: $\tilde{T}_{\sigma,\sigma'} \neq \tilde{T}_{\sigma',\sigma}$.
QAtlas uses the symmetric split construction, which distributes the
horizontal bond energy equally between the two rows:

$$T_{\sigma,\sigma'} = e^{\beta J E_h(\sigma)/2} \cdot e^{\beta J E_v(\sigma,\sigma')} \cdot e^{\beta J E_h(\sigma')/2}$$

This ensures $T = T^\top$, which is important for numerical
stability and AD compatibility.

See [calc/transfer-matrix-symmetric-split.md](../../calc/transfer-matrix-symmetric-split.md)
for the full derivation and symmetry proof.

---

## Why `tr(T^Lx)`, Not Eigenvalues

In principle, the partition function can be computed more efficiently
via eigenvalue decomposition:

$$Z = \sum_{i} \lambda_i^{L_x}$$

However, LAPACK's eigenvalue routines (`DSYEV`, etc.) do **not**
support `ForwardDiff.Dual` inputs. The internal Fortran code performs
comparisons and branches incompatible with dual-number arithmetic.

In contrast, matrix multiplication and the trace operation are purely
algebraic and propagate dual numbers correctly. Therefore QAtlas
computes $Z = \mathrm{tr}(T^{L_x})$ via repeated squaring
($O(\log L_x)$ matrix multiplications), enabling
[automatic differentiation](../automatic-differentiation/index.md)
of all thermodynamic quantities.

---

## Thermodynamic Quantities

Once $Z(\beta, J)$ is available as a differentiable function, all
equilibrium observables follow from derivatives of $\ln Z$:

$$\langle E \rangle = -\frac{\partial \ln Z}{\partial \beta}, \qquad
C_v = \beta^2 \frac{\partial^2 \ln Z}{\partial \beta^2}$$

See [calc/ad-thermodynamics-from-z.md](../../calc/ad-thermodynamics-from-z.md)
for the full list of thermodynamic identities.

---

## Applications in QAtlas

| Model | $T$ size | Result |
|-------|----------|--------|
| IsingSquare ($L_y = 4$) | $16 \times 16$ | $Z(\beta, J)$, $\langle E \rangle$, $C_v$, $\langle\sigma\sigma\rangle_{\mathrm{nn}}$ |

---

## References

- R. J. Baxter, *Exactly Solved Models in Statistical Mechanics*
  (Academic Press, 1982), Ch. 7 --- transfer matrix formalism.
- H. A. Kramers, G. H. Wannier, Phys. Rev. **60**, 252 (1941) ---
  transfer matrix for the Ising model.

# Jordan-Wigner Transformation of the TFIM → BdG Matrix

## Setup

Starting Hamiltonian (1D TFIM, OBC, $N$ sites):

$$H = -J\sum_{i=1}^{N-1}\sigma^z_i\sigma^z_{i+1} - h\sum_{i=1}^{N}\sigma^x_i$$

**Goal**: Reduce $H$ to a quadratic (free) fermion form and extract
the BdG quasiparticle energies $\{\Lambda_n\}$.

## Step 1: Kramers-Wannier Duality

The $\sigma^z\sigma^z$ convention does not directly yield a free-fermion
Hamiltonian under Jordan-Wigner (it produces quartic terms). We first
apply the Kramers-Wannier duality to obtain the dual representation.

Define dual spin operators $\tau$ on the **bonds** of the original chain:

$$\tau^x_i = \sigma^z_i\sigma^z_{i+1}, \qquad \tau^z_i\tau^z_{i+1} = \sigma^x_{i+1}$$

The Hamiltonian becomes:

$$H = -J\sum_i \tau^x_i - h\sum_i \tau^z_i\tau^z_{i+1}$$

Now $\tau^z\tau^z$ is the interaction term and $\tau^x$ is the field
— this is the form amenable to Jordan-Wigner.

## Step 2: Jordan-Wigner on the Dual Spins

Apply the JW transformation (see [general JW](../methods/jordan-wigner/index.md)):

$$d_i = \left(\prod_{j<i}\tau^z_j\right)\tau^-_i$$

Key results:

$$\tau^x_i = 1 - 2d^\dagger_i d_i$$

$$\tau^z_i\tau^z_{i+1} = -(d^\dagger_i - d_i)(d^\dagger_{i+1} + d_{i+1})$$

Substituting:

$$H = -J\sum_i(1 - 2d^\dagger_i d_i) + h\sum_i(d^\dagger_i - d_i)(d^\dagger_{i+1} + d_{i+1})$$

$$= 2J\sum_i d^\dagger_i d_i - h\sum_i(d^\dagger_i d_{i+1} + d^\dagger_i d^\dagger_{i+1} + \text{h.c.}) + \text{const}$$

## Step 3: BdG Matrix

Collecting the quadratic terms in the Nambu basis
$(d^\dagger_1, \ldots, d^\dagger_N, d_1, \ldots, d_N)$:

$$H = \frac{1}{2}\begin{pmatrix}d^\dagger & d\end{pmatrix}\begin{pmatrix}A & B \\ -B & -A\end{pmatrix}\begin{pmatrix}d \\ d^\dagger\end{pmatrix} + \text{const}$$

where:

$$A_{ii} = 2h, \quad A_{i,i+1} = A_{i+1,i} = -J$$

$$B_{i,i+1} = J, \quad B_{i+1,i} = -J$$

**Key property**: $H_\text{BdG} = \begin{pmatrix}A & B \\ -B & -A\end{pmatrix}$
is **real symmetric**. Proof: $A$ is symmetric and $B$ is antisymmetric,
so $H_\text{BdG}^T = \begin{pmatrix}A^T & -B^T \\ B^T & -A^T\end{pmatrix} = \begin{pmatrix}A & B \\ -B & -A\end{pmatrix} = H_\text{BdG}$.

## Result

The eigenvalues of $H_\text{BdG}$ come in $\pm\Lambda_n$ pairs.
The positive eigenvalues $\Lambda_n > 0$ are the quasiparticle energies.

$$E_0 = -\frac{1}{2}\sum_n \Lambda_n, \qquad \langle H\rangle(\beta) = -\sum_n \frac{\Lambda_n}{2}\tanh\!\left(\frac{\beta\Lambda_n}{2}\right)$$

## References

- P. Pfeuty, Ann. Phys. **57**, 79 (1970) — TFIM exact solution.
- H. A. Kramers, G. H. Wannier, Phys. Rev. **60**, 252 (1941) — duality.

## Used by

- [TFIM](../models/quantum/tfim.md) — ground-state energy, thermal observables
- [Jordan-Wigner method](../methods/jordan-wigner/index.md) — as application example

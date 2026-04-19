# Transfer Matrix: Symmetric-Split Construction

## Main result

For the 2D classical Ising model on an $L_x \times L_y$ square
lattice with periodic boundary conditions in both directions,

$$H(\{s\}) = -J\,\sum_{\langle u, v\rangle} s_u s_v,
\qquad s_v \in \{\pm 1\},$$

the partition function $Z = \sum_{\{s\}} e^{-\beta H}$ admits the
row-by-row **transfer-matrix representation**

$$\boxed{\;
Z \;=\; \mathrm{Tr}\bigl(T^{L_x}\bigr),
\qquad
T \in \mathbb{R}^{2^{L_y}\times 2^{L_y}},
\;}$$

where $T$ is the $2^{L_y}\times 2^{L_y}$ transfer matrix whose
$(\sigma, \sigma')$ entry is the Boltzmann weight for adding a new
row with spin configuration $\sigma'$ adjacent to a row with
configuration $\sigma$.

**Symmetric-split construction.** Splitting the in-row horizontal
bond energy equally between the two rows connected by the vertical
transfer step gives

$$\boxed{\;
T_{\sigma,\sigma'} \;=\; \exp\!\Bigl[\,
    \tfrac{K}{2}\,E_h(\sigma)
  \;+\; K\,E_v(\sigma, \sigma')
  \;+\; \tfrac{K}{2}\,E_h(\sigma')\,\Bigr],
\;}$$

where $K = \beta J$,
$E_h(\sigma) = \sum_{j=1}^{L_y} \sigma_j\sigma_{j+1}$ is the
horizontal (in-row) bond energy with PBC $\sigma_{L_y + 1} \equiv
\sigma_1$, and
$E_v(\sigma, \sigma') = \sum_{j=1}^{L_y} \sigma_j \sigma'_j$ is the
vertical (inter-row) bond energy. Since
$E_h(\sigma) = E_h(\sigma)$ and $E_v(\sigma, \sigma') =
E_v(\sigma', \sigma)$ by inspection, the matrix is

$$\boxed{\;T_{\sigma,\sigma'} = T_{\sigma', \sigma}\;}$$

— **real symmetric**. This is the form QAtlas uses at
[`src/models/classical/IsingSquare`](../models/classical/ising-square.md);
the symmetry enables automatic-differentiation-compatible evaluation
of thermodynamic quantities via $Z = \mathrm{Tr}(T^{L_x})$
without an LAPACK eigendecomposition (see Step 6 below).

---

## Setup

### 2D classical Ising model

Spins $s_v \in \{\pm 1\}$ on the $L_x \times L_y$ square lattice
(sites indexed by $(i, j)$ with $i \in \{0, 1, \dots, L_x - 1\}$,
$j \in \{1, 2, \dots, L_y\}$), nearest-neighbour coupling $J > 0$,
periodic boundary conditions in both directions. Reduced coupling
$K \equiv \beta J$.

### Row-labelled configurations

Let $\sigma^{(i)} \in \{\pm 1\}^{L_y}$ denote the spin
configuration of row $i$; its $j$-th entry $\sigma^{(i)}_j = s_{i,
j}$. A full lattice configuration is the tuple
$(\sigma^{(0)}, \sigma^{(1)}, \dots, \sigma^{(L_x - 1)})$ with PBC
$\sigma^{(L_x)} \equiv \sigma^{(0)}$.

### Bond decomposition

Every bond is either **horizontal** (within a single row, $i$
fixed) or **vertical** (between two adjacent rows, $j$ fixed):

$$E(\{s\}) \;=\; -\,J\sum_{\text{bonds}}s_{u}s_{v}
 \;=\; -\,J\,\sum_{i = 0}^{L_x - 1}E_h(\sigma^{(i)})
       \;-\; J\,\sum_{i = 0}^{L_x - 1} E_v(\sigma^{(i)}, \sigma^{(i + 1)}),
\tag{1}$$

with the row-energy functions

$$E_h(\sigma) \;=\; \sum_{j=1}^{L_y}\sigma_{j}\sigma_{j + 1}
\qquad(\sigma_{L_y + 1}\equiv \sigma_{1},\text{ PBC in }y),$$

$$E_v(\sigma, \sigma') \;=\; \sum_{j=1}^{L_y}\sigma_{j}\sigma'_{j}.$$

There are $L_x \cdot L_y$ horizontal bonds and $L_x \cdot L_y$
vertical bonds, totaling $2 L_x L_y$ bonds — consistent with a
4-neighbour coordination (4 NN per site × $L_x L_y$ sites ÷ 2 for
double counting).

### Goal

Package (1) into the form $Z = \mathrm{Tr}(T^{L_x})$ with a
real symmetric matrix $T$, so that the partition function can be
computed by $L_x$ successive matrix multiplications.

---

## Calculation

### Step 1 — Row-by-row factorisation of $Z$

Substitute (1) into $Z = \sum_{\{s\}} e^{-\beta E}$ and factorise
the exponential:

$$Z \;=\; \sum_{\sigma^{(0)}, \dots, \sigma^{(L_x - 1)}}
   \exp\!\Bigl[\,K\sum_{i = 0}^{L_x - 1}E_h(\sigma^{(i)})
               + K\sum_{i = 0}^{L_x - 1}
                 E_v(\sigma^{(i)}, \sigma^{(i + 1)})\,\Bigr].$$

Group the exponential by rows. Using the identity
$\sum_i f(\sigma^{(i)}) = \tfrac{1}{2}\sum_i[f(\sigma^{(i)}) +
f(\sigma^{(i)})]$ and then, exploiting the PBC
$\sigma^{(L_x)} \equiv \sigma^{(0)}$,

$$\sum_{i = 0}^{L_x - 1}E_h(\sigma^{(i)})
 = \sum_{i = 0}^{L_x - 1}\tfrac{1}{2}\bigl[E_h(\sigma^{(i)}) + E_h(\sigma^{(i + 1)})\bigr],$$

where the right-hand side has been re-indexed using
$\sigma^{(L_x)} \equiv \sigma^{(0)}$ to pair up the "leftover" term.
This little trick — splitting each horizontal-energy sum into two
halves and reassigning one half to the adjacent row — is what makes
the symmetric transfer matrix possible.

Combining with the vertical energy,

$$Z \;=\; \sum_{\{\sigma^{(i)}\}}
    \prod_{i = 0}^{L_x - 1}
    \exp\!\Bigl[\,\tfrac{K}{2}E_h(\sigma^{(i)})
            + K\,E_v(\sigma^{(i)}, \sigma^{(i + 1)})
            + \tfrac{K}{2}E_h(\sigma^{(i + 1)})\,\Bigr].
\tag{2}$$

Each factor in the product depends only on the pair $(\sigma^{(i)},
\sigma^{(i + 1)})$ of adjacent-row configurations, so (2) matches
the template of a transfer-matrix trace.

### Step 2 — Define the transfer matrix $T$

Read off from (2):

$$\boxed{\;
T_{\sigma, \sigma'}
 \;=\; \exp\!\Bigl[\,
    \tfrac{K}{2}\,E_h(\sigma)
  \;+\; K\,E_v(\sigma, \sigma')
  \;+\; \tfrac{K}{2}\,E_h(\sigma')\,\Bigr].
\;}
\tag{3}$$

This is a $2^{L_y}\times 2^{L_y}$ matrix indexed by the $2^{L_y}$
spin configurations $\sigma \in \{\pm 1\}^{L_y}$. The entry
$T_{\sigma, \sigma'}$ is the Boltzmann weight of placing row
$\sigma$ adjacent to row $\sigma'$ in the "half-half" bond
assignment.

### Step 3 — $Z = \mathrm{Tr}(T^{L_x})$

Substitute (3) back into (2):

$$Z \;=\; \sum_{\sigma^{(0)}, \dots, \sigma^{(L_x - 1)}}
         \prod_{i = 0}^{L_x - 1} T_{\sigma^{(i)}, \sigma^{(i + 1)}}.$$

Using the PBC identification $\sigma^{(L_x)} \equiv \sigma^{(0)}$,
the sum over $\sigma^{(1)}, \dots, \sigma^{(L_x - 1)}$ is a product
of matrices contracted cyclically:

$$\sum_{\sigma^{(1)}} T_{\sigma^{(0)}, \sigma^{(1)}}\,
 \sum_{\sigma^{(2)}} T_{\sigma^{(1)}, \sigma^{(2)}}\dotsb
 T_{\sigma^{(L_x - 1)}, \sigma^{(0)}}
 = (T^{L_x})_{\sigma^{(0)}, \sigma^{(0)}}.$$

Summing over $\sigma^{(0)}$ gives the diagonal trace:

$$\boxed{\;
Z \;=\; \sum_{\sigma^{(0)}} (T^{L_x})_{\sigma^{(0)}, \sigma^{(0)}}
 \;=\; \mathrm{Tr}\bigl(T^{L_x}\bigr).
\;}
\tag{4}$$

### Step 4 — Symmetry of $T$

From (3):

$$T_{\sigma',\sigma}
 = \exp\!\Bigl[\tfrac{K}{2}E_h(\sigma') + K E_v(\sigma', \sigma)
               + \tfrac{K}{2}E_h(\sigma)\Bigr].$$

Compare to $T_{\sigma,\sigma'}$ using the two obvious symmetries

$$E_h(\sigma') = E_h(\sigma')\quad(\text{trivial}),\qquad
E_v(\sigma', \sigma) = E_v(\sigma, \sigma'),$$

the second following from $\sum_j \sigma_j' \sigma_j = \sum_j
\sigma_j \sigma_j'$. Hence

$$T_{\sigma',\sigma}
 = \exp\!\Bigl[\tfrac{K}{2}E_h(\sigma) + K E_v(\sigma, \sigma')
               + \tfrac{K}{2}E_h(\sigma')\Bigr]
 = T_{\sigma,\sigma'}. \tag{5}$$

$T$ is **real symmetric**. Every entry is positive (exponential of
a real number), so $T$ is also strictly positive and the
Perron–Frobenius theorem guarantees a unique largest eigenvalue
$\lambda_{\max} > 0$.

### Step 5 — Contrast with the asymmetric transfer matrix

A naive row-by-row construction assigns **all** horizontal bonds to
one row instead of half-and-half:

$$\widetilde{T}_{\sigma,\sigma'}
 = \exp\!\Bigl[\,K\,E_h(\sigma')\;+\;K\,E_v(\sigma, \sigma')\,\Bigr].
\tag{6}$$

This matrix has the same trace $\mathrm{Tr}(\widetilde{T}^{L_x}) = Z$
(provable by the same cyclic-trace argument, since
$\prod_i \exp[K E_h(\sigma^{(i+1)})] = \exp[K \sum_i E_h(\sigma^{(i+1)})]
= \exp[K \sum_i E_h(\sigma^{(i)})]$ by PBC reindexing). But (6) is
**not symmetric**:

$$\widetilde{T}_{\sigma',\sigma}
 = \exp\!\Bigl[\,K E_h(\sigma)
               + K E_v(\sigma', \sigma)\,\Bigr]
 = \exp\!\Bigl[\,K E_h(\sigma)
               + K E_v(\sigma, \sigma')\,\Bigr]
 \neq \widetilde{T}_{\sigma,\sigma'}$$

unless $E_h(\sigma) = E_h(\sigma')$, which is not generic.

Symmetric (3) and asymmetric (6) are **similar matrices**: one can
check $T = D^{1/2}\widetilde{T}\,D^{-1/2}$ with $D =
\mathrm{diag}_\sigma \exp[K E_h(\sigma)]$, so they share the
same eigenvalues. The traces of all their powers agree
($\mathrm{Tr}\widetilde{T}^{L_x} = \mathrm{Tr}(D^{1/2}
\widetilde{T} D^{-1/2})^{L_x} = \mathrm{Tr} T^{L_x}$). Both
give the same partition function, but the symmetric form has three
distinct advantages (Step 6).

### Step 6 — Why symmetric: eigenvalues, AD, and numerical stability

**(i) Spectral theorem.** Because $T$ is real symmetric, the
spectral theorem guarantees a complete orthonormal eigenbasis with
real eigenvalues $\lambda_1 \ge \lambda_2 \ge \dots \ge
\lambda_{2^{L_y}}$. The partition function can thus be written as

$$Z \;=\; \sum_{i = 1}^{2^{L_y}} \lambda_i^{L_x}
 \;\xrightarrow[L_x\to\infty]{}\; \lambda_{\max}^{L_x},$$

giving immediate access to the **thermodynamic-limit free energy**
$f = -\tfrac{1}{\beta}\lim_{L_x\to\infty}\tfrac{1}{L_x L_y}\ln Z =
-\tfrac{1}{\beta L_y}\ln\lambda_{\max}$.

**(ii) Automatic-differentiation compatibility.** Computing
$\lambda_{\max}$ via `eigvals(Symmetric(T))` routes through the
LAPACK Hermitian-eigenvalue solver (SYEVR / SYEVD). LAPACK operates
on raw `Float64` arrays and does **not** support Julia's
`ForwardDiff.Dual` number type. By contrast,
`Z = tr(T^Lx)` uses only `*` (matrix multiplication) and `tr`,
both of which dispatch over any real-number ring — including dual
numbers — without LAPACK. QAtlas exploits this to compute

$$F(\beta) = -\tfrac{1}{\beta}\ln Z(\beta),\qquad
S(\beta) = -\partial_T F,\qquad C_v(\beta) = -T\partial_T^{2} F$$

as first- and second-order forward-mode derivatives of
$\mathrm{Tr}(T^{L_x})$ with respect to $\beta$ — see
[`ad-thermodynamics-from-z`](ad-thermodynamics-from-z.md) for the
derivation of the chain-rule formulas.

**(iii) Numerical stability.** The asymmetric $\widetilde{T}$ has
condition number $\kappa(\widetilde{T}) = \|\widetilde{T}\|\cdot
\|\widetilde{T}^{-1}\| \sim e^{2K L_y}$ in the low-temperature
limit $K \to \infty$, because $E_h(\sigma)$ ranges over $[-L_y,
L_y]$ and the prefactor $e^{K E_h(\sigma')}$ varies by $e^{2K
L_y}$ across $\sigma'$. The symmetric split halves this range
(the $\tfrac{K}{2}$ factors on *both* ends are bounded by $e^{K
L_y/2}$), giving $\kappa(T) \sim e^{K L_y}$. The better condition
number matters for $L_y \gtrsim 20$ at low temperatures.

### Step 7 — Limiting-case checks

**(i) $\beta = 0$ (infinite temperature).** All bonds contribute
$\exp(0) = 1$, so $T_{\sigma, \sigma'} = 1$ for every pair. The
matrix is the $2^{L_y}\times 2^{L_y}$ all-ones matrix;
$\mathrm{rank}(T) = 1$ with single non-zero eigenvalue
$\lambda_{\max} = 2^{L_y}$ and eigenvector the uniform vector
$(1, 1, \dots, 1)^T$. Then

$$Z = \mathrm{Tr}(T^{L_x}) = \lambda_{\max}^{L_x} = 2^{L_y L_x},$$

the correct infinite-temperature result $Z = 2^{\#\text{spins}}$.

**(ii) $\beta \to \infty$ (zero temperature).** The dominant
configurations are uniform $\sigma = (+, +, \dots, +)$ or
$(-, -, \dots, -)$; each has $E_h = L_y$ and $E_v = L_y$, giving
$T_{\sigma, \sigma} = e^{K L_y}\cdot e^{K L_y} = e^{2 K L_y}$. For
the mixed pair $\sigma = -\sigma'$, $E_v = -L_y$ (anti-aligned
rows), so $T_{\sigma, -\sigma} = e^{K L_y - K L_y + K L_y} = e^{K
L_y}$ — smaller by $e^{K L_y}$. As $K \to \infty$, the
$(\pm +\dots +)$-diagonal entries dominate and $\lambda_{\max}
\to e^{2 K L_y}$, giving $Z \to 2\,e^{2 K L_x L_y}$ in the
ordered-phase free-energy limit.

**(iii) Small-system check $L_y = 1$, $L_x = 2$.** With $L_y = 1$
there is no horizontal bond structure: $E_h(\pm) = 0$ (PBC on a
1-site row gives the bond $\sigma_1\sigma_1 = 1$ — wait, $L_y = 1$
with PBC is a self-loop, which is degenerate. Take $L_y = 2$
instead.) With $L_y = 2$, $\sigma = (s_1, s_2)$ and $E_h(\sigma) =
s_1 s_2 + s_2 s_1 = 2 s_1 s_2$ (PBC gives two identical bonds on
a 2-site ring; this is the standard double-counting artefact of
small-$L_y$ PBC lattices noted in the repository
`CONTRIBUTING.md`).  For $L_x = 2$ the $4\times 4$ matrix $T$ has
$\mathrm{Tr}(T^2) = Z$. The explicit brute-force sum over
the $2^{4} = 16$ spin configurations on the $2\times 2$ torus
matches $\mathrm{Tr}(T^2)$ to machine precision in every
QAtlas test (`test/verification/test_ising_2x2_classical.jl`).

---

## References

- E. W. Montroll, *Statistical mechanics of nearest neighbor
  systems*, J. Chem. Phys. **9**, 706 (1941). One of the earliest
  row-transfer-matrix formulations.
- L. Onsager, Phys. Rev. **65**, 117 (1944). The transfer-matrix
  form used to solve the 2D Ising model.
- R. J. Baxter, *Exactly Solved Models in Statistical Mechanics*
  (Academic Press, 1982), Ch. 7. Standard reference for row and
  symmetric transfer matrices; the symmetric-split formulation is
  eq. (7.2.5) there.
- N. G. van Kampen, *Stochastic Processes in Physics and
  Chemistry*, 3rd ed. (Elsevier, 2007), Ch. III.3. Concise
  treatment of the symmetric transfer matrix for Markov chains;
  same formal structure as the Ising case.

## Used by

- [IsingSquare model page](../models/classical/ising-square.md) —
  `fetch(IsingSquare(; J, Lx, Ly), PartitionFunction(); β)` uses
  (3) and computes (4) by matrix multiplication; the resulting $Z$
  is differentiable in $\beta$ via ForwardDiff.
- [Transfer-matrix method](../methods/transfer-matrix/index.md) —
  general framework; this note is the canonical worked example.
- [AD thermodynamics note](ad-thermodynamics-from-z.md) — uses (4)
  to derive $F, S, C_v$ as higher-order forward-mode derivatives.
- [Kramers–Wannier duality note](kramers-wannier-duality.md) — the
  classical duality is most naturally expressed in the transfer-
  matrix formalism used here.

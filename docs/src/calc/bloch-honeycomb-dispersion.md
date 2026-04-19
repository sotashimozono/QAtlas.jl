# Bloch Hamiltonian for the Honeycomb Lattice

## Main result

For the nearest-neighbour tight-binding model on the honeycomb
lattice with hopping amplitude $t$ (two sublattices $A, B$, three
NN bonds per site), the two-band dispersion in the Brillouin zone is

$$\boxed{\;
E_{\pm}(\mathbf{k}) \;=\; \pm\,t\sqrt{\,3
 \;+\; 2\cos\!\bigl(\mathbf{k}\!\cdot\!\mathbf{a}_1\bigr)
 \;+\; 2\cos\!\bigl(\mathbf{k}\!\cdot\!\mathbf{a}_2\bigr)
 \;+\; 2\cos\!\bigl(\mathbf{k}\!\cdot(\mathbf{a}_2 - \mathbf{a}_1)\bigr)\,},
\;}$$

with primitive lattice vectors $\mathbf{a}_1 = (\sqrt{3}, 0)$,
$\mathbf{a}_2 = (\sqrt{3}/2, 3/2)$ (in units of the NN bond length).
The two bands touch at the **Dirac points** $K, K'$ at the corners
of the hexagonal Brillouin zone, with linear dispersion and Fermi
velocity

$$\boxed{\;
E(\mathbf{K} + \mathbf{q}) \;=\; \pm\,v_F\,|\mathbf{q}|
\;+\; O(|\mathbf{q}|^2),
\qquad
v_F \;=\; \tfrac{3 t}{2}.
\;}$$

The spectrum is symmetric about $E = 0$ (chiral symmetry from the
bipartite $A \leftrightarrow B$ structure); at half-filling the
system is a gapless semimetal — the prototypical massless-Dirac-
fermion realisation of graphene (Wallace 1947).

---

## Setup

### Lattice geometry

Two sublattices $A, B$ on the honeycomb lattice, with each $A$-site
having three nearest neighbours all on the $B$ sublattice. Take the
NN bond length as the unit of length. Primitive lattice vectors:

$$\mathbf{a}_1 \;=\; (\sqrt{3},\,0),
\qquad
\mathbf{a}_2 \;=\; \bigl(\tfrac{\sqrt{3}}{2},\,\tfrac{3}{2}\bigr).$$

Writing an $A$-site at the origin and its three $B$-neighbour
displacements,

$$\boldsymbol{\delta}_1 \;=\; (0,\,1),\qquad
  \boldsymbol{\delta}_2 \;=\; \bigl(\tfrac{\sqrt{3}}{2},\,-\tfrac{1}{2}\bigr),\qquad
  \boldsymbol{\delta}_3 \;=\; \bigl(-\tfrac{\sqrt{3}}{2},\,-\tfrac{1}{2}\bigr).
\tag{1}$$

Check consistency: $|\boldsymbol{\delta}_i| = 1$ for $i = 1, 2, 3$
(NN bond length, e.g. $|\boldsymbol{\delta}_2|^2 = 3/4 + 1/4 = 1$),
and $\boldsymbol{\delta}_1 + \boldsymbol{\delta}_2 +
\boldsymbol{\delta}_3 = (0 + \sqrt{3}/2 + (-\sqrt{3}/2),\;
1 + (-1/2) + (-1/2)) = (0, 0)$ (the three $B$-neighbours of a given
$A$-site are symmetric under $C_3$ rotation about the $A$-site).

The honeycomb lattice is **bipartite**: every bond connects $A$ to
$B$, so the adjacency matrix has no $A$-$A$ or $B$-$B$ entries.

### Brillouin zone

Reciprocal lattice vectors satisfy
$\mathbf{a}_i \cdot \mathbf{b}_j = 2\pi\,\delta_{ij}$. Solving,

$$\mathbf{b}_1 \;=\; \frac{2\pi}{\sqrt{3}}\,\bigl(1,\,-\tfrac{1}{\sqrt{3}}\bigr),
\qquad
\mathbf{b}_2 \;=\; \frac{2\pi}{\sqrt{3}}\,\bigl(0,\,\tfrac{2}{\sqrt{3}}\bigr).$$

The first Brillouin zone is a hexagon whose corners are the six
**$K$-points** of the lattice — see Step 4.

### Hamiltonian

Real-space tight-binding,

$$H \;=\; -t\sum_{\langle i, j\rangle}
        \bigl(c^{\dagger}_i c_j + c^{\dagger}_j c_i\bigr),
\tag{2}$$

with $\langle i, j\rangle$ over NN $A$-$B$ bonds. Equivalently, with
$\mathbf{R}$ indexing unit cells and $A_{\mathbf{R}}, B_{\mathbf{R}}$
the sublattice fermions,

$$H \;=\; -t\sum_{\mathbf{R}}\sum_{i=1}^{3}
        \bigl(A^{\dagger}_{\mathbf{R}}\,B_{\mathbf{R} + \boldsymbol{\delta}_i^{(c)}}
              \;+\; \text{h.c.}\bigr),$$

where $\boldsymbol{\delta}_i^{(c)}$ is the **cell-offset** vector for
the $i$-th NN bond (see Step 1 for the explicit list).

### Goal

Derive the Bloch Hamiltonian from (2), solve the $2\times 2$
eigenproblem at each $\mathbf{k}$, verify the Dirac-point
structure and linear expansion.

---

## Calculation

### Step 1 — Fourier transform to Bloch form

Of the three NN vectors $\boldsymbol{\delta}_i$ in (1),
$\boldsymbol{\delta}_1$ connects an $A$-site to the $B$-site in the
**same** unit cell, while $\boldsymbol{\delta}_2$ and
$\boldsymbol{\delta}_3$ connect to $B$-sites in **neighbouring**
unit cells offset by $-\mathbf{a}_1$ and $-\mathbf{a}_2$
respectively. Explicitly, with the convention that each unit cell
contains one $A$-site at position $\mathbf{R}$ and one $B$-site at
$\mathbf{R} + \boldsymbol{\delta}_1$, the three NN hoppings from the
$A$-site at $\mathbf{R}$ are

$$A_{\mathbf{R}} \leftrightarrow
\begin{cases}
B_{\mathbf{R}} & \text{(same cell, offset }\boldsymbol{\delta}_1\text{)} \\
B_{\mathbf{R} - \mathbf{a}_1} & \text{(cell }-\mathbf{a}_1\text{, offset }\boldsymbol{\delta}_2\text{)} \\
B_{\mathbf{R} - \mathbf{a}_2} & \text{(cell }-\mathbf{a}_2\text{, offset }\boldsymbol{\delta}_3\text{)}.
\end{cases}
\tag{3}$$

(One can verify (3) geometrically by drawing a honeycomb lattice and
labelling cells by $\mathbf{R} = n_1 \mathbf{a}_1 + n_2 \mathbf{a}_2$.)

Introduce Bloch operators

$$A_{\mathbf{k}} = \frac{1}{\sqrt{N_c}}\sum_{\mathbf{R}}
                    e^{-i\,\mathbf{k}\cdot\mathbf{R}}\,A_{\mathbf{R}},
\qquad
B_{\mathbf{k}} = \frac{1}{\sqrt{N_c}}\sum_{\mathbf{R}}
                    e^{-i\,\mathbf{k}\cdot\mathbf{R}}\,B_{\mathbf{R}},$$

where $N_c = L_x L_y$ is the number of unit cells and $\mathbf{k}$
runs over the Brillouin zone. Substitute (3) into (2) and use
$\sum_{\mathbf{R}} e^{i\mathbf{k}\cdot\mathbf{R}}\,e^{-i\mathbf{k}'
\cdot\mathbf{R}} = N_c\,\delta_{\mathbf{k}, \mathbf{k}'}$:

$$A^{\dagger}_{\mathbf{R}}\,B_{\mathbf{R} - \mathbf{a}_i}
 = \frac{1}{N_c}\sum_{\mathbf{k}, \mathbf{k}'}
   e^{i\mathbf{k}\cdot\mathbf{R}}\,e^{-i\mathbf{k}'\cdot(\mathbf{R} - \mathbf{a}_i)}\,
   A^{\dagger}_{\mathbf{k}}\,B_{\mathbf{k}'},$$

and summing over $\mathbf{R}$ gives $N_c\,\delta_{\mathbf{k},
\mathbf{k}'}$, leaving

$$\sum_{\mathbf{R}}A^{\dagger}_{\mathbf{R}}\,B_{\mathbf{R} - \mathbf{a}_i}
 = \sum_{\mathbf{k}} e^{-i\mathbf{k}\cdot(-\mathbf{a}_i)}
   A^{\dagger}_{\mathbf{k}}\,B_{\mathbf{k}}
 = \sum_{\mathbf{k}} e^{i\mathbf{k}\cdot\mathbf{a}_i}
   A^{\dagger}_{\mathbf{k}}\,B_{\mathbf{k}}.$$

For the $\boldsymbol{\delta}_1$ bond ($B$ in the same cell),
$\mathbf{a}_i = 0$ and the phase is $1$.

Assembling all three NN contributions,

$$H = -t\sum_{\mathbf{k}}\Bigl[
   \bigl(1 + e^{i\mathbf{k}\cdot\mathbf{a}_1} + e^{i\mathbf{k}\cdot\mathbf{a}_2}\bigr)
   A^{\dagger}_{\mathbf{k}}\,B_{\mathbf{k}} \;+\;\text{h.c.}\Bigr]
\;=\;
\sum_{\mathbf{k}}
\begin{pmatrix}A^{\dagger}_{\mathbf{k}} & B^{\dagger}_{\mathbf{k}}\end{pmatrix}
\begin{pmatrix} 0 & f(\mathbf{k})\\ f^{*}(\mathbf{k}) & 0\end{pmatrix}
\begin{pmatrix}A_{\mathbf{k}} \\ B_{\mathbf{k}}\end{pmatrix},$$

with the off-diagonal function

$$\boxed{\;
f(\mathbf{k}) \;=\; -t\,\bigl(1 + e^{i\mathbf{k}\cdot\mathbf{a}_1}
                                  + e^{i\mathbf{k}\cdot\mathbf{a}_2}\bigr).
\;}
\tag{4}$$

The diagonal entries vanish — this is the algebraic expression of
the bipartite (chiral) symmetry: the adjacency matrix connects only
$A \leftrightarrow B$, so the on-site block is zero.

### Step 2 — Eigenvalues: $E(\mathbf{k}) = \pm|f(\mathbf{k})|$

The $2\times 2$ Bloch Hamiltonian
$H(\mathbf{k}) = \begin{pmatrix}0 & f\\ f^{*} & 0\end{pmatrix}$ has
characteristic polynomial $\det(H - E\mathbb{I}) = E^2 - |f|^2$, so
the two bands are

$$E_{\pm}(\mathbf{k}) \;=\; \pm\,|f(\mathbf{k})|.
\tag{5}$$

The eigenvalues are symmetric about $E = 0$ for every
$\mathbf{k}$ — this is the chiral / sublattice symmetry
$\sigma^z H(\mathbf{k})\sigma^z = -H(\mathbf{k})$ with $\sigma^z =
\operatorname{diag}(1, -1)$ acting in the $A/B$ basis.

### Step 3 — Expand $|f(\mathbf{k})|^2$

From (4),

$$|f(\mathbf{k})|^2
 = t^2\,\bigl(1 + e^{i\mathbf{k}\cdot\mathbf{a}_1}
                  + e^{i\mathbf{k}\cdot\mathbf{a}_2}\bigr)\,
       \bigl(1 + e^{-i\mathbf{k}\cdot\mathbf{a}_1}
                  + e^{-i\mathbf{k}\cdot\mathbf{a}_2}\bigr).$$

Expand the product term by term. The diagonal $1 \cdot 1$,
$e^{i\theta_1}\cdot e^{-i\theta_1}$, $e^{i\theta_2}\cdot e^{-i\theta_2}$
each give $1$, contributing $3$. The off-diagonal products come in
conjugate pairs:

- $1 \cdot e^{-i\theta_1} + e^{i\theta_1}\cdot 1 = 2\cos\theta_1$,
- $1 \cdot e^{-i\theta_2} + e^{i\theta_2}\cdot 1 = 2\cos\theta_2$,
- $e^{i\theta_1}\cdot e^{-i\theta_2} + e^{i\theta_2}\cdot e^{-i\theta_1}
  = 2\cos(\theta_1 - \theta_2)$.

Introducing $\theta_j \equiv \mathbf{k}\cdot\mathbf{a}_j$ and using
$\cos(\theta_1 - \theta_2) = \cos(\theta_2 - \theta_1)$,

$$|f(\mathbf{k})|^2
 \;=\; t^2\,\bigl[\,3
        \;+\; 2\cos\theta_1 \;+\; 2\cos\theta_2
        \;+\; 2\cos(\theta_2 - \theta_1)\,\bigr].
\tag{6}$$

Taking the square root yields the boxed Main-result dispersion (5).

### Step 4 — Locate the Dirac points

$E_{\pm}(\mathbf{k}) = 0$ iff $f(\mathbf{k}) = 0$, iff
$1 + e^{i\theta_1} + e^{i\theta_2} = 0$. Interpreting this as the sum
of three unit complex numbers vanishing, the three phases must form
an equilateral triangle on the unit circle, i.e. be $120°$ apart:

$$\{0,\,\theta_1,\,\theta_2\} \;=\; \{0,\,2\pi/3,\,4\pi/3\}
\pmod{2\pi}.$$

Two solutions:

$$\theta_1 = \tfrac{2\pi}{3},\;\theta_2 = -\tfrac{2\pi}{3}
\qquad\text{(point }K\text{)},$$
$$\theta_1 = -\tfrac{2\pi}{3},\;\theta_2 = \tfrac{2\pi}{3}
\qquad\text{(point }K'\text{)}.$$

Invert $\theta_j = \mathbf{k}\cdot\mathbf{a}_j$ using the
reciprocal-lattice relation $\mathbf{k} = (\theta_1\,\mathbf{b}_1 +
\theta_2\,\mathbf{b}_2)/(2\pi)$:

$$\mathbf{K}
 \;=\; \tfrac{1}{3}\,\mathbf{b}_1 - \tfrac{1}{3}\,\mathbf{b}_2
 \;=\; \frac{2\pi}{3}\,\Bigl(\tfrac{1}{\sqrt{3}},\,-1\Bigr)
 \;\cdot\;(-1)\text{ sign by choice of unit cell},$$

or equivalently, using the NN vectors (1),

$$\mathbf{K} \;=\; \frac{2\pi}{3}\,\bigl(\tfrac{1}{\sqrt{3}},\,1\bigr),
\qquad
\mathbf{K'} \;=\; \frac{2\pi}{3}\,\bigl(\tfrac{1}{\sqrt{3}},\,-1\bigr).
\tag{7}$$

The other four corners of the hexagonal BZ are related to $\mathbf{K}$
and $\mathbf{K'}$ by reciprocal-lattice translations, so there are
only two inequivalent Dirac points. (To check (7), evaluate $\mathbf{K}
\cdot \mathbf{a}_1 = (2\pi/3)(\tfrac{1}{\sqrt{3}}\cdot\sqrt{3} +
1\cdot 0) = 2\pi/3$ and $\mathbf{K}\cdot\mathbf{a}_2 =
(2\pi/3)(\tfrac{1}{\sqrt{3}}\cdot\tfrac{\sqrt{3}}{2} + 1\cdot
\tfrac{3}{2}) = (2\pi/3)(1/2 + 3/2) = (2\pi/3)\cdot 2 = 4\pi/3 \equiv
-2\pi/3\;(\mathrm{mod}\;2\pi)$. ✓)

### Step 5 — Linear expansion near the Dirac point

Set $\mathbf{k} = \mathbf{K} + \mathbf{q}$ with $|\mathbf{q}| \ll
|\mathbf{K}|$, and Taylor-expand $f(\mathbf{K} + \mathbf{q})$.
From (4),

$$f(\mathbf{K} + \mathbf{q})
 \;=\; -t\,\bigl(1 + e^{i(\mathbf{K} + \mathbf{q})\cdot\mathbf{a}_1}
                    + e^{i(\mathbf{K} + \mathbf{q})\cdot\mathbf{a}_2}\bigr).$$

Using $f(\mathbf{K}) = 0$ (from Step 4) and Taylor-expanding each
exponential to first order in $\mathbf{q}$,

$$e^{i(\mathbf{K} + \mathbf{q})\cdot\mathbf{a}_j}
 = e^{i\mathbf{K}\cdot\mathbf{a}_j}\bigl(1 + i\mathbf{q}\cdot\mathbf{a}_j
 + O(|\mathbf{q}|^2)\bigr),$$

so

$$f(\mathbf{K} + \mathbf{q})
 = -t\,\bigl[\underbrace{1 + e^{i\theta_1^K} + e^{i\theta_2^K}}_{=\,0}
   \;+\; i\,(e^{i\theta_1^K}\,\mathbf{q}\cdot\mathbf{a}_1
              + e^{i\theta_2^K}\,\mathbf{q}\cdot\mathbf{a}_2)
   \;+\; O(|\mathbf{q}|^2)\bigr],$$

with $\theta_j^K = \mathbf{K}\cdot\mathbf{a}_j$ from Step 4:
$\theta_1^K = 2\pi/3$, $\theta_2^K = -2\pi/3$. Substituting
$e^{i\,2\pi/3} = -\tfrac{1}{2} + i\tfrac{\sqrt{3}}{2}$ and
$e^{-i\,2\pi/3} = -\tfrac{1}{2} - i\tfrac{\sqrt{3}}{2}$, and the
explicit $\mathbf{a}_1 = (\sqrt{3}, 0)$, $\mathbf{a}_2 =
(\sqrt{3}/2, 3/2)$:

$$e^{i\theta_1^K}\,\mathbf{a}_1
 = \bigl(-\tfrac{1}{2} + i\tfrac{\sqrt{3}}{2}\bigr)\,(\sqrt{3}, 0)
 = \bigl(-\tfrac{\sqrt{3}}{2} + i\tfrac{3}{2},\;0\bigr),$$
$$e^{i\theta_2^K}\,\mathbf{a}_2
 = \bigl(-\tfrac{1}{2} - i\tfrac{\sqrt{3}}{2}\bigr)\,\bigl(\tfrac{\sqrt{3}}{2},\,\tfrac{3}{2}\bigr)
 = \bigl(-\tfrac{\sqrt{3}}{4} - i\tfrac{3}{4},\;
        -\tfrac{3}{4} - i\tfrac{3\sqrt{3}}{4}\bigr).$$

Sum component-by-component:

$$e^{i\theta_1^K}\mathbf{a}_1 + e^{i\theta_2^K}\mathbf{a}_2
 \;=\; \bigl(-\tfrac{3\sqrt{3}}{4} + i\tfrac{3}{4},\;
              -\tfrac{3}{4} - i\tfrac{3\sqrt{3}}{4}\bigr).$$

Multiply by $-it$ (from $f = -t\cdot i(\ldots)$):

$$f(\mathbf{K} + \mathbf{q})
 \;\approx\; (-t) \cdot i \cdot \bigl(-\tfrac{3\sqrt{3}}{4}
              + i\tfrac{3}{4},\;
              -\tfrac{3}{4} - i\tfrac{3\sqrt{3}}{4}\bigr)\cdot\mathbf{q}$$
$$= \bigl(\tfrac{3}{4} + i\tfrac{3\sqrt{3}}{4},\;
          -\tfrac{3\sqrt{3}}{4} + i\tfrac{3}{4}\bigr) \cdot t\mathbf{q}\cdot(-i\cdot\frac{-1}{1})$$

Let me redo more cleanly. Let me write $f(\mathbf{K} + \mathbf{q})
= i t (\alpha\,q_x + \beta\,q_y) + O(q^2)$ with complex coefficients
$\alpha, \beta$. From the calculation above,

$$\alpha = -\,(\,e^{i\theta_1^K}\,a_{1,x} + e^{i\theta_2^K}\,a_{2,x}\,)
 = -\bigl(-\tfrac{\sqrt{3}}{2} + i\tfrac{3}{2} \;+\;
         -\tfrac{\sqrt{3}}{4} - i\tfrac{3}{4}\bigr)
 = \tfrac{3\sqrt{3}}{4} - i\tfrac{3}{4},$$

$$\beta = -\,(\,e^{i\theta_1^K}\,a_{1,y} + e^{i\theta_2^K}\,a_{2,y}\,)
 = -\bigl(0 \;+\; -\tfrac{3}{4} - i\tfrac{3\sqrt{3}}{4}\bigr)
 = \tfrac{3}{4} + i\tfrac{3\sqrt{3}}{4}.$$

The extra sign is from the outer $-t$ in (4). So

$$f(\mathbf{K} + \mathbf{q})
 \;\approx\; i\,t\,\bigl(\alpha\,q_x + \beta\,q_y\bigr).$$

Compute $|f|^2 = t^2\,|\alpha\,q_x + \beta\,q_y|^2$. Using

$$|\alpha|^2 = \bigl(\tfrac{3\sqrt{3}}{4}\bigr)^2 + \bigl(\tfrac{3}{4}\bigr)^2
 = \tfrac{27}{16} + \tfrac{9}{16} = \tfrac{36}{16} = \tfrac{9}{4},$$
$$|\beta|^2 = \bigl(\tfrac{3}{4}\bigr)^2 + \bigl(\tfrac{3\sqrt{3}}{4}\bigr)^2
 = \tfrac{9}{4},$$

and the cross term

$$\alpha^{*}\,\beta = \bigl(\tfrac{3\sqrt{3}}{4} + i\tfrac{3}{4}\bigr)
                     \bigl(\tfrac{3}{4} + i\tfrac{3\sqrt{3}}{4}\bigr)
 = \tfrac{9\sqrt{3}}{16} + i\tfrac{27}{16}
 + i\tfrac{9}{16} + i^2\,\tfrac{9\sqrt{3}}{16}
 = 0 + i\,\tfrac{36}{16} = i\,\tfrac{9}{4},$$

so $2\operatorname{Re}(\alpha^{*}\beta) = 0$. Therefore

$$|\alpha\,q_x + \beta\,q_y|^2
 = |\alpha|^2\,q_x^2 + 2\operatorname{Re}(\alpha^{*}\beta)\,q_x q_y
   + |\beta|^2\,q_y^2
 = \tfrac{9}{4}\,(q_x^2 + q_y^2) = \tfrac{9}{4}\,|\mathbf{q}|^2,$$

and

$$|f(\mathbf{K} + \mathbf{q})|^2
 \;\approx\; t^2\,\tfrac{9}{4}\,|\mathbf{q}|^2,
\qquad
|f(\mathbf{K} + \mathbf{q})|
 \;\approx\; \tfrac{3 t}{2}\,|\mathbf{q}|.$$

The dispersion near the $K$-point is therefore

$$\boxed{\;
E_{\pm}(\mathbf{K} + \mathbf{q})
 \;=\; \pm\,v_F\,|\mathbf{q}| \;+\; O(|\mathbf{q}|^2),
\qquad
v_F \;=\; \frac{3 t}{2}.
\;}
\tag{8}$$

The same expansion at $\mathbf{K}'$ gives the same $v_F$ with the
chirality of the Dirac cone inverted (Castro Neto et al. 2009
eq. (2.9)). At half-filling (one electron per site) the lower band
is filled and the upper band is empty; the chemical potential sits
exactly at the Dirac point and the system is a semimetal.

### Step 6 — Finite-size spectrum (PBC)

On an $L_x \times L_y$ lattice with PBC, allowed momenta are

$$\mathbf{k}_{m, n}
 = \frac{m}{L_x}\,\mathbf{b}_1 + \frac{n}{L_y}\,\mathbf{b}_2,
\qquad m = 0,\dots,L_x - 1,\;\; n = 0,\dots,L_y - 1,$$

and $\mathbf{k}_{m, n}\cdot\mathbf{a}_1 = 2\pi m/L_x$,
$\mathbf{k}_{m, n}\cdot\mathbf{a}_2 = 2\pi n/L_y$. Substitute into
(6):

$$E_{m, n}
 = \pm t\sqrt{3 + 2\cos\!\bigl(\tfrac{2\pi m}{L_x}\bigr)
              + 2\cos\!\bigl(\tfrac{2\pi n}{L_y}\bigr)
              + 2\cos\!\bigl(\tfrac{2\pi n}{L_y} - \tfrac{2\pi m}{L_x}\bigr)}.
\tag{9}$$

This is the formula used by the QAtlas test utility `bloch_tb_spectrum`
in `test/util/bloch.jl`; the agreement between (9) and a real-space ED
at every $(L_x, L_y)$ is one of the package's standing cross-checks
(`test/verification/test_bloch_generic.jl`).

### Step 7 — Limiting-case checks

**(i) Half-filling / Dirac point.** At half-filling, the lower band
$E_-(\mathbf{k}) = -|f(\mathbf{k})|$ is completely filled. The
density of states vanishes linearly at $E = 0$ (a hallmark of 2D
Dirac fermions), so the chemical potential sits exactly at the
Dirac point and transport is gapless. For the QAtlas
`Honeycomb` model, the $(m, n)$ momentum grid includes a zero
eigenvalue exactly when $(L_x, L_y)$ commensurates the $K$-point
position, i.e. when $L_x$ and $L_y$ are both multiples of $3$. On
non-commensurate lattices the gap closes only in the $L_x, L_y \to
\infty$ limit.

**(ii) Band edges.** The bandwidth is $W = 2\max_{\mathbf{k}} |f| =
6 t$ (at $\mathbf{k} = 0$, from $|f(0)| = |-t\cdot 3| = 3 t$). So
the spectrum spans $[-3 t, +3 t]$ symmetrically.

**(iii) Bipartite / chiral check.** Under
$(A_{\mathbf{k}}, B_{\mathbf{k}}) \to (A_{\mathbf{k}},
-B_{\mathbf{k}})$, $H(\mathbf{k}) \to -H(\mathbf{k})$ (the
off-diagonal element flips sign). This is the algebraic statement
of bipartiteness, and it forces the $E_{+}(\mathbf{k}) =
-E_{-}(\mathbf{k})$ pairing that we derived from the $2\times 2$
structure directly.

---

## References

- P. R. Wallace, *The band theory of graphite*, Phys. Rev. **71**,
  622 (1947). Original tight-binding derivation of the honeycomb
  dispersion; eqs. (1)–(7) are equivalent to (4)–(8) here.
- A. H. Castro Neto, F. Guinea, N. M. R. Peres, K. S. Novoselov,
  A. K. Geim, *The electronic properties of graphene*, Rev. Mod.
  Phys. **81**, 109 (2009). Comprehensive review; eq. (2.5) is the
  Bloch Hamiltonian (4), eq. (2.9) is the Dirac-cone expansion (8).
- N. W. Ashcroft and N. D. Mermin, *Solid State Physics* (Saunders,
  1976), Ch. 8. Textbook derivation of Bloch theorem and unit-cell
  Fourier transforms.

## Used by

- [Honeycomb tight-binding model](../models/quantum/tightbinding/honeycomb.md) —
  `fetch(Honeycomb(; t, Lx, Ly), TightBindingSpectrum(), …)` returns
  the finite-size spectrum (9).
- [Bloch Hamiltonian method page](../methods/bloch-hamiltonian/index.md) —
  this note is the canonical worked example for 2D multi-sublattice
  tight-binding bands.
- [Lieb flat-band note](bloch-lieb-flat-band.md) — contrast: both
  honeycomb and Lieb are bipartite, but only Lieb has an odd
  sublattice-imbalance and therefore a flat band at $E = 0$.
- [Kagome flat-band note](bloch-kagome-flat-band.md) — contrast:
  kagome has a flat band from destructive interference on a
  non-bipartite three-sublattice unit cell.

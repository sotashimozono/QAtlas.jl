# Yang's Spontaneous Magnetization via Toeplitz Determinant

## Main result

For the isotropic 2D classical Ising model on the square lattice at
reduced coupling $K = \beta J$ below the critical coupling
$K_c = \tfrac{1}{2}\ln(1 + \sqrt{2})$, the spontaneous
magnetization — defined as the long-distance limit of the
spin-spin correlation function —

$$M(T) \;\equiv\; \lim_{|\mathbf{r}|\to\infty}\sqrt{\,\bigl\langle \sigma_{\mathbf{0}}\,\sigma_{\mathbf{r}}\bigr\rangle\,}$$

has the exact closed form

$$\boxed{\;
M(T) \;=\; \Bigl(\,1 - \sinh^{-4}(2\beta J)\,\Bigr)^{1/8},
\qquad T < T_c.
\;}$$

The $(T_c - T)^{1/8}$ vanishing at criticality identifies the
order-parameter critical exponent

$$\boxed{\;\beta = \tfrac{1}{8}.\;}$$

The derivation proceeds via three ingredients:

1. **Pfaffian reduction** (Kaufman 1949, Montroll–Potts–Ward 1963):
   the spin-spin correlation function is an $n\times n$
   **Toeplitz determinant** $D_n[a]$ with a known symbol
   $a(\theta)$.
2. **Strong Szegő limit theorem** (Szegő 1915; sharp form
   Kac 1954): the $n\to\infty$ limit of $D_n[a]$ is
   $\exp\!\bigl[\sum_{k\ge 1} k\,c_k c_{-k}\bigr]$, with
   $c_k$ the Fourier coefficients of $\log a(\theta)$.
3. **Wiener–Hopf factorisation** of $a(\theta)$ on the
   symmetric-split form, followed by explicit Fourier-coefficient
   evaluation.

Substituting in a closed-form relation $\alpha^{2} =
\sinh^{-4}(2\beta J)$ between the Wiener–Hopf parameter $\alpha$
and the Ising coupling $K = \beta J$ (Yang 1952 §5) yields the
boxed result.

---

## Setup

### 2D classical Ising model

On an infinite square lattice with reduced coupling $K = \beta J$
and PBC,

$$Z = \sum_{\{s\}} \exp\!\Bigl(K\sum_{\langle i, j\rangle} s_i s_j\Bigr),
\qquad s_v \in \{\pm 1\}.$$

For $K > K_c = \tfrac{1}{2}\ln(1 + \sqrt{2})$ (equivalently
$\sinh(2K) > 1$, $T < T_c$) the system is in the ordered
ferromagnetic phase; the spontaneous magnetization $M(T)$ is the
standard order parameter.

### Correlation-function definition of $M(T)$

By translational symmetry and cluster decomposition, in the ordered
phase

$$\langle\sigma_{\mathbf{0}}\,\sigma_{\mathbf{r}}\rangle
 \;\xrightarrow[|\mathbf{r}|\to\infty]{}\; \langle\sigma\rangle^{2}
 \;=\; M(T)^{2},$$

so

$$M(T) \;=\; \sqrt{\lim_{|\mathbf{r}|\to\infty}
             \langle\sigma_{\mathbf{0}}\,\sigma_{\mathbf{r}}\rangle}.
\tag{1}$$

### Goal

Reduce the right-hand side of (1) to a Toeplitz determinant $D_n[a]$
with an explicit symbol, take the $n\to\infty$ limit via Szegő's
theorem, evaluate the resulting Fourier-coefficient sum in closed
form, and identify $M(T)$ as a function of $\beta J$.

---

## Calculation

### Step 1 — Correlation function as a Toeplitz determinant

The 2D Ising correlator along a diagonal direction
(Kaufman 1949, MPW 1963) is expressible as the determinant of an
$n\times n$ matrix whose $(j, k)$ entry depends only on $j - k$:

$$\bigl\langle\sigma_{0,0}\,\sigma_{n,n}\bigr\rangle
 \;=\; D_{n}[a] \;:=\; \det\bigl[\,a_{j - k}\,\bigr]_{j, k = 0}^{n - 1}.
\tag{2}$$

A matrix of this form is called a **Toeplitz matrix**, and its
determinant $D_n[a]$ is a **Toeplitz determinant** with
*symbol* $a(\theta) = \sum_{m\in\mathbb{Z}} a_m\,e^{-i m\theta}$.

The reduction (2) from the original Ising spin variables to a
fermion-determinant form goes through:

(i) Jordan–Wigner transformation of the Kaufman transfer matrix,
    producing a free-fermion quadratic form. The symbol of the
    resulting covariance matrix is
    $a(\theta) = \cos\phi(\theta)$ with
    $\phi(\theta) = \phi(\theta; K_1, K_2)$ the Bogoliubov angle of
    the $k$-th mode.

(ii) Wick contraction of the spin–spin correlator at sites
     $(0, 0)$ and $(n, n)$ using the free-fermion 2-point function.
     The result is a Pfaffian, which for a suitable ordering of
     Majorana indices collapses to a determinant of Toeplitz form
     (Montroll–Potts–Ward 1963 Appendix B).

The full reduction is a tour-de-force spread across Kaufman (1949),
Yang (1952), and Montroll–Potts–Ward (1963), totalling some
80 journal pages. For the purpose of this note we **take (2) as the
starting point** and focus on the Szegő evaluation, which is where
the critical exponent $\beta = 1/8$ ultimately emerges. The reader
consulting (e.g.) McCoy–Wu (1973), Ch. VII, finds the reduction
worked out end-to-end.

### Step 2 — The symbol $a(\theta)$ for $T < T_c$

After the Wiener–Hopf factorisation of the raw Bogoliubov-angle
symbol, the effective Toeplitz symbol for the diagonal correlator on
the isotropic lattice at $T < T_c$ has the remarkably clean form
(McCoy–Wu 1973, eq. VII.2.36a)

$$\boxed{\;
a(\theta) \;=\; \left(\frac{1 - \alpha\,e^{+i\theta}}
                            {1 - \alpha\,e^{-i\theta}}\right)^{1/2},
\qquad 0 \le \alpha < 1,
\;}
\tag{3}$$

with the Wiener–Hopf parameter

$$\boxed{\;
\alpha^{2} \;=\; \sinh^{-4}(2\beta J),
\qquad \alpha \in [0, 1)\text{ for }T < T_c.
\;}
\tag{4}$$

Note that $\alpha \to 0$ as $T \to 0$ ($\sinh(2\beta J) \to \infty$)
and $\alpha \to 1^{-}$ as $T \to T_c^{-}$
($\sinh(2\beta_c J) = 1$). The symbol is **unimodular**
($|a(\theta)| = 1$ for all real $\theta$), which implies $G \equiv
\exp(c_0) = 1$ in the Szegő theorem below.

### Step 3 — Strong Szegő limit theorem

**Theorem (Szegő 1915; sharp form Kac 1954).** Let $a(\theta)$ be
a non-vanishing function on the unit circle with smooth logarithm,
and let

$$\log a(\theta) \;=\; \sum_{k\in\mathbb{Z}} c_k\,e^{i k\theta},
\qquad c_k = \frac{1}{2\pi}\int_{-\pi}^{\pi}
                e^{-i k\theta}\,\log a(\theta)\,d\theta.$$

Then as $n\to\infty$,

$$\boxed{\;
D_n[a] \;=\; G^{n}\,E\,\bigl(1 + o(1)\bigr),\qquad
G = e^{c_0},\qquad
E = \exp\!\Bigl[\sum_{k=1}^{\infty} k\,c_k\,c_{-k}\Bigr].
\;}
\tag{5}$$

The sum $\sum k\,c_k c_{-k}$ converges absolutely whenever
$\sum_{k} |k|\,|c_k|^{2} < \infty$ — guaranteed here by the smooth
rational form (3), which has $c_k$ decaying geometrically as
$\alpha^{|k|}/|k|$ (Step 4).

(References: the original Szegő 1915 paper proves $D_n[a] \sim
G^n$; the $E$-constant $\exp(\sum k c_k c_{-k})$ is due to
Szegő 1952 and in its sharp form Kac 1954. Modern pedagogical
treatments: Böttcher–Silbermann 1999, *Analysis of Toeplitz
Operators*, Ch. 1; Simon 2005, *Orthogonal Polynomials on the Unit
Circle*, Ch. 6.)

For the Ising symbol (3), $a$ is unimodular so $c_0 = 0$ and
$G = 1$. Therefore $D_n[a] \to E$ as $n\to\infty$.

### Step 4 — Fourier coefficients of $\log a(\theta)$

From (3),

$$\log a(\theta)
 \;=\; \tfrac{1}{2}\log(1 - \alpha e^{+i\theta})
      \;-\; \tfrac{1}{2}\log(1 - \alpha e^{-i\theta}).$$

Each logarithm is expanded via $\log(1 - x) = -\sum_{k\ge 1}
x^{k}/k$ (convergent for $|\alpha| < 1$):

$$\log(1 - \alpha e^{+i\theta})
 \;=\; -\sum_{k=1}^{\infty}\frac{\alpha^{k}}{k}\,e^{+i k\theta},$$

$$\log(1 - \alpha e^{-i\theta})
 \;=\; -\sum_{k=1}^{\infty}\frac{\alpha^{k}}{k}\,e^{-i k\theta}.$$

Therefore

$$\log a(\theta)
 \;=\; -\tfrac{1}{2}\sum_{k=1}^{\infty}\frac{\alpha^{k}}{k}\,e^{+i k\theta}
       \;+\; \tfrac{1}{2}\sum_{k=1}^{\infty}\frac{\alpha^{k}}{k}\,e^{-i k\theta},$$

reading off the Fourier coefficients:

$$\boxed{\;
c_{+k} = -\frac{\alpha^{k}}{2 k},\qquad
c_{-k} = +\frac{\alpha^{k}}{2 k}\qquad (k \ge 1),\qquad c_{0} = 0.
\;}
\tag{6}$$

The zero-th coefficient $c_0 = 0$ confirms $G = 1$ (geometric mean
of a unimodular symbol equals unity). Fourier coefficients decay as
$\alpha^{|k|}/(2|k|)$; since $|\alpha| < 1$ the series
$\sum k |c_k|^{2} \sim \sum \alpha^{2k}/k$ converges.

### Step 5 — Evaluate the Szegő constant $E$

Substitute (6) into the constant $E$ of (5):

$$\sum_{k=1}^{\infty} k\,c_{k}\,c_{-k}
 \;=\; \sum_{k=1}^{\infty} k\,\Bigl(-\frac{\alpha^{k}}{2 k}\Bigr)\Bigl(+\frac{\alpha^{k}}{2 k}\Bigr)
 \;=\; -\,\frac{1}{4}\sum_{k=1}^{\infty}\frac{\alpha^{2k}}{k}.$$

The inner sum is the standard Taylor series for $-\log(1 - x)$ at
$x = \alpha^{2}$:

$$\sum_{k=1}^{\infty}\frac{\alpha^{2k}}{k} \;=\; -\log(1 - \alpha^{2}).$$

Therefore

$$\sum_{k=1}^{\infty} k\,c_{k}\,c_{-k}
 \;=\; -\,\frac{1}{4}\cdot\bigl(-\log(1 - \alpha^{2})\bigr)
 \;=\; \tfrac{1}{4}\log(1 - \alpha^{2}),$$

and the Szegő constant is

$$\boxed{\;
E \;=\; \exp\!\Bigl[\,\tfrac{1}{4}\log(1 - \alpha^{2})\,\Bigr]
     \;=\; (1 - \alpha^{2})^{1/4}.
\;}
\tag{7}$$

### Step 6 — Assemble $M(T)$

From Step 3 and (7),

$$\lim_{n\to\infty} D_n[a] \;=\; E \;=\; (1 - \alpha^{2})^{1/4}.$$

The defining relation (1) expresses the magnetization as the square
root of this limit:

$$M(T)^{2}
 \;=\; \lim_{|\mathbf{r}|\to\infty}
    \bigl\langle\sigma_{\mathbf{0}}\,\sigma_{\mathbf{r}}\bigr\rangle
 \;=\; (1 - \alpha^{2})^{1/4}.$$

Substituting (4),

$$M(T)^{2} \;=\; \bigl(1 - \sinh^{-4}(2\beta J)\bigr)^{1/4}.$$

Taking the positive square root (magnetization is non-negative for
$T < T_c$),

$$M(T) \;=\; \bigl(1 - \sinh^{-4}(2\beta J)\bigr)^{1/8}.
\tag{8}$$

This is the Main-result formula.

### Step 7 — Critical exponent $\beta = 1/8$

Near criticality $T \lesssim T_c$, expand $\alpha^{2} =
\sinh^{-4}(2K)$ in small $\delta K \equiv K - K_c > 0$ (recall
$K \propto \beta \propto 1/T$, so $\delta K > 0 \Leftrightarrow
T < T_c$). Using

$$\sinh(2K) \;=\; \sinh(2 K_c + 2\,\delta K)
 \;=\; \sinh(2 K_c)\cosh(2\,\delta K)
    + \cosh(2 K_c)\sinh(2\,\delta K),$$

and the critical-point identities $\sinh(2 K_c) = 1$,
$\cosh^{2}(2 K_c) = 1 + \sinh^{2}(2 K_c) = 2$, so
$\cosh(2 K_c) = \sqrt{2}$:

$$\sinh(2K) \;=\; 1\cdot\bigl(1 + 2(\delta K)^{2} + \dots\bigr)
                   + \sqrt{2}\cdot\bigl(2\,\delta K + O((\delta K)^{3})\bigr)
 \;=\; 1 + 2\sqrt{2}\,\delta K + O((\delta K)^{2}).$$

Therefore

$$\sinh^{-4}(2K)
 \;=\; \bigl(1 + 2\sqrt{2}\,\delta K + \dots\bigr)^{-4}
 \;=\; 1 - 8\sqrt{2}\,\delta K + O((\delta K)^{2}),$$

and

$$1 - \alpha^{2}
 \;=\; 1 - \sinh^{-4}(2 K)
 \;=\; 8\sqrt{2}\,\delta K + O((\delta K)^{2}).$$

Converting to temperature: $K = J/(k_B T)$, so $\delta K =
-J\,\delta T/(k_B T^{2}) \approx -J(T - T_c)/(k_B T_c^{2})$, and
$\delta K > 0 \Leftrightarrow T < T_c$. Hence

$$1 - \alpha^{2} \;=\; \frac{8\sqrt{2}\,J}{k_B T_c^{2}}(T_c - T)
                      + O((T_c - T)^{2}).$$

From (8),

$$M(T) \;\sim\; (T_c - T)^{1/8}
\qquad (T \to T_c^{-}),$$

identifying the order-parameter critical exponent as $\beta = 1/8$.

### Step 8 — Limiting-case checks

**(i) $T = 0$ saturation.** At $T \to 0$, $K \to \infty$ and
$\sinh(2 K) \to \infty$, so $\alpha^{2} = \sinh^{-4}(2K) \to 0$ and
$M \to 1^{-}$. In fact $M(0) = 1$ is the saturated ferromagnetic
limit — all spins aligned in the ground state, $\langle\sigma\rangle
= \pm 1$.

**(ii) $T = T_c$ criticality.** At $K = K_c$,
$\sinh(2 K_c) = 1$, so $\alpha^{2} = 1$ and $M(T_c) =
(1 - 1)^{1/8} = 0$. Continuity of $M$ at $T_c$ from below confirms
the second-order nature of the phase transition.

**(iii) Universality-class cross-check.** The 2D Ising CFT has
primary-operator scaling dimensions
$\Delta_{\sigma} = 1/8$ and $\Delta_{\varepsilon} = 1$, and a
short-distance OPE
$\sigma(\mathbf{r})\,\sigma(\mathbf{0}) \sim
C/|\mathbf{r}|^{2\Delta_{\sigma}}$ at the critical point. The
exponent $\eta = 2\Delta_{\sigma} = 1/4$ (Fisher exponent for the
critical correlator) matches our $M \propto (T_c - T)^{1/8}$
through the scaling relation

$$\beta = \frac{\nu}{2}(d - 2 + \eta)$$

with $\nu = 1$ (correlation-length exponent of 2D Ising),
$d = 2$, $\eta = 1/4$: $\beta = (1/2)(0 + 1/4) = 1/8$. ✓

Cross-check with Ising-scaling relations derived in
[`ising-scaling-relations`](ising-scaling-relations.md): Widom
$\gamma = \beta(\delta - 1)$ with $\delta = 15$ gives
$\gamma = 1/8 \cdot 14 = 7/4$, matching the Onsager-susceptibility
result.

---

## Remark on why $\beta = 1/8$

The exponent $\beta = 1/8$ traces back to the $(1/2) \times (1/2) ×
(1/2) = 1/8$ structure of Step 6: the `1/2` from taking the square
root in (1) compounds with the `1/4` from the Szegő exponent (7),
which is itself `1/2 × 1/2` from the Wiener–Hopf square-root
factorisation of the symbol (3). In the closely related $c = 1/2$
Ising minimal-model CFT interpretation, the same factor
$2\Delta_{\sigma} = 1/4$ appears as the scaling dimension of the
two-point function $\langle\sigma\sigma\rangle$ — see
[`ising-cft-primary-operators`](ising-cft-primary-operators.md)
and [`ising-scaling-relations`](ising-scaling-relations.md) for
the CFT-based derivation.

---

## References

- C. N. Yang, *The spontaneous magnetization of a two-dimensional
  Ising model*, Phys. Rev. **85**, 808 (1952). Original derivation.
  Sections III–V give the Toeplitz reduction; §VI evaluates the
  limit via what is now called the strong Szegő theorem.
- L. Onsager, *Discussion remarks*, Nuovo Cimento Supp. **6**, 261
  (1949). The announcement of the result $M = (1 -
  \sinh^{-4}(2K))^{1/8}$ predates Yang's derivation; Onsager never
  published his proof.
- B. Kaufman, *Crystal statistics II: Partition function evaluated
  by spinor analysis*, Phys. Rev. **76**, 1232 (1949). Transfer-
  matrix diagonalisation via JW / Clifford algebra; underlies the
  Pfaffian reduction of (2).
- E. W. Montroll, R. B. Potts, J. C. Ward, *Correlations and
  spontaneous magnetization of the two-dimensional Ising model*,
  J. Math. Phys. **4**, 308 (1963). Systematic derivation of the
  Toeplitz form (2) starting from Kaufman's result.
- B. M. McCoy, T. T. Wu, *The Two-Dimensional Ising Model*,
  Harvard University Press (1973). Chapter VII is a textbook
  treatment of the Toeplitz reduction and the Wiener–Hopf
  factorisation. Equation (VII.2.36a) is our (3), and eq. (5.13)
  is our (4).
- G. Szegő, *Ein Grenzwertsatz über die Toeplitzschen
  Determinanten einer reellen positiven Funktion*, Math. Ann.
  **76**, 490 (1915). Original weak Szegő theorem (leading $G^n$).
- G. Szegő, *On certain Hermitian forms associated with the
  Fourier series of a positive function*, Comm. Sém. Math. Univ.
  Lund, tome supplémentaire (1952), 228. The sharp form including
  the constant $E = \exp(\sum k c_k c_{-k})$.
- M. Kac, *Toeplitz matrices, translation kernels and a related
  problem in probability theory*, Duke Math. J. **21**, 501 (1954).
  Independent rigorous proof of the sharp Szegő theorem.
- A. Böttcher, B. Silbermann, *Analysis of Toeplitz Operators*,
  2nd ed. (Springer, 2006), Ch. 1–3. Modern pedagogical treatment.

## Used by

- [IsingSquare model page](../models/classical/ising-square.md) —
  `fetch(IsingSquare(; J, Lx, Ly), SpontaneousMagnetization();
  β)` returns exactly (8), with `β = 1/8` as the order-parameter
  exponent.
- [Ising universality class](../universalities/ising.md) — the
  exponent $\beta = 1/8$ is one of the four independent 2D Ising
  critical exponents. The other three ($\alpha = 0\text{ (log)}$,
  $\gamma = 7/4$, $\nu = 1$) follow by scaling relations
  ([`ising-scaling-relations`](ising-scaling-relations.md)).
- [Kramers–Wannier duality note](kramers-wannier-duality.md) —
  $T_c$ location at the self-dual point; the exponent calculation
  here confirms the critical point is a second-order transition.

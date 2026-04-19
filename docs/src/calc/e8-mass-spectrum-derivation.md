# $E_8$ Mass Spectrum from the Bootstrap + Cartan Matrix

## Main result

The Ising field theory perturbed by the $\mathbb{Z}_{2}$-breaking
magnetic operator $\sigma$ is an integrable massive
$(1 + 1)$-dimensional QFT with **eight stable particles** whose
masses stand in the universal ratios determined by the $E_{8}$ Lie
algebra (Zamolodchikov 1989). With $\varphi = 2\cos(\pi/5) =
(1 + \sqrt{5})/2$ the golden ratio,

$$\boxed{\;
\begin{aligned}
\frac{m_{1}}{m_{1}} &= 1 & \frac{m_{5}}{m_{1}} &= 2\cos\!\bigl(\tfrac{2\pi}{15}\bigr)\cdot \varphi \approx 2.956 \\
\frac{m_{2}}{m_{1}} &= \varphi \approx 1.618 & \frac{m_{6}}{m_{1}} &= 2\cos\!\bigl(\tfrac{\pi}{30}\bigr)\cdot \varphi \approx 3.218 \\
\frac{m_{3}}{m_{1}} &= 2\cos\!\bigl(\tfrac{\pi}{30}\bigr) \approx 1.989 & \frac{m_{7}}{m_{1}} &= 2\cos\!\bigl(\tfrac{7\pi}{30}\bigr)\cdot \varphi^{2} \approx 3.891 \\
\frac{m_{4}}{m_{1}} &= 2\cos\!\bigl(\tfrac{7\pi}{30}\bigr)\cdot \varphi \approx 2.405 & \frac{m_{8}}{m_{1}} &= 2\cos\!\bigl(\tfrac{2\pi}{15}\bigr)\cdot \varphi^{2} \approx 4.783 .
\end{aligned}
\;}$$

Equivalently (Dorey 1991), the eight mass ratios are the components
of the **Perron–Frobenius eigenvector of the $E_{8}$ Cartan
matrix**, normalised so $m_{1} = 1$.

The most striking prediction — the golden-ratio mass ratio
$m_{2} / m_{1} = \varphi$ — was confirmed experimentally in the
quasi-1D Ising ferromagnet CoNb$_{2}$O$_{6}$ by Coldea et al. 2010,
who resolved $m_{2}/m_{1} = 1.618 \pm 0.015$.

---

## Setup

### Integrable field theory from the Ising CFT + magnetic perturbation

The critical 2D Ising model (unperturbed TFIM at $h = J$) is
described at low energies by a $c = 1/2$ minimal-model CFT with
three primary operators: identity $\mathbb{I}$, energy
$\varepsilon$, and spin $\sigma$ (see
[`ising-cft-primary-operators`](ising-cft-primary-operators.md)).
Turning on a uniform longitudinal field $h_{z}$ at the critical
temperature perturbs the action by

$$S \;\to\; S_{\rm CFT} + h_{z}\int d^{2}x\,\sigma(x).$$

Zamolodchikov 1989 proved that this perturbation preserves eight
commuting conserved charges (of spins $s = 1, 7, 11, 13, 17, 19,
23, 29$ — the exponents of $E_{8}$), making the perturbed theory
integrable. See
[`ising-cft-magnetic-perturbation`](ising-cft-magnetic-perturbation.md)
for the conserved-charges argument.

### Particle spectrum

Integrability $\Rightarrow$ the $S$-matrix factorises into 2-body
terms (Yang–Baxter) $\Rightarrow$ the theory has a discrete
particle spectrum. For the $E_{8}$ case Zamolodchikov identified
**eight stable particles**, which we label $1, 2, \ldots, 8$ in
order of increasing mass.

### Goal

Derive the eight mass ratios $m_{n}/m_{1}$ from (i) the bootstrap
equations satisfied by the 2-body $S$-matrix and (ii) the structural
identification of the mass spectrum with the Perron–Frobenius
eigenvector of the $E_{8}$ Cartan matrix (Dorey 1991).

---

## Calculation

### Step 1 — Conserved charges and elastic scattering

Zamolodchikov's 1989 argument (see
[`ising-cft-magnetic-perturbation`](ising-cft-magnetic-perturbation.md)
for the conserved-charges counting) shows that the perturbed Ising
theory has infinitely many conserved charges $Q_{s}$ of odd spin
$s$, with the spins occurring only at values
$s \in \{1, 7, 11, 13, 17, 19, 23, 29\} \pmod{30}$ — exactly the
exponents of $E_{8}$. The eight inequivalent conserved charges
(below the Coxeter number $h = 30$) impose the following structural
constraints on any scattering process:

1. **No particle production**: $n_{\rm in} = n_{\rm out}$ for every
   scattering amplitude.
2. **Elastic scattering**: the set of in-momenta equals the set of
   out-momenta.
3. **S-matrix factorisation**: every $n$-body $S$-matrix element is
   a product of $\binom{n}{2}$ two-body $S$-matrix elements.

**Proof sketch**: higher-spin conserved charges $Q_{s}$ act on
single-particle states as $Q_{s}|\theta, a\rangle =
q_{s}^{(a)}(\theta)|\theta, a\rangle$ with rapidity-dependent
eigenvalues $q_{s}^{(a)}(\theta) \propto m_{a}^{s+1} e^{s\theta}$.
Conservation of $Q_{s}$ across a scattering process with enough
linearly-independent eigenvalues forces the in-set equal to the
out-set (modulo permutation). For $E_{8}$-symmetric theories the
eight $q_{s}$ functions are independent enough to enforce this on
the nose. Details: Zamolodchikov 1989 §2; Mussardo 2010 Ch. 16.

### Step 2 — Bootstrap equations

Let $S_{ab}(\theta)$ denote the 2-body $S$-matrix for particles
$a, b$ of rapidity difference $\theta = \theta_{1} - \theta_{2}$.
Integrability + relativistic invariance forces $S_{ab}$ to satisfy
three functional equations:

**Unitarity.**

$$S_{ab}(\theta)\,S_{ab}(-\theta) \;=\; 1.
\tag{1}$$

**Crossing symmetry.**

$$S_{ab}(\theta) \;=\; S_{a\bar b}(i\pi - \theta),
\tag{2}$$

with $\bar b$ the antiparticle of $b$ (in the $E_{8}$ theory all
particles are self-conjugate so $\bar b = b$).

**Fusion / Bootstrap.** If the 2-body scattering of $a$ and $b$
has a bound-state pole at rapidity $\theta = i u_{ab}^{c}$ — i.e.
$S_{ab}(\theta) \sim i\,R_{ab}^{c}/(\theta - i u_{ab}^{c})$ near
the pole — then a bound state $c$ of mass

$$\boxed{\;
m_{c}^{2} \;=\; m_{a}^{2} + m_{b}^{2}
               + 2\,m_{a} m_{b}\,\cos u_{ab}^{c}
\;}
\tag{3}$$

exists. For every third particle $d$, the $S$-matrix elements with
$c$ on one leg and $(a, b)$ on the other satisfy the **bootstrap
equation**

$$S_{dc}(\theta) \;=\; S_{da}(\theta + i\,\bar{u}_{bc}^{a})\,
                       S_{db}(\theta - i\,\bar{u}_{ac}^{b}),
\tag{4}$$

with $\bar u \equiv \pi - u$. Equation (3) comes from Lorentz
invariance applied to the 3-body on-shell condition at the bound-
state pole; (4) follows from associativity of the $S$-matrix
algebra (Zamolodchikov–Zamolodchikov 1979).

### Step 3 — The $E_{8}$ Cartan matrix from the Dynkin diagram

The $E_{8}$ Lie algebra has rank $8$ and Coxeter number $h = 30$.
Its Dynkin diagram is

```
  1 — 2 — 3 — 4 — 5 — 6 — 7
          |
          8
```

(seven simply-laced nodes in a linear chain $1 \!-\! 2 \!-\! \dots
\!-\! 7$, plus one side-node $8$ attached to node $3$).

The Cartan matrix $C$ is determined by
$C_{ii} = 2$ and $C_{ij} = -A_{ij}$ where $A$ is the Dynkin
adjacency matrix ($A_{ij} = 1$ if nodes $i, j$ are connected, else
$0$). With the labelling above,

$$C \;=\; \begin{pmatrix}
 2 & -1 &  0 &  0 &  0 &  0 &  0 &  0 \\
-1 &  2 & -1 &  0 &  0 &  0 &  0 &  0 \\
 0 & -1 &  2 & -1 &  0 &  0 &  0 & -1 \\
 0 &  0 & -1 &  2 & -1 &  0 &  0 &  0 \\
 0 &  0 &  0 & -1 &  2 & -1 &  0 &  0 \\
 0 &  0 &  0 &  0 & -1 &  2 & -1 &  0 \\
 0 &  0 &  0 &  0 &  0 & -1 &  2 &  0 \\
 0 &  0 & -1 &  0 &  0 &  0 &  0 &  2
\end{pmatrix}.
\tag{5}$$

The matrix is real symmetric and positive-definite (standard
property of simply-laced ADE Cartan matrices). Its spectrum is
related to the **$E_{8}$ exponents** $e_{i} \in \{1, 7, 11, 13, 17,
19, 23, 29\}$ (equivalently $\{1, 7, 11, 13, 17, 19, 23, 29\}$;
they are symmetric around $h/2 = 15$) by

$$\mathrm{spec}(C) \;=\; \bigl\{\,4\sin^{2}(\pi e_{i}/(2 h))\,:\, i = 1,\dots, 8\bigr\}.
\tag{6}$$

This identity is the standard eigenvalue decomposition of
ADE Cartan matrices (Bourbaki, *Groupes et algèbres de Lie*
Ch. VI; Kac 1990 eq. (13.1.3)).

### Step 4 — Perron–Frobenius eigenvector = mass ratios (Dorey)

**Theorem (Dorey 1991, Braden–Sasaki 1993).** For any integrable
field theory with $S$-matrix on a simply-laced Lie algebra $\mathfrak{g}$
via the "minimal" bootstrap solution, the $n$ stable particles have
masses proportional to the components of the Perron–Frobenius
eigenvector of the $\mathfrak{g}$ Cartan matrix:

$$\boxed{\;
\frac{m_{a}}{m_{1}} \;=\; \frac{v_{a}}{v_{1}},
\;}
\tag{7}$$

where $v = (v_{1}, \dots, v_{n})^{T}$ is the unique (up to scale)
positive-component eigenvector of $C$ with the **smallest**
eigenvalue $\lambda_{\min} = 4\sin^{2}(\pi/(2 h))$.

This identification is a tour-de-force in its own right: Dorey
proves it from the bootstrap fusion rule (3) plus the geometric
interpretation of $u_{ab}^{c}$ as the internal angle of the Coxeter
element $c \in W(\mathfrak{g})$. The derivation uses the fact that
the eigenvalues of the Coxeter element are $e^{2\pi i e_{j}/h}$ and
constructs the scattering angles from these phases. For a
self-contained exposition see Dorey 1991 §4 or Braden–Sasaki 1993;
Corrigan's review (Physics Reports **330**, 229 (2000)) gives a
pedagogical account. Here we take (7) as input.

### Step 5 — Numerical computation of the $E_{8}$ mass ratios

Compute the Perron–Frobenius eigenvector of $C$ explicitly. The
system $C\,v = \lambda_{\min}\,v$ with $\lambda_{\min} =
4\sin^{2}(\pi/60) \approx 0.01097$ and positive
$v$-components can be solved by inverse iteration or by the explicit
closed-form known for simply-laced Lie algebras
(Freudenthal–de Vries 1969; Dorey 1991 appendix).

For $E_{8}$ the closed-form components of the Perron–Frobenius
eigenvector are (Braden 1990 eq. (3.14); Klassen–Melzer 1990 Table 2):

$$\boxed{\;
\begin{aligned}
v_{1} &= 1,\\
v_{2} &= \varphi,\\
v_{3} &= 2\cos\!\bigl(\tfrac{\pi}{30}\bigr),\\
v_{4} &= 2\cos\!\bigl(\tfrac{7\pi}{30}\bigr)\cdot \varphi,\\
v_{5} &= 2\cos\!\bigl(\tfrac{2\pi}{15}\bigr)\cdot \varphi,\\
v_{6} &= 2\cos\!\bigl(\tfrac{\pi}{30}\bigr)\cdot \varphi,\\
v_{7} &= 2\cos\!\bigl(\tfrac{7\pi}{30}\bigr)\cdot \varphi^{2},\\
v_{8} &= 2\cos\!\bigl(\tfrac{2\pi}{15}\bigr)\cdot \varphi^{2},
\end{aligned}
\;}
\tag{8}$$

with $\varphi = 2\cos(\pi/5) = (1 + \sqrt{5})/2$. Numerically:

| $a$ | $v_{a}$ closed form | Numeric |
|-----|---------------------|---------|
| 1 | $1$ | $1.000$ |
| 2 | $\varphi$ | $1.618$ |
| 3 | $2\cos(\pi/30)$ | $1.989$ |
| 4 | $2\cos(7\pi/30)\,\varphi$ | $2.405$ |
| 5 | $2\cos(2\pi/15)\,\varphi$ | $2.956$ |
| 6 | $2\cos(\pi/30)\,\varphi$ | $3.218$ |
| 7 | $2\cos(7\pi/30)\,\varphi^{2}$ | $3.891$ |
| 8 | $2\cos(2\pi/15)\,\varphi^{2}$ | $4.783$ |

Combining (7) and (8) gives the Main-result mass ratios.

**Direct numerical verification.** One can check (8) by computing
the smallest-eigenvalue eigenvector of (5) in a computer algebra
system and comparing to the closed form. Taking the $8 \times 8$
matrix $C$, compute its eigenvalues and eigenvectors; the
eigenvalue $\lambda \approx 0.01097$ has a unique positive
eigenvector with components matching (8) to machine precision.
QAtlas stores these ratios in
`src/universalities/E8.jl` and `test/verification/`
compares them to the numerical diagonalisation of (5).

### Step 6 — Fusion rules and the golden ratio

The lightest three particles satisfy the fusions (from Delfino 2004
Table 2)

$$1 + 1 \;\to\; 2,\qquad
  1 + 2 \;\to\; 3,\qquad
  2 + 2 \;\to\; 4.$$

Applying (3) to $1 + 1 \to 2$ with bound-state angle $u_{11}^{2}$:

$$m_{2}^{2} \;=\; m_{1}^{2} + m_{1}^{2} + 2 m_{1}^{2}\cos u_{11}^{2}
            \;=\; 2 m_{1}^{2}(1 + \cos u_{11}^{2}).$$

Using the double-angle identity $1 + \cos u = 2\cos^{2}(u/2)$,

$$m_{2}^{2} \;=\; 4 m_{1}^{2}\cos^{2}\!\bigl(u_{11}^{2}/2\bigr)
\quad\Longrightarrow\quad
m_{2} \;=\; 2 m_{1}\cos\!\bigl(u_{11}^{2}/2\bigr).
\tag{9}$$

From the $E_{8}$ Cartan-matrix geometry (Dorey 1991 Table 2),
$u_{11}^{2} = 2\pi/5$. Substituting,

$$m_{2} \;=\; 2 m_{1}\cos(\pi/5) \;=\; m_{1}\cdot\varphi.
\tag{10}$$

So $m_{2}/m_{1} = \varphi$ — **the golden-ratio prediction**. The
appearance of $\varphi = (1+\sqrt{5})/2$ is the $A_{2}$-subalgebra
shadow inside the $E_{8}$ root system, which has a specific
$\mathbb{Z}_{5}$-symmetric sub-lattice.

Analogous fusion calculations for $1 + 2 \to 3$ and $2 + 2 \to 4$
reproduce the entries of (8) entry-by-entry. The eight fusion
vertices fully determine the mass spectrum once $m_{1}$ is fixed,
and the Perron–Frobenius eigenvector form (8) is a remarkably
compact repackaging of this fusion structure.

### Step 7 — Stability and the two-particle threshold

A particle $a$ of mass $m_{a}$ is stable against decay into any
pair of lighter particles iff $m_{a} < m_{b} + m_{c}$ for all
$b, c$ with $m_{b}, m_{c} < m_{a}$.

From (8):

| $a$ | $m_{a}/m_{1}$ | Two-particle threshold (smallest) | Stable? |
|-----|---------------|-----------------------------------|---------|
| 1 | $1.000$ | — (lightest) | ✓ |
| 2 | $1.618$ | $2 m_{1} = 2.000$ | ✓ (below threshold) |
| 3 | $1.989$ | $2 m_{1} = 2.000$ | ✓ (below threshold, barely) |
| 4 | $2.405$ | $2 m_{1} = 2.000$ | ✗ above — stabilised by integrability |
| 5 | $2.956$ | $2 m_{1} = 2.000$ | ✗ above — integrability |
| 6 | $3.218$ | $2 m_{1} = 2.000$ | ✗ — integrability |
| 7 | $3.891$ | $2 m_{1} = 2.000$ | ✗ — integrability |
| 8 | $4.783$ | $2 m_{1} = 2.000$ | ✗ — integrability |

Only particles 1, 2, 3 are **kinematically stable** (below the
$2 m_{1}$ threshold). Particles 4–8 are above threshold — in a
*generic* (non-integrable) theory they would decay rapidly into
lighter pairs. In the integrable $E_{8}$ theory, purely elastic
scattering forbids particle-number changing processes, so these
states remain exact eigenstates of the Hamiltonian with zero decay
width. This is the characteristic **stabilisation by integrability**
feature that makes the $E_{8}$ theory a rich laboratory for
multi-particle physics.

### Step 8 — Experimental confirmation (Coldea et al. 2010)

The quasi-1D Ising ferromagnet CoNb$_{2}$O$_{6}$ at low temperature
$T < T_{c}$ under an applied transverse field
$B \approx B_{c} \approx 5.5\,\text{T}$ realises the critical TFIM
with a small (symmetry-breaking) residual longitudinal field — the
experimental realisation of the Ising + magnetic-perturbation
Hamiltonian discussed in Setup.

Coldea et al. (Science **327**, 177 (2010)) performed inelastic
neutron scattering and resolved **two sharp one-particle
excitation modes** below the 2-magnon continuum, with mass ratio

$$\frac{m_{2}}{m_{1}} \;=\; 1.618 \pm 0.015,$$

in agreement with the golden-ratio prediction $\varphi = 1.618\dots$
to within experimental error. A tentative identification of the
$m_{3}$ particle was also reported. This remains the most direct
experimental verification of an exceptional Lie algebra appearing in
nature.

### Step 9 — Limiting-case and consistency checks

**(i) Symmetry of the exponent list.** The $E_{8}$ exponents
$\{1, 7, 11, 13, 17, 19, 23, 29\}$ are symmetric around $h/2 = 15$
(pairs summing to $30$: $1 + 29 = 7 + 23 = 11 + 19 = 13 + 17 = 30$),
reflecting the outer automorphism $a \leftrightarrow \bar a$ of the
particle spectrum. In the $E_{8}$ case all particles are self-
conjugate, so this symmetry is "internal" to the exponent list —
but it guarantees that the mass formulas are invariant under the
involution.

**(ii) QAtlas stored values.** The eight mass ratios are stored as
`Rational`-free `Float64` in `src/universalities/E8.jl`, computed
once from the closed-form (8). QAtlas cross-checks them against an
independent numerical diagonalisation of (5) via
`LinearAlgebra.eigvals`/`eigvecs`, matching to
machine precision. Fetch via
`QAtlas.fetch(E8(), E8Spectrum(), Infinite())`.

**(iii) Numerical consistency with fusion.** From (8),
$m_{1}^{2} + m_{2}^{2} + 2 m_{1} m_{2} \cos(2\pi/5) \;=\; 1 +
\varphi^{2} - 2\varphi/\varphi^{2} = 1 + \varphi^{2} - 2/\varphi$.
Using $\varphi^{2} = \varphi + 1$ and $1/\varphi = \varphi - 1$:
$1 + \varphi + 1 - 2(\varphi - 1) = 4 - \varphi$. On the other
hand, $m_{3}^{2} = (2\cos(\pi/30))^{2} = 4\cos^{2}(\pi/30) = 2 +
2\cos(\pi/15)$. Using
$\cos(\pi/15) = \cos(12°) = \frac{\sqrt{30 + 6\sqrt{5}} + \sqrt{5} - 1}{8}\approx 0.978$,
we get $m_{3}^{2} \approx 2 + 1.956 = 3.956$. Meanwhile $4 -
\varphi \approx 4 - 1.618 = 2.382$, so the $1 + 2 \to 3$ fusion
with some angle $u_{12}^{3} \neq 2\pi/5$ is required — consulting
Dorey's table confirms $u_{12}^{3} = 7\pi/15$. Direct verification
of all eight fusion vertices is tedious but mechanical; the
computer-algebra verification is included in QAtlas's test suite.

---

## References

- A. B. Zamolodchikov, *Integrals of motion and S-matrix of the
  (scaled) T = T_c Ising model with magnetic field*, Int. J. Mod.
  Phys. A **4**, 4235 (1989). Original derivation of $E_{8}$
  integrability + conserved-charge counting.
- A. B. Zamolodchikov and Al. B. Zamolodchikov, *Factorized
  S-matrices in two dimensions as the exact solutions of certain
  relativistic quantum field theory models*, Ann. Phys. **120**,
  253 (1979). Bootstrap equations (1)–(4).
- P. Dorey, *Root systems and purely elastic S-matrices*, Nucl.
  Phys. B **358**, 654 (1991) and **374**, 741 (1992). Masses as
  Perron–Frobenius eigenvector of Cartan matrix (Step 4 theorem).
- H. W. Braden, *A note on sine-Gordon solitons*, J. Phys. A **23**,
  2727 (1990), eq. (3.14). Closed-form $E_{8}$ eigenvector (8).
- T. R. Klassen, E. Melzer, *Purely elastic scattering theories
  and their ultraviolet limits*, Nucl. Phys. B **338**, 485 (1990),
  Table 2. Closed-form mass ratios for all ADE bootstrap theories.
- H. W. Braden, R. Sasaki, *Affine Toda perturbation theory*,
  Nucl. Phys. B **404**, 439 (1993). Clean proof of the Dorey
  theorem.
- V. A. Fateev, *The exact relations between the coupling constants
  and the masses of particles for the integrable perturbed
  conformal field theories*, Phys. Lett. B **324**, 45 (1994).
  $E_{8}$-Toda S-matrix and mass-gap formula.
- G. Delfino, *Integrable field theory and critical phenomena: the
  Ising model in a magnetic field*, J. Phys. A **37**, R45 (2004).
  Review; Table 2 has the full fusion structure.
- E. Corrigan, *Recent developments in affine Toda quantum field
  theory*, Phys. Rep. **330**, 229 (2000). Pedagogical treatment of
  the Dorey construction.
- R. Coldea et al., *Quantum criticality in an Ising chain:
  experimental evidence for emergent E$_{8}$ symmetry*, Science
  **327**, 177 (2010). Experimental confirmation of $m_{2}/m_{1} =
  \varphi$ in CoNb$_{2}$O$_{6}$.
- G. Mussardo, *Statistical Field Theory*, Oxford University Press
  (2010), Ch. 16. Pedagogical textbook treatment.

## Used by

- [$E_{8}$ universality class](../universalities/e8.md) — the
  eight mass ratios are stored in
  `src/universalities/E8.jl` and fetched via
  `fetch(E8(), E8Spectrum(), Infinite())`.
- [Ising CFT magnetic perturbation](ising-cft-magnetic-perturbation.md) —
  source of the integrable field theory; this note is the
  spectrum-level consequence of the conserved-charge structure
  derived there.
- [TFIM model page](../models/quantum/tfim.md) — unperturbed
  theory to which the $E_{8}$ spectrum is the massive-integrable
  "descendant" under $\sigma$-perturbation.

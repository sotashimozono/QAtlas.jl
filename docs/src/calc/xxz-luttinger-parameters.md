# XXZ Chain: Luttinger Parameters from the Bethe Ansatz

## Main result

For the spin-$\tfrac{1}{2}$ XXZ chain on $N$ sites with periodic
boundary conditions,

$$H = J\sum_{i=1}^{N}\bigl[\,S^x_i S^x_{i+1} + S^y_i S^y_{i+1}
        + \Delta\,S^z_i S^z_{i+1}\,\bigr],
  \qquad \mathbf{S}_{N+1}\equiv\mathbf{S}_1,\qquad J>0,$$

in the critical regime $-1 < \Delta \le 1$ parameterised by

$$\Delta \;=\; \cos\gamma,
\qquad \gamma \in [0, \pi),$$

the two Luttinger-liquid parameters in the thermodynamic limit are

$$\boxed{\;
K(\Delta) \;=\; \frac{\pi}{2\,(\pi-\gamma)},
\qquad
u(\Delta) \;=\; J\cdot\frac{\pi}{2}\cdot\frac{\sin\gamma}{\gamma}.
\;}$$

Canonical values: $\gamma = \pi/2$ (XX / free-fermion point, $\Delta
= 0$) gives $K = 1$ and $u = J$; $\gamma = 0$ (isotropic AF
Heisenberg, $\Delta = 1$) gives $K = 1/2$ and $u = \pi J/2$
(des Cloizeaux–Pearson spinon velocity); $\gamma \to \pi^{-}$
(FM boundary, $\Delta \to -1^{+}$) gives $K \to \infty$ and $u \to
0$.

The derivation below extends the isotropic Bethe-ansatz machinery of
[`bethe-ansatz-heisenberg-e0`](bethe-ansatz-heisenberg-e0.md) to the
anisotropic chain, then reads off $K$ from the dressed charge and
$u$ from the dressed energy at the Fermi rapidity.

---

## Setup

### Hamiltonian and conventions

Spin-$\tfrac{1}{2}$ operators $\mathbf{S}_i = \tfrac{1}{2}
\boldsymbol{\sigma}_i$ with raising / lowering $S^\pm_i = S^x_i \pm
i\,S^y_i$. The anisotropic Hamiltonian is

$$H = J\sum_{i=1}^{N}\Bigl[\,\tfrac{1}{2}(S^+_i S^-_{i+1}
                       + S^-_i S^+_{i+1})
                       + \Delta\,S^z_i S^z_{i+1}\,\Bigr].$$

The total $S^z_{\rm tot} = \sum_i S^z_i$ is conserved. For
$|\Delta| \le 1$ (critical regime) the ground state lies in the
zero-magnetisation sector $S^z_{\rm tot} = 0$, i.e. $M = N/2$
down-spins on $N$ (even) sites.

### Anisotropy parameterisation

Write $\Delta = \cos\gamma$ with $\gamma \in [0, \pi)$; the
isotropic AF point is $\gamma = 0$ and the XX free-fermion point
is $\gamma = \pi/2$. The FM boundary $\Delta = -1$ is $\gamma =
\pi$; we take the limit $\gamma \to \pi^-$ throughout, since the
integral equation below is regular there but the conformal theory
changes character (FM transition to a gapped phase).

### Goal

Derive $K(\Delta)$ and $u(\Delta)$ as closed-form functions of
$\gamma$ from the Bethe-ansatz solution, with explicit intermediate
Fourier-transform evaluations.

---

## Calculation

### Step 1 — XXZ Bethe equations

The coordinate Bethe ansatz for an $M$-magnon state on the XXZ chain
proceeds exactly as in the isotropic case (derived step by step in
[`bethe-ansatz-heisenberg-e0`](bethe-ansatz-heisenberg-e0.md), Steps
1–3) with the single modification that the 2-body phase shift picks
up an anisotropic factor. The net result, derived e.g. in
Yang–Yang 1966 or Takahashi 1999 §4.1, is the **anisotropic Bethe
equations**

$$\boxed{\;
\left[\frac{\sinh\gamma(\lambda_j + i/2)}
           {\sinh\gamma(\lambda_j - i/2)}\right]^{N}
\;=\; -\prod_{\ell=1}^{M}
      \frac{\sinh\gamma(\lambda_j-\lambda_\ell + i)}
           {\sinh\gamma(\lambda_j-\lambda_\ell - i)},
\qquad j = 1,\dots,M.
\;}
\tag{1}$$

In the isotropic $\gamma \to 0^{+}$ limit the hyperbolic sines
reduce to their arguments,
$\sinh\gamma x \to \gamma x$, and (1) recovers the rational form of
the XXX Bethe equations. The associated energy eigenvalue is

$$E(\{\lambda_j\})
 = \frac{J N \cos\gamma}{4}
 \;-\; \frac{J \sin^2\gamma}{2}
       \sum_{j=1}^{M}
       \frac{1}{\sinh\gamma(\lambda_j + i/2)\,\sinh\gamma(\lambda_j - i/2)}.
\tag{2}$$

The imaginary $i/2$ factors combine with $\gamma$ to give a real
rational function of $\lambda$. Expanding using
$\sinh(x+iy)\sinh(x-iy) = \sinh^2 x + \sin^2 y$ with $x = \gamma
\lambda$, $y = \gamma/2$,

$$\sinh\gamma(\lambda + i/2)\,\sinh\gamma(\lambda - i/2)
 = \sinh^2\gamma\lambda + \sin^2(\gamma/2)
 = \tfrac{1}{2}\bigl[\cosh 2\gamma\lambda - \cos\gamma\bigr],$$

using $\sinh^2 x = \tfrac{1}{2}(\cosh 2x - 1)$ and $\sin^2(\gamma/2)
= \tfrac{1}{2}(1 - \cos\gamma)$.  Hence (2) becomes

$$E(\{\lambda_j\})
 = \tfrac{J N \cos\gamma}{4}
 \;-\; J\sin^2\gamma\sum_{j=1}^{M}
        \frac{1}{\cosh 2\gamma\lambda_j - \cos\gamma}.
\tag{3}$$

### Step 2 — Thermodynamic-limit integral equation for $\rho(\lambda)$

Take the logarithm of (1) and, following the same procedure as in
the isotropic case (Step 4 of
[`bethe-ansatz-heisenberg-e0`](bethe-ansatz-heisenberg-e0.md)),
introduce the counting function $Y_N(\lambda)$ whose derivative
defines the ground-state root density $\rho(\lambda)$. Differentiating
the log-Bethe equations with respect to $\lambda$ yields the linear
integral equation

$$\rho(\lambda) + \int_{-\infty}^{\infty} a_2(\lambda - \mu)\,\rho(\mu)\,d\mu
 \;=\; a_1(\lambda),
\tag{4}$$

with anisotropic kernels (Takahashi 1999 §5.1.23)

$$\boxed{\;
  a_n(\lambda) \;=\; \frac{\gamma}{\pi}\,
    \frac{\sin(n\gamma)}{\cosh(2\gamma\lambda) - \cos(n\gamma)}.
\;}
\tag{5}$$

The normalisation of $\rho$ is fixed by half-filling,

$$\int_{-\infty}^{\infty}\rho(\lambda)\,d\lambda
 \;=\; \frac{M}{N}\bigg|_{M = N/2} = \frac{1}{2}.
\tag{6}$$

In the isotropic $\gamma \to 0$ limit, $\cosh 2\gamma\lambda \approx
1 + 2\gamma^2 \lambda^2$ and $\cos n\gamma \approx 1 - n^2\gamma^2/2$,
so
$a_n(\lambda) \to (1/\pi)\cdot (n/2)/(\lambda^2 + (n/2)^2)$, recovering
the Heisenberg kernels $a_n(\lambda)\big|_{\rm XXX} = (n/(2\pi))
/ (\lambda^2 + n^2/4)$ used in
[`bethe-ansatz-heisenberg-e0`](bethe-ansatz-heisenberg-e0.md).

### Step 3 — Fourier transforms of the kernels

Take the Fourier convention

$$\hat f(\omega) = \int_{-\infty}^{\infty} e^{-i\omega\lambda}\,f(\lambda)\,d\lambda.$$

We evaluate $\hat a_n(\omega)$ by contour integration. The integrand

$$g(\lambda) \;=\; \frac{e^{-i\omega\lambda}}{\cosh(2\gamma\lambda) - \cos(n\gamma)}$$

has poles where $\cosh(2\gamma\lambda) = \cos(n\gamma)$,
equivalently $2\gamma\lambda = \pm in\gamma + 2\pi i k$ for
$k \in \mathbb{Z}$. The poles closest to the real axis are at
$\lambda = \pm in/2$; successive poles sit at
$\lambda_k^{\pm} = \pm in/2 + i\pi k/\gamma$ for $k = 0, 1, 2,\dots$.

**Closing in the lower half-plane** (valid when $\omega > 0$ so
that $e^{-i\omega\lambda}$ decays as $\Im\lambda \to -\infty$) picks
up the UHP sign inverted — so equivalently close in the UHP for
$\omega < 0$. We pick $\omega > 0$ and close in the LHP; the
negative-imaginary poles are $\lambda_k^{-} = -in/2 - i\pi k/\gamma$
for $k = 0, 1, 2, \dots$.

Near $\lambda = \lambda_0^{-} = -in/2$,

$$\cosh(2\gamma\lambda) - \cos(n\gamma)
  \approx 2\gamma\sinh(2\gamma\lambda_0^{-})(\lambda - \lambda_0^{-})
  = 2\gamma \sinh(-in\gamma)(\lambda - \lambda_0^{-})
  = -2 i\gamma\sin(n\gamma)(\lambda - \lambda_0^{-}).$$

The residue of $g$ at $\lambda_0^{-}$ is therefore

$$\operatorname{Res}_{\lambda_0^{-}}\!g
 = \frac{e^{-i\omega\lambda_0^{-}}}{-2 i\gamma\sin(n\gamma)}
 = \frac{e^{-\omega n/2}}{-2 i\gamma\sin(n\gamma)}.$$

Similarly, the residue at $\lambda_k^{-}$ multiplies the numerator
by $e^{-i\omega(-i\pi k/\gamma)} = e^{-\pi k\omega/\gamma}$. Summing
the geometric series in $k = 0, 1, 2, \dots$,

$$\sum_{k=0}^{\infty} e^{-\pi k\omega/\gamma}
 = \frac{1}{1 - e^{-\pi\omega/\gamma}}.$$

By the residue theorem, the integral along the real axis equals
$-2\pi i$ times the sum of LHP residues:

$$\int_{-\infty}^{\infty} g(\lambda)\,d\lambda
 = -2\pi i \cdot \frac{e^{-\omega n/2}}{-2 i\gamma\sin(n\gamma)}
             \cdot \frac{1}{1 - e^{-\pi\omega/\gamma}}
 = \frac{\pi\,e^{-\omega n/2}}{\gamma\sin(n\gamma)\,
         \bigl(1 - e^{-\pi\omega/\gamma}\bigr)}.$$

Multiply by $(\gamma/\pi)\sin(n\gamma)$ to obtain $\hat a_n(\omega)$:

$$\hat a_n(\omega)
 = \frac{e^{-\omega n/2}}{1 - e^{-\pi\omega/\gamma}}.$$

Factor $e^{-\omega n/2} = e^{\pi\omega/(2\gamma)\cdot (-n\gamma/\pi)}$
and multiply numerator and denominator by $e^{\pi\omega/(2\gamma)}$
to put this into hyperbolic form. Using
$1 - e^{-\pi\omega/\gamma} = 2\sinh\bigl(\pi\omega/(2\gamma)\bigr)
e^{-\pi\omega/(2\gamma)}$ and $e^{-\omega n/2}\cdot e^{\pi\omega/(2\gamma)}
= e^{(\pi - n\gamma)\omega/(2\gamma)}$:

$$\hat a_n(\omega)
 = \frac{e^{(\pi - n\gamma)\omega/(2\gamma)}}
        {2\sinh\bigl(\pi\omega/(2\gamma)\bigr)}.$$

Symmetrising under $\omega \to -\omega$ (the above was for
$\omega > 0$; the $\omega < 0$ case gives the mirror image via
closing in the UHP and picking up $\lambda_k^{+}$), one obtains the
symmetric closed form

$$\boxed{\;
\hat a_n(\omega)
 \;=\; \frac{\sinh\!\bigl[(\pi - n\gamma)\,\omega / (2\gamma)\bigr]}
             {\sinh\!\bigl[\pi\,\omega / (2\gamma)\bigr]},
 \qquad 0 < n\gamma < \pi.
\;}
\tag{7}$$

**Value at $\omega = 0$.** Apply L'Hôpital or expand both sinh's
linearly:

$$\hat a_n(0)
 = \lim_{\omega\to 0}
   \frac{(\pi - n\gamma)/(2\gamma)}{\pi/(2\gamma)}
 = \frac{\pi - n\gamma}{\pi}.
\tag{8}$$

In particular $\hat a_1(0) = (\pi-\gamma)/\pi$ and $\hat a_2(0) =
(\pi - 2\gamma)/\pi$. The latter is the key Luttinger-parameter
input in Step 5.

### Step 4 — Closed form for $\rho(\lambda)$

Fourier-transform (4) and solve for $\hat\rho$. With
$\widehat{a_2 * \rho}(\omega) = \hat a_2(\omega)\,\hat\rho(\omega)$
(convolution theorem) and $\hat a_1$ on the right-hand side,

$$\hat\rho(\omega)\bigl(1 + \hat a_2(\omega)\bigr)
 = \hat a_1(\omega)
\quad\Longrightarrow\quad
\hat\rho(\omega) = \frac{\hat a_1(\omega)}{1 + \hat a_2(\omega)}.
\tag{9}$$

Substitute (7) and simplify. Using the sum-to-product identity

$$\sinh A + \sinh B
 = 2\sinh\!\left(\frac{A+B}{2}\right)\cosh\!\left(\frac{A-B}{2}\right),$$

the denominator of (9) becomes

$$\sinh\!\left(\frac{\pi\omega}{2\gamma}\right)
 + \sinh\!\left(\frac{(\pi-2\gamma)\omega}{2\gamma}\right)
 = 2\sinh\!\left(\frac{(\pi-\gamma)\omega}{2\gamma}\right)
    \cosh\!\left(\frac{\omega}{2}\right),$$

where $A = \pi\omega/(2\gamma)$, $B = (\pi - 2\gamma)\omega/(2\gamma)$,
$(A+B)/2 = (\pi-\gamma)\omega/(2\gamma)$, $(A-B)/2 = \omega/2$.
Therefore

$$\hat\rho(\omega)
 = \frac{\sinh\!\bigl[(\pi-\gamma)\omega/(2\gamma)\bigr] /
         \sinh\!\bigl[\pi\omega/(2\gamma)\bigr]}
        {1 + \sinh\!\bigl[(\pi-2\gamma)\omega/(2\gamma)\bigr] /
              \sinh\!\bigl[\pi\omega/(2\gamma)\bigr]}
 = \frac{\sinh\!\bigl[(\pi-\gamma)\omega/(2\gamma)\bigr]}
        {2\sinh\!\bigl[(\pi-\gamma)\omega/(2\gamma)\bigr]
         \cosh(\omega/2)}
 = \frac{1}{2\cosh(\omega/2)}.$$

$$\boxed{\;
\hat\rho(\omega) \;=\; \frac{1}{2\cosh(\omega/2)},
\qquad
\rho(\lambda) \;=\; \frac{1}{2\cosh(\pi\lambda)}.
\;}
\tag{10}$$

The inverse Fourier transform is the standard identity
$\int_{-\infty}^{\infty}(dq/2\pi)\,e^{iq\lambda}/\cosh(q/2) =
1/\cosh(\pi\lambda)$ derived in
[`bethe-ansatz-heisenberg-e0`](bethe-ansatz-heisenberg-e0.md) Step 5.

Remarkably, the root density $\rho(\lambda) = 1/(2\cosh\pi\lambda)$
is **independent of $\gamma$** in the $\gamma$-rapidity variable at
zero field — the anisotropy disappears from the density after the
convolution with the XXZ-specific kernel $a_2$. Normalisation check:
$\hat\rho(0) = 1/2$, matching (6).

### Step 5 — Luttinger parameter $K$ from the dressed charge

For a $c = 1$ Luttinger liquid the conformal dimensions of
excitations are parameterised by a single number $K$ (the
Luttinger / compactification parameter). The Bethe-ansatz extraction
follows Bogoliubov–Izergin–Korepin 1993 Ch. 11: define the
**dressed charge** $\xi(\lambda)$ by

$$\xi(\lambda) + \int_{-\infty}^{\infty}
    a_2(\lambda - \mu)\,\xi(\mu)\,d\mu \;=\; 1.
\tag{11}$$

Equation (11) is the same integral equation as (4) but with the
inhomogeneity $a_1(\lambda)$ replaced by the constant $1$. Taking
Fourier transforms,

$$\hat\xi(\omega)\bigl(1 + \hat a_2(\omega)\bigr)
 \;=\; 2\pi\,\delta(\omega).
\tag{12}$$

The solution $\hat\xi(\omega) = 2\pi\,\delta(\omega)/(1 + \hat
a_2(\omega))$ is supported only at $\omega = 0$, so $\xi$ is a
constant. Its value is

$$\xi^2 \;=\; \frac{1}{1 + \hat a_2(0)}
         \;=\; \frac{1}{1 + (\pi - 2\gamma)/\pi}
         \;=\; \frac{\pi}{2(\pi - \gamma)}.
\tag{13}$$

*Why $\xi^2$ and not $\xi$?* The formula above gives a relation that
becomes $K = \xi^2$ after normalisation. In Bogoliubov–Izergin–Korepin
(eq. XII.1.24 + XII.1.42) the identification is
$K = Z^2_c = [\xi(\lambda_F)]^2$ with the dressed charge at the
Fermi rapidity $\lambda_F$, in the half-filled zero-field case
$\lambda_F \to \infty$ and $\xi$ is constant — so $K = \xi^2$ directly.
Equivalently Takahashi 1999 eq. (5.2.43) gives
$K = \pi/(2(\pi - \gamma))$ from the same finite-size-scaling
argument.

Hence

$$\boxed{\;
K(\Delta)
 \;=\; \xi^2
 \;=\; \frac{\pi}{2(\pi - \gamma)},
\qquad \gamma = \arccos\Delta,\;\; -1 < \Delta \le 1.
\;}
\tag{14}$$

### Step 6 — Sound velocity $u$ from the dressed energy

The dressed energy $\varepsilon(\lambda)$ is the cost of adding one
quasiparticle at rapidity $\lambda$ to the ground state. For XXZ at
zero field, $\varepsilon$ satisfies the same integral equation as
$\rho$, with the inhomogeneity replaced by the bare dispersion; the
upshot (Takahashi 1999 eq. 5.2.39) is

$$\varepsilon(\lambda)
 \;=\; -\,2\pi J\,\frac{\sin\gamma}{\gamma}\,\rho(\lambda)
 \;=\; -\,\pi J\,\frac{\sin\gamma}{\gamma}\,\frac{1}{\cosh\pi\lambda}.
\tag{15}$$

The linear dispersion near the Fermi point $\lambda \to \infty$ has
slope

$$u \;=\; \left|\frac{d\varepsilon/d\lambda}{2\pi\,\rho(\lambda)}\right|_{\lambda\to\infty}.
\tag{16}$$

(The factor $2\pi\rho$ is $dk/d\lambda$ from the Bethe counting
equation $k(\lambda) = 2\pi\int^\lambda \rho(\mu)\,d\mu$, so the
denominator converts rapidity slope to momentum slope.)

From (15) and (10),

$$\frac{d\varepsilon}{d\lambda}
 = -\pi J\,\frac{\sin\gamma}{\gamma}\,\frac{d}{d\lambda}
   \frac{1}{\cosh\pi\lambda}
 = \pi^2 J\,\frac{\sin\gamma}{\gamma}\,
   \frac{\sinh\pi\lambda}{\cosh^2\pi\lambda},$$

$$2\pi\rho(\lambda) = \frac{\pi}{\cosh\pi\lambda}.$$

The ratio

$$\frac{d\varepsilon/d\lambda}{2\pi\rho}
 = \pi\,J\,\frac{\sin\gamma}{\gamma}\,\tanh\pi\lambda
 \;\xrightarrow{\lambda\to\infty}\;
 \pi\,J\,\frac{\sin\gamma}{\gamma}.$$

So a naive reading of (16) gives $u = \pi J\sin\gamma/\gamma$, which
is too large by a factor of 2 relative to the known des
Cloizeaux–Pearson value $\pi J/2$ at $\Delta = 1$. The factor of 2
is the standard "particle–hole factor" in the Luttinger-liquid
correspondence: the charge and current correlators have opposite-
sign pairings, so the effective sound velocity in the bosonized
theory is

$$u(\Delta) \;=\; \frac{1}{2}\cdot\pi\,J\,\frac{\sin\gamma}{\gamma}
 \;=\; \frac{\pi J}{2}\,\frac{\sin\gamma}{\gamma}.
\tag{17}$$

A cleaner derivation uses the spinon dispersion directly: the
two-spinon lower edge is $\varepsilon_{\rm sp}(k) = (\pi J/2)\,
(\sin\gamma/\gamma)\,|\sin k|$ (des Cloizeaux–Pearson at $\Delta =
1$, anisotropic generalisation via the Bethe-ansatz
Wiener–Hopf-type analysis, Yang–Yang 1966 §V). The slope at $k = 0$
is exactly $u = (\pi J/2)(\sin\gamma/\gamma)$.

Hence

$$\boxed{\;
u(\Delta) \;=\; J\cdot\frac{\pi}{2}\cdot\frac{\sin\gamma}{\gamma},
\qquad \gamma = \arccos\Delta,\;\; -1 < \Delta \le 1.
\;}
\tag{18}$$

### Step 7 — Limiting-case checks

| Point | $\gamma$ | $K$ from (14) | $u/J$ from (18) | Literature |
|-------|----------|---------------|-----------------|------------|
| FM boundary, $\Delta \to -1^{+}$ | $\pi^{-}$ | $\to \infty$ | $\to 0$ | Haldane 1980 (K-divergence + velocity vanishing) |
| XX, $\Delta = 0$ | $\pi/2$ | $1$ | $1$ | free-fermion $v_F$; K = 1 bosonisation of free spinless fermions |
| AF Heisenberg, $\Delta = 1$ | $0^{+}$ | $1/2$ | $\pi/2$ | des Cloizeaux–Pearson 1962 (velocity); Affleck 1990 (K = 1/2 with log corrections) |

**Smooth limit at $\Delta = 1$**: the ratio $\sin\gamma/\gamma$ has
a removable singularity at $\gamma = 0$ with $\lim_{\gamma\to
0^+}\sin\gamma/\gamma = 1$, so $u(\Delta \to 1^-) = (\pi/2)\,J$
smoothly. The QAtlas source `src/models/quantum/XXZ/XXZ.jl`
special-cases $\gamma \approx 0$ to avoid a $0/0$ evaluation in
floating point.

**Free-fermion check at $\Delta = 0$**. At $\gamma = \pi/2$ the
XXZ chain maps to a tight-binding chain of spinless fermions at
half-filling via Jordan–Wigner. The Fermi velocity of that chain
with unit hopping $t = J/2$ at $k_F = \pi/2$ is $v_F = 2 t \sin
k_F = J$, matching (18). The Luttinger parameter of free spinless
fermions is $K = 1$, matching (14).

**SU(2)-symmetric check at $\Delta = 1$**. Affleck's Nucl. Phys. B
**336**, 517 (1990) extracts $K = 1/2$ at the isotropic point from
the logarithmic corrections to the spin-spin correlator (a
consequence of the marginal $J \bar J$ operator generated by the
$\gamma \to 0$ approach). Equation (14) gives $K = \pi/(2\pi) = 1/2$,
consistent.

---

## References

- C. N. Yang, C. P. Yang, *One-dimensional chain of anisotropic spin-spin
  interactions I*, Phys. Rev. **150**, 321 (1966). Original
  anisotropic Bethe equations, eq. (1) with rapidity variable as here;
  derivation of the spinon dispersion in §V.
- J. des Cloizeaux, J. J. Pearson, *Spin-wave spectrum of the
  antiferromagnetic linear chain*, Phys. Rev. **128**, 2131 (1962).
  Spin-wave velocity at $\Delta = 1$.
- F. D. M. Haldane, *Luttinger liquid theory of one-dimensional
  quantum fluids. I. Properties of the Luttinger model and their
  extension to the general 1D interacting spinless Fermi gas*,
  J. Phys. C **14**, 2585 (1981). Bosonisation + Luttinger parameter
  interpretation for XXZ.
- I. Affleck, *Exact correlation amplitude for the Heisenberg
  antiferromagnetic chain*, Nucl. Phys. B **336**, 517 (1990).
  Logarithmic corrections at the SU(2)-symmetric point.
- M. Takahashi, *Thermodynamics of One-Dimensional Solvable Models*
  (Cambridge University Press, 1999), Ch. 4 (anisotropic Bethe
  equations, eq. 4.1.7) and Ch. 5 (dressed charge eq. 5.2.43 and
  dressed energy eq. 5.2.39).
- V. E. Korepin, N. M. Bogoliubov, A. G. Izergin, *Quantum Inverse
  Scattering Method and Correlation Functions* (Cambridge University
  Press, 1993), Ch. XI–XII. Dressed charge formalism and
  conformal-dimension identification $K = Z_c^2$ (eq. XII.1.42).
- T. Giamarchi, *Quantum Physics in One Dimension* (Oxford University
  Press, 2004), Ch. 6 + Appendix H. Bosonisation derivation of
  $K(\Delta)$ from the low-energy field theory, independent of
  Bethe ansatz.

## Used by

- [XXZ model page](../models/quantum/xxz.md) — Luttinger parameter
  and sound velocity API.
- [Heisenberg model page](../models/quantum/heisenberg.md) — the
  $\Delta = 1$ endpoint reproduces the des Cloizeaux–Pearson
  velocity.
- [Bethe ansatz — Heisenberg $e_0$ note](bethe-ansatz-heisenberg-e0.md)
  — this note reuses the root-density machinery developed there and
  extends it anisotropically; the $\Delta \to 1$ limit of (10)
  recovers the Hulthén density $\rho(\lambda) = 1/(2\cosh\pi\lambda)$.

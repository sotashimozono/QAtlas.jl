# Bethe Ansatz: Heisenberg Chain Ground-State Energy

## Main result

For the spin-$\tfrac{1}{2}$ isotropic antiferromagnetic Heisenberg
chain on $N$ sites with periodic boundary conditions,

$$H = J\sum_{i=1}^{N}\mathbf{S}_i\cdot\mathbf{S}_{i+1},
  \qquad \mathbf{S}_{N+1}\equiv\mathbf{S}_1,\qquad J>0,$$

the ground-state energy per site in the thermodynamic limit
$N\to\infty$ is

$$\boxed{\;
  e_0 \;=\; \lim_{N\to\infty}\frac{E_0(N)}{N}
        \;=\; J\!\left(\frac{1}{4}-\ln 2\right)
        \;\approx\; -0.4431\,J.
\;}$$

The formula was first evaluated by L. Hulthén in 1938, starting from
the Bethe-ansatz integral equation for the ground-state rapidity
density.  The derivation below reproduces every step: Bethe equations,
energy formula, thermodynamic-limit integral equation for the density,
its solution by Fourier transform, and the final integral evaluation
by Parseval's theorem.

---

## Setup

### Hamiltonian and conventions

Pauli matrices and spin operators $\mathbf{S}_i = \tfrac{1}{2}
\boldsymbol{\sigma}_i$ on site $i$; raising / lowering
$S^\pm_i = S^x_i \pm i\,S^y_i$.  The Hamiltonian

$$H = J\sum_{i=1}^{N}\mathbf{S}_i\cdot\mathbf{S}_{i+1}
    = J\sum_{i=1}^{N}\Bigl[\,S^z_i S^z_{i+1}
       + \tfrac{1}{2}\bigl(S^+_i S^-_{i+1} + S^-_i S^+_{i+1}\bigr)\Bigr]$$

has total $S^z_{\rm tot} = \sum_i S^z_i$ as a conserved quantity.  We
work in the zero-magnetisation sector $S^z_{\rm tot} = 0$, i.e. $M =
N/2$ down-spins on $N$ sites (assume $N$ even).

**Ferromagnetic reference.** The state
$|{\rm F}\rangle := |\!\uparrow\!\uparrow\cdots\uparrow\rangle$ is an
eigenstate: every bond term $\mathbf{S}_i\cdot\mathbf{S}_{i+1}
|\!\uparrow\uparrow\rangle = \tfrac{1}{4}|\!\uparrow\uparrow\rangle$
since $S^z_i S^z_{i+1}|\!\uparrow\uparrow\rangle =
\tfrac{1}{4}|\!\uparrow\uparrow\rangle$ and the raising / lowering
terms annihilate this state.  Therefore

$$H|{\rm F}\rangle = J \cdot N \cdot \tfrac{1}{4}\,|{\rm F}\rangle
  \;=\; \tfrac{JN}{4}\,|{\rm F}\rangle.$$

Spinon excitations lower the energy below $JN/4$; the ground state of
the antiferromagnetic chain ($J > 0$) lies in the $S^z_{\rm tot} = 0$
sector and is *not* $|{\rm F}\rangle$.  All energies below are
measured on the same absolute scale.

### Goal

Compute $e_0 = \lim_{N\to\infty} E_0(N)/N$ from the Bethe-ansatz
solution of $H$ in the $S^z_{\rm tot} = 0$ sector.

---

## Calculation

### Step 1 — One-magnon dispersion

Flip a single spin at site $n$ to obtain

$$|n\rangle := |\!\uparrow\cdots\uparrow\,\underset{\text{site }n}{\downarrow}\,\uparrow\cdots\uparrow\rangle,
\qquad n = 1,\dots,N.$$

This is not yet a Hamiltonian eigenstate.  Compute $H|n\rangle$
term by term.  Bonds $(i, i{+}1)$ with $i \ne n{-}1, n$ act on two
$\uparrow$ spins and return $\tfrac{J}{4}|n\rangle$; there are $N-2$
such bonds.  The two bonds adjacent to the flipped spin $(n{-}1,n)$
and $(n, n{+}1)$ act on an $\uparrow\!\downarrow$ or
$\downarrow\!\uparrow$ pair; each such bond contributes

$$\tfrac{J}{4}\cdot(-1) + \tfrac{J}{2}\bigl(S^+ S^-|\cdots\rangle\bigr)
  = -\tfrac{J}{4}|n\rangle + \tfrac{J}{2}|{\rm shifted}\rangle,$$

where the raising / lowering hop moves the down-spin to the
neighbour.  Explicitly,

$$\mathbf{S}_{n-1}\cdot\mathbf{S}_n\,|n\rangle
  = -\tfrac{1}{4}|n\rangle + \tfrac{1}{2}|n{-}1\rangle,\qquad
  \mathbf{S}_n\cdot\mathbf{S}_{n+1}\,|n\rangle
  = -\tfrac{1}{4}|n\rangle + \tfrac{1}{2}|n{+}1\rangle.$$

Summing,

$$H|n\rangle
  \;=\; J\Bigl[\tfrac{N-2}{4} - \tfrac{1}{2}\Bigr]|n\rangle
        + \tfrac{J}{2}|n{-}1\rangle + \tfrac{J}{2}|n{+}1\rangle
  \;=\; \bigl(\tfrac{JN}{4} - J\bigr)|n\rangle
        + \tfrac{J}{2}\bigl(|n{-}1\rangle + |n{+}1\rangle\bigr).
\tag{1}$$

Diagonalise by a Fourier transform.  With
$|k\rangle = N^{-1/2}\sum_{n=1}^{N} e^{ikn}|n\rangle$ and PBC
$|n+N\rangle = |n\rangle$, the allowed momenta are $k = 2\pi m / N$
for $m = 0, 1, \dots, N-1$.  Acting with $H$,

$$H|k\rangle = \bigl(\tfrac{JN}{4} - J\bigr)|k\rangle
             + \tfrac{J}{2}(e^{ik} + e^{-ik})|k\rangle
             = \Bigl[\tfrac{JN}{4} - J(1 - \cos k)\Bigr]|k\rangle.$$

Thus the one-magnon eigenenergy is

$$E_{1}(k) \;=\; \tfrac{JN}{4} - J(1-\cos k)
          \;=\; \tfrac{JN}{4} - 2 J\sin^{2}(k/2).
\tag{2}$$

The quantity $\varepsilon(k) \equiv -2 J\sin^{2}(k/2)$ is the single-
magnon dispersion *relative to the ferromagnetic reference*.  It is
negative for every $k \ne 0$ — magnons lower the energy of the AF
chain.

### Step 2 — Two-magnon Bethe ansatz and phase shift

For $M = 2$ down-spins at positions $(n_1, n_2)$ with $n_1 < n_2$,
let $|n_1 n_2\rangle$ denote the corresponding basis state and write
the eigenstate ansatz

$$|\psi\rangle = \sum_{1 \le n_1 < n_2 \le N}
  \Psi(n_1, n_2)\,|n_1 n_2\rangle.$$

The same kinetic-energy bookkeeping as in Step 1 gives, when the two
down-spins are *not adjacent* ($n_2 > n_1 + 1$),

$$[H\Psi](n_1, n_2)
 = \bigl(\tfrac{JN}{4} - 2 J\bigr)\Psi(n_1, n_2)
 + \tfrac{J}{2}\bigl[\Psi(n_1{-}1, n_2) + \Psi(n_1{+}1, n_2)
                    + \Psi(n_1, n_2{-}1) + \Psi(n_1, n_2{+}1)\bigr].
\tag{3}$$

When the two down-spins are adjacent ($n_2 = n_1 + 1$), the bond
$(n_1, n_2)$ has two down-spins, so its diagonal contribution is
$+\tfrac{J}{4}$ (not $-\tfrac{J}{4}$) and its off-diagonal part
vanishes.  This modifies (3) into a *boundary condition* on
$\Psi(n_1, n_2)$ when $n_2 = n_1 + 1$.  Bethe's idea is to absorb the
boundary condition into a single "2-body phase shift" by writing

$$\Psi(n_1, n_2) = A\,e^{i(k_1 n_1 + k_2 n_2)}
                 + A'\,e^{i(k_2 n_1 + k_1 n_2)}.
\tag{4}$$

Substituting (4) into (3) diagonalises the bulk operator with
eigenvalue

$$E_2 = \tfrac{JN}{4} - J(1 - \cos k_1) - J(1 - \cos k_2)
     = \tfrac{JN}{4} + \varepsilon(k_1) + \varepsilon(k_2).
\tag{5}$$

The adjacency boundary condition determines the ratio $A'/A$.
Carrying out the algebra (see e.g. Karbach–Muller 1997 eqs. (7)–(10);
the calculation is standard and amounts to matching the two sides of
the scattering equation at $n_2 = n_1 + 1$) gives

$$\frac{A'}{A}
 \;=\; -\,\frac{e^{i(k_1 + k_2)} + 1 - 2\,e^{i k_1}}
              {e^{i(k_1 + k_2)} + 1 - 2\,e^{i k_2}}
 \;\equiv\; e^{i\theta(k_1, k_2)}.
\tag{6}$$

Taking the log of (6) after simplifying the trigonometry yields the
Bethe 2-body phase shift

$$\theta(k_1, k_2) \;=\; 2\arctan\!\bigl(\cot(k_1/2) - \cot(k_2/2)\bigr)
\tag{7}$$

(a unique branch choice is fixed by requiring $\theta \to 0$ when
$k_1 \to k_2$).  Changing variables to **rapidities**

$$\lambda_j \;\equiv\; \tfrac{1}{2}\cot(k_j/2),
\qquad \Longleftrightarrow\qquad
  e^{i k_j} = \frac{\lambda_j + i/2}{\lambda_j - i/2},
\tag{8}$$

equation (7) becomes

$$\theta(k_1, k_2) = 2\arctan(\lambda_1 - \lambda_2)
 = -i\ln\frac{\lambda_1 - \lambda_2 + i}{\lambda_1 - \lambda_2 - i}.
\tag{9}$$

### Step 3 — Bethe equations for $M$ magnons

For a state with $M$ magnons (down-spins at positions $n_1 < \dots <
n_M$), the ansatz generalises to a superposition over all $M!$
permutations of $\{k_1, \dots, k_M\}$:

$$\Psi(n_1,\dots,n_M)
 = \sum_{P \in S_M} A_P\,
   \exp\!\Bigl(i\sum_{j=1}^{M} k_{P(j)} n_j\Bigr),
\tag{10}$$

with the ratios $A_P/A_{P'}$ fixed by products of the 2-body phase
shifts (9) across transpositions — this is the *Yang–Baxter
factorisation* specific to integrable models.  Imposing periodic
boundary conditions $\Psi(n_1,\dots,n_M) = \Psi(n_2,\dots,n_M,
n_1 + N)$ on the ansatz gives, for each $j = 1,\dots,M$,

$$e^{i k_j N}
 = \prod_{\ell \ne j}(-e^{-i\theta(k_j, k_\ell)}).$$

Substituting (8) and (9) and simplifying,

$$\boxed{\;
  \left(\frac{\lambda_j + i/2}{\lambda_j - i/2}\right)^{N}
  = \prod_{\ell \ne j}^{M}\frac{\lambda_j - \lambda_\ell + i}
                                {\lambda_j - \lambda_\ell - i},
  \qquad j = 1,\dots, M.
\;}
\tag{11}$$

These are the **Bethe equations**.  Any solution $\{\lambda_j\}$
of (11) defines a Hamiltonian eigenstate via (10) with eigenvalue

$$E(\{\lambda_j\})
 \;=\; \tfrac{JN}{4} + \sum_{j=1}^{M}\varepsilon(k_j).
\tag{12}$$

Using the rapidity map (8), $\sin^2(k_j/2) = 1/(1 + \cot^2(k_j/2)) =
1/(1 + 4\lambda_j^2) = 1/[4(\lambda_j^2 + 1/4)]$, hence

$$\varepsilon(k_j) = -2 J\sin^{2}(k_j/2)
 = -\frac{J}{2(\lambda_j^{2} + 1/4)},$$

and (12) takes its rapidity form

$$\boxed{\;
E(\{\lambda_j\})
 \;=\; \tfrac{JN}{4}
       \;-\; \frac{J}{2}\sum_{j=1}^{M}\frac{1}{\lambda_j^{2} + 1/4}.
\;}
\tag{13}$$

The factor $\tfrac{1}{2}$ in front of the sum is essential — it is
the source of the factor of $2$ that distinguishes the Hulthén
answer $1/4 - \ln 2$ from the naive $1/4 - 2\ln 2$ one would get by
omitting it.

### Step 4 — Ground-state rapidity distribution

Taking the logarithm of (11) and using
$\tfrac{\lambda + i/2}{\lambda - i/2} = e^{-i\,2\arctan(2\lambda)}
\cdot e^{i\cdot(\text{integer multiple of } 2\pi)}$,

$$-2\arctan(2\lambda_j) \cdot N
 \;=\; -\sum_{\ell \ne j}^{M} 2\arctan(\lambda_j - \lambda_\ell)
       + 2\pi I_j,
\tag{14}$$

where $I_j$ are distinct half-integers (their parity fixed by
$M$).  Divide by $N$ and define the *counting function*

$$Y_N(\lambda) \;:=\; -\frac{1}{\pi}\arctan(2\lambda)
 + \frac{1}{\pi N}\sum_{\ell=1}^{M}\arctan(\lambda - \lambda_\ell).$$

Equation (14) reads $Y_N(\lambda_j) = I_j/N$.  In the ground state
the quantum numbers $I_j$ fill consecutively the smallest possible
integers symmetric about zero, so as $N \to \infty$ the
$\{\lambda_j\}$ become a dense, symmetric set on the real axis with
density $\rho(\lambda) > 0$ defined by

$$\rho(\lambda)\,d\lambda \;=\; \lim_{N\to\infty}
  \#\{j\;:\;\lambda_j \in [\lambda, \lambda + d\lambda]\}/N.$$

Differentiating $Y_N(\lambda_j) = I_j/N$ w.r.t. $j$ gives $Y_N'(\lambda)
\rho_N(\lambda) = 1/N$ on one side, and on the other — after replacing
the sum by an integral $(1/N)\sum_\ell \to \int \rho(\mu) d\mu$ in
the $N \to \infty$ limit — the derivative of $Y_N$ reads

$$Y'(\lambda)
 = -\frac{1}{\pi}\cdot\frac{2}{1 + 4\lambda^2}
   + \frac{1}{\pi}\int\frac{1}{1 + (\lambda - \mu)^2}\rho(\mu)\,d\mu.$$

Matching $Y'(\lambda) \cdot \rho \to \rho \cdot \rho$ at leading order
in $1/N$ (for the ground-state sector this simplifies to the linear
equation below; a careful derivation of the sign and factor of two is
in Takahashi 1999 §3.2) gives the **linear integral equation** for
$\rho$:

$$\boxed{\;
2\pi\,\rho(\lambda)
 \;=\; \frac{1}{\lambda^{2} + 1/4}
       \;-\; 2\int_{-\infty}^{\infty}\frac{\rho(\mu)}{(\lambda-\mu)^{2} + 1}\,d\mu.
\;}
\tag{15}$$

Equation (15) is the cornerstone of the Hulthén calculation.  The
normalisation constraint is the half-filling condition

$$\int_{-\infty}^{\infty}\rho(\lambda)\,d\lambda
 \;=\; \frac{M}{N}\Bigr|_{M = N/2}
 \;=\; \frac{1}{2}.
\tag{16}$$

### Step 5 — Solution by Fourier transform

Adopt the convention

$$\hat{f}(q) \;=\; \int_{-\infty}^{\infty} e^{-i q\lambda}\,f(\lambda)\,d\lambda,
\qquad
f(\lambda) \;=\; \frac{1}{2\pi}\int_{-\infty}^{\infty} e^{i q\lambda}\,\hat{f}(q)\,dq.$$

Two Fourier transforms we shall need (both are elementary contour
integrals; see e.g. Whittaker–Watson 1927 §6):

$$\int_{-\infty}^{\infty}\frac{e^{-i q\lambda}}{\lambda^{2}+a^{2}}\,d\lambda
 \;=\; \frac{\pi}{a}\,e^{-a|q|},
\qquad a > 0.
\tag{17}$$

With $a = \tfrac{1}{2}$, the Fourier transform of the inhomogeneous
term in (15) is $2\pi\,e^{-|q|/2}$.  With $a = 1$, the kernel
$K(\lambda) := 2/(\lambda^2 + 1)$ transforms to $\hat{K}(q) = 2\pi\,
e^{-|q|}$.  Fourier-transforming (15) and using the convolution
theorem on the $\rho \star K$ term,

$$2\pi\,\hat{\rho}(q)
 \;=\; 2\pi\,e^{-|q|/2}
       \;-\; 2\pi\,e^{-|q|}\,\hat{\rho}(q).$$

Solve for $\hat{\rho}$:

$$\hat{\rho}(q)\,\bigl(1 + e^{-|q|}\bigr)
 \;=\; e^{-|q|/2}
\quad\Longrightarrow\quad
\hat{\rho}(q) \;=\; \frac{e^{-|q|/2}}{1 + e^{-|q|}}
 \;=\; \frac{1}{e^{|q|/2} + e^{-|q|/2}}
 \;=\; \frac{1}{2\cosh(q/2)}.
\tag{18}$$

**Normalisation check.** $\hat{\rho}(0) = 1/(2\cosh 0) = 1/2$, and
$\int\rho\,d\lambda = \hat{\rho}(0) = 1/2$.  Equation (16) is
satisfied.

**Inverting** (18) uses the companion Fourier transform

$$\int_{-\infty}^{\infty}\frac{e^{-i q\lambda}}{\cosh(a\lambda)}\,d\lambda
 \;=\; \frac{\pi}{a\,\cosh(\pi q/(2 a))},
\qquad a > 0.
\tag{19}$$

This is the standard result obtained by closing the contour in the
lower half-plane and summing residues over $\lambda = -i(n+1/2)\pi/a$;
the sign pattern $(-1)^n$ telescopes to the $\cosh$ in the
denominator.  Applied with $a = \pi/2$,

$$\int_{-\infty}^{\infty}\frac{e^{-i q\lambda}}{\cosh(\pi\lambda/2)}\,d\lambda
 \;=\; \frac{2}{\cosh q}.$$

Rescaling $\lambda \to 2\lambda$ and $q \to q/2$ converts this to
the inverse Fourier transform of (18):

$$\rho(\lambda)
 \;=\; \frac{1}{2\pi}\int_{-\infty}^{\infty}
        \frac{e^{i q\lambda}}{2\cosh(q/2)}\,dq
 \;=\; \frac{1}{2\cosh(\pi\lambda)}.$$

Hence

$$\boxed{\;
\rho(\lambda) \;=\; \frac{1}{2\cosh(\pi\lambda)},
\qquad -\infty < \lambda < \infty.
\;}
\tag{20}$$

### Step 6 — Evaluation of the energy integral

From (13), the ground-state energy per site in the thermodynamic
limit is

$$e_0 \;=\; \frac{E_0(N)}{N}\bigg|_{N\to\infty}
 \;=\; \frac{J}{4} \;-\; \frac{J}{2}\cdot
        \lim_{N\to\infty}\frac{1}{N}\sum_{j=1}^{M}\frac{1}{\lambda_j^2 + 1/4}.$$

The continuum replacement $(1/N)\sum_j \to \int \rho(\lambda)\,
d\lambda$ converts the sum to

$$e_0 \;=\; \frac{J}{4} - \frac{J}{2}\int_{-\infty}^{\infty}
            \frac{\rho(\lambda)}{\lambda^{2} + 1/4}\,d\lambda.
\tag{21}$$

Write $I := \int \rho(\lambda)/(\lambda^2 + 1/4)\,d\lambda$ and
evaluate it by Parseval's theorem,

$$\int_{-\infty}^{\infty} f(\lambda)\,g(\lambda)\,d\lambda
 \;=\; \frac{1}{2\pi}\int_{-\infty}^{\infty}
        \hat{f}(q)\,\overline{\hat{g}(q)}\,dq,$$

with $f = \rho$ (so $\hat{f}(q) = 1/(2\cosh(q/2))$ by (18)) and
$g(\lambda) = 1/(\lambda^2 + 1/4)$ (so $\hat{g}(q) = 2\pi\,e^{-|q|/2}$
by (17)).  Both are real and even; complex conjugation leaves them
unchanged.  Therefore

$$I \;=\; \frac{1}{2\pi}\int_{-\infty}^{\infty}
         \frac{1}{2\cosh(q/2)}\cdot 2\pi\,e^{-|q|/2}\,dq
 \;=\; \int_{-\infty}^{\infty}\frac{e^{-|q|/2}}{2\cosh(q/2)}\,dq.$$

The integrand is even in $q$, so

$$I \;=\; 2\int_{0}^{\infty}\frac{e^{-q/2}}{2\cosh(q/2)}\,dq
 \;=\; \int_{0}^{\infty}\frac{e^{-q/2}}{\cosh(q/2)}\,dq.$$

Using $\cosh(q/2) = (e^{q/2} + e^{-q/2})/2$,

$$\frac{e^{-q/2}}{\cosh(q/2)}
 \;=\; \frac{2 e^{-q/2}}{e^{q/2} + e^{-q/2}}
 \;=\; \frac{2}{e^{q} + 1}.$$

So

$$I \;=\; 2\int_{0}^{\infty}\frac{dq}{e^{q} + 1}.$$

Substitute $u = e^{-q}$, $du = -e^{-q}\,dq$, $dq = -du/u$; limits
$q\in[0,\infty)$ map to $u\in(0,1]$.  In the new variable,
$1/(e^{q} + 1) = u/(1 + u)$ and $dq = -du/u$, so

$$\int_{0}^{\infty}\frac{dq}{e^{q} + 1}
 \;=\; \int_{1}^{0}\frac{u}{1 + u}\cdot\Bigl(-\frac{du}{u}\Bigr)
 \;=\; \int_{0}^{1}\frac{du}{1 + u}
 \;=\; \bigl[\ln(1 + u)\bigr]_{0}^{1}
 \;=\; \ln 2.$$

Therefore

$$\boxed{\;I \;=\; 2\ln 2.\;}
\tag{22}$$

### Step 7 — Assembly of $e_0$

Substitute (22) into (21):

$$e_0 \;=\; \frac{J}{4} - \frac{J}{2}\cdot 2\ln 2
 \;=\; \frac{J}{4} - J\ln 2
 \;=\; J\!\left(\frac{1}{4} - \ln 2\right).$$

This is the Main-result value.  Numerically $1/4 - \ln 2 \approx
0.25 - 0.6931 = -0.4431$, i.e. $e_0 \approx -0.4431\,J$.

### Step 8 — Limiting-case and finite-size checks

**(i) Dimer $N = 2$.** With $N = 2$ and PBC, the two bonds are
identical, $H = 2 J\,\mathbf{S}_1\cdot\mathbf{S}_2$.  The two-spin
Hilbert space decomposes into singlet ($S_{\rm tot} = 0$, energy
$2 J\cdot(-3/4) = -3 J/2$) and triplet ($S_{\rm tot} = 1$, energy
$+ J/2$); see the self-contained calculation at
[Heisenberg dimer: singlet–triplet splitting](heisenberg-dimer-singlet-triplet.md).
The per-bond ground-state energy is $-3J/4 = -0.75\,J$.  For $N = 2$
this overshoots the thermodynamic value $-0.4431\,J$: finite-size
corrections scale as $\sim -\pi^2/(6 N^2)$, and at $N = 2$ this
estimate is of order unity, consistent with the large overshoot.

**(ii) Finite-size PBC trend.** Numerical diagonalisation of the
$2^N$-dimensional chain for $N = 4, 6, 8, 10, \dots$ converges to
$e_0 = J(1/4 - \ln 2)$ as $1/N^2$ (for PBC, in the singlet sector).
The test suite verifies this in
`test/verification/test_universality_cross_check.jl`:
ED at $N = 4, 6, 8$ PBC cross-checks that
$|E_0(N)/N - e_0| < 5\%$ already at $N = 8$.

**(iii) Independent route — XXZ $\Delta\to 1^-$ limit.**  The same
Main-result value is the $\Delta = 1$ endpoint of the XXZ Luttinger-
parameter family derived in
[XXZ Luttinger parameters](xxz-luttinger-parameters.md): the
Luttinger velocity $u(\Delta) = J(\pi/2) \sin\gamma/\gamma$ with
$\gamma = \arccos\Delta$ tends to $\pi J/2$ as $\gamma\to 0$ (the
des Cloizeaux–Pearson spinon velocity), and the same Fourier machinery
that produces (18) in the isotropic case produces the XXZ dressed-
charge $Z^2 = \pi/(2(\pi - \gamma)) \to 1/2$ at $\Delta = 1$.  The
ground-state energy per site in the Heisenberg limit of that XXZ
calculation reproduces $J(1/4 - \ln 2)$ (Takahashi 1999 §4.3).

---

## References

- H. Bethe, *Zur Theorie der Metalle. I. Eigenwerte und Eigenfunktionen
  der linearen Atomkette*, Z. Physik **71**, 205 (1931).  Original
  Bethe ansatz.
- L. Hulthén, *Über das Austauschproblem eines Kristalles*,
  Ark. Mat. Astron. Fys. **26A**, No. 11 (1938).  Evaluation of the
  ground-state energy per site.
- C. N. Yang and C. P. Yang, *One-dimensional chain of anisotropic spin-spin
  interactions I*, Phys. Rev. **150**, 321 (1966).  Systematic treatment
  of the rapidity-density equation.
- M. Karbach and G. Müller, *Introduction to the Bethe ansatz I*,
  Comput. Phys. **11**, 36 (1997).  Pedagogical derivation of (11) and
  (13) including the 2-body phase shift.
- M. Takahashi, *Thermodynamics of One-Dimensional Solvable Models*
  (Cambridge University Press, 1999), §3.2 (derivation of (15))
  and §4.3 (Heisenberg limit as $\Delta = 1$ of XXZ).
- E. T. Whittaker and G. N. Watson, *A Course of Modern Analysis*,
  4th ed. (Cambridge University Press, 1927), §6 (contour integrals
  (17) and (19)).

## Used by

- [Heisenberg model page](../models/quantum/heisenberg.md) —
  thermodynamic-limit ground-state energy density.
- [XXZ model page](../models/quantum/xxz.md) — $\Delta = 1$ isotropic
  point.
- [XXZ Luttinger parameters note](xxz-luttinger-parameters.md) —
  Heisenberg limit as a cross-check of the same Fourier machinery.

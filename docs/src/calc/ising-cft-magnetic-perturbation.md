# Ising CFT + Magnetic Perturbation: Surviving Conserved Charges

## Main result

Consider the $c = 1/2$ Ising CFT (minimal model $\mathcal{M}(3, 4)$
— see [`ising-cft-primary-operators`](ising-cft-primary-operators.md))
perturbed by the spin primary operator $\sigma$ of conformal
weight $h_{\sigma} = 1/16$:

$$\mathcal{A} \;=\; \mathcal{A}_{\rm Ising\,CFT}
 \;+\; \lambda\,\int d^{2}x\,\sigma(x),
\qquad \lambda \in \mathbb{R}.$$

Zamolodchikov 1989 showed that of the infinitely many conserved
charges $Q_{s}$ of the unperturbed CFT (one for every odd
integer spin $s \ge 1$), **exactly eight** survive as exact
conservation laws of the perturbed theory:

$$\boxed{\;
\text{Surviving spins} \;=\; \{1,\,7,\,11,\,13,\,17,\,19,\,23,\,29\}.
\;}$$

These are the **exponents of the $E_{8}$ Lie algebra** (Coxeter
number $h = 30$, so the exponents run from $1$ up to $h - 1 = 29$
and come in complementary pairs summing to $h$).

Consequences:

1. The massive perturbed theory has a factorised elastic
   $S$-matrix (from eight conserved charges + Zamolodchikov's
   theorem).
2. The bootstrap equations with this particular constraint set have
   a unique minimal solution: the $E_{8}$ affine Toda field theory
   (see [`e8-mass-spectrum-derivation`](e8-mass-spectrum-derivation.md)).
3. The eight stable-particle mass ratios are the
   Perron–Frobenius eigenvector of the $E_{8}$ Cartan matrix, with
   $m_{2}/m_{1} = \varphi = (1 + \sqrt{5})/2$ the famous
   golden-ratio prediction.
4. The mass gap scales with the perturbation as
   $m_{1} \propto |\lambda|^{8/15}$ by dimensional analysis.

The contrasting perturbation by the energy density $\varepsilon$
(conformal weight $1/2$, thermal direction) preserves only
**infinitely many charges via a trivial mechanism**: the perturbed
theory is a free massive Majorana fermion, which is integrable but
has a single-particle spectrum with no $E_{8}$ structure.

---

## Setup

### Conventions

Euclidean 2D, complex coordinates $z = x^{1} + i x^{2}$,
$\bar z = x^{1} - i x^{2}$. Holomorphic and antiholomorphic
derivatives $\partial \equiv \partial_{z}$, $\bar\partial \equiv
\partial_{\bar z}$. For the Ising CFT we have $c = 1/2$ with
primary operators $\{\mathbb{I}, \sigma, \varepsilon\}$ of
conformal weights $(h, \bar h) = (0, 0), (1/16, 1/16), (1/2, 1/2)$.

### Stress tensor and conserved currents

The CFT stress tensor $T(z)$ (holomorphic) has conformal weight
$h = 2$ and is conserved: $\bar\partial T = 0$. Similarly
$\bar T(\bar z)$ is antiholomorphic-conserved. From $T$ we build a
tower of **conserved higher-spin currents**

$$T_{2} = T,\qquad T_{4} = \,:T^{2}:,\qquad T_{6} = \,:T^{3}:\,+\,
      \text{(descendant corrections)},\qquad \dots,$$

with $T_{2 k}$ of conformal weight $2 k$. The "descendant
corrections" are levels of the Virasoro descendants needed to make
each $T_{2 k}$ a primary of the $\mathcal{W}$-algebra (Zamolodchikov
1985); for our purposes we only need that such conserved
higher-spin tensors exist at every even weight.

Holomorphic conservation $\bar\partial T_{s} = 0$ implies, by Stokes,
the **conserved charge**

$$Q_{s} \;=\; \oint\frac{dz}{2\pi i}\,T_{s}(z),$$

one for every odd spin $s = 1, 3, 5, 7, \dots$ (the spin $s = s_{\rm
current} - 1$ — the "weight $2k$, spin $2k - 1$" convention for a
spin-$s$ current).  Infinitely many conservation laws.

### Goal

Determine which of the infinite sequence of $Q_{s}$ survives in the
*perturbed* theory $\mathcal{A} \to \mathcal{A} + \lambda\int
\sigma$, and derive the distinctive list of eight spins
$\{1, 7, 11, 13, 17, 19, 23, 29\}$.

---

## Calculation

### Step 1 — Modification of conservation laws under perturbation

Under $\lambda$-perturbation by an operator $\Phi$ of conformal
weight $(h_{\Phi}, \bar h_{\Phi})$, the Euclidean action gains a
term $\lambda\int d^{2}z\,\Phi$, and the holomorphic conservation
law $\bar\partial T_{s} = 0$ is modified to

$$\bar\partial T_{s}(z) \;=\; \lambda\cdot R_{s}[\Phi](z)
 \;+\; O(\lambda^{2}),
\tag{1}$$

where $R_{s}[\Phi]$ is determined by the OPE of $T_{s}(z)$ with the
perturbing operator integrated over its insertion point. Using
the perturbative formula (Zamolodchikov 1989 §2; Cardy 1996 §11)

$$R_{s}[\Phi](z)
 \;=\; \underset{w\to z}{\mathrm{Res}}\,\,T_{s}(z)\,\Phi(w) \quad
 \text{(singular OPE residue)}.
\tag{2}$$

The charge $Q_{s} = \oint (dz/2\pi i)\,T_{s}$ remains conserved
*iff* the right-hand side of (1) is a **total derivative** with
respect to $z$:

$$R_{s}[\Phi](z) \;=\; \partial\,\bigl(\text{something}\bigr)_{z}
\quad \Longrightarrow\quad \oint\frac{dz}{2\pi i}\,R_{s} = 0
\quad \Longrightarrow\quad Q_{s}\text{ is conserved at }O(\lambda).
\tag{3}$$

The **survival condition** (3) is the key criterion. It reduces to
a combinatorial statement about the OPE coefficients of $T_{s}$
with $\Phi$.

### Step 2 — $R_{s}[\Phi]$ as a descendant of $\Phi$

By general CFT operator-product-expansion (Di Francesco–Mathieu–
Sénéchal 1997 §6.6), the singular part of $T_{s}(z)\,\Phi(w)$ takes
the form

$$T_{s}(z)\,\Phi(w)
 \;=\; \sum_{k \ge 0}\frac{c_{s, k}\,L_{-k}^{\Phi}(w)}{(z - w)^{s - k}}
       \;+\; \text{regular},
\tag{4}$$

where $L_{-k}^{\Phi}$ are Virasoro descendants at level $k$ of
$\Phi$ (the original primary), and $c_{s, k}$ are OPE coefficients
determined by the conformal symmetry.

Taking the residue at $w \to z$ — i.e. extracting the coefficient
of $1/(z - w)$ — yields

$$R_{s}[\Phi](z) \;=\; c_{s,\,s - 1}\,L_{-(s - 1)}^{\Phi}(z).
\tag{5}$$

For (3) to hold, this single Virasoro descendant of $\Phi$ at level
$s - 1$ must be a **total $\partial$-derivative**, equivalently a
descendant of the form $L_{-1}^{\Phi}\,\mathcal{O}$ for some local
operator $\mathcal{O}$. (Recall $L_{-1}$ is the translation
generator, so $L_{-1} \mathcal{O} = \partial \mathcal{O}$.)

Therefore **$Q_{s}$ survives iff $L_{-(s - 1)}^{\Phi}$ is a
$\partial$-derivative**, i.e. iff $L_{-(s - 1)}^{\Phi}\in
L_{-1}\cdot(\text{level}\,(s - 2)\text{-descendants of }\Phi)$ at
level $s - 1$ of the Verma module $V_{\Phi}$.

### Step 3 — Null-vector constraint for $\Phi = \sigma$

For a generic primary $\Phi$ with no null vectors, $L_{-(s-1)}\Phi$
is linearly independent of $L_{-1}\cdot(\text{descendants})$ — the
survival condition fails at every spin except the trivial $s = 1$
(where $L_{0}\Phi = h_{\Phi}\,\Phi$ is the total-derivative
tautology). No nontrivial conserved charges survive.

The Ising $\sigma$ is special because it has a **null descendant**
at level 2 (from the Kac table: $\sigma = \phi_{2, 1}$ has a
level-$r s = 2 \cdot 1 = 2$ null vector). Explicitly (Di Francesco–
Mathieu–Sénéchal 1997 §8.1):

$$\bigl(L_{-2} \;-\; \tfrac{3}{4}\,L_{-1}^{2}\bigr)\sigma \;=\; 0
\quad\text{modulo total descendants.}
\tag{6}$$

This null relation implies that $L_{-2}\sigma$ is a
$\partial$-derivative: $L_{-2}\sigma = \tfrac{3}{4}L_{-1}^{2}\sigma
= \tfrac{3}{4}\partial^{2}\sigma$.

More generally, the null vector (6) propagates up the Verma module:
$L_{-1}^{n}$ acting on the level-2 null gives nullifications at
higher levels, which create "gaps" in the space of
non-derivative descendants. The survival condition (3) is
**non-trivially satisfied** at precisely those levels $s - 1$ where
the null-vector cascade produces a $\partial$-derivative-equivalent
class.

### Step 4 — Spin-counting via null cascades (result)

The survival condition (3) must be checked spin-by-spin. The
relevant object is the generating function of "surviving descendants"
— those in the Verma module $V_{\sigma}$ modulo $\partial$-exact
forms and the null-vector ideal generated by (6):

$$\chi_{\sigma}^{\rm surv}(q) \;=\; \sum_{n \ge 0}
 \dim\!\bigl[\,V_{\sigma}^{(n)} / \bigl(L_{-1}\,V_{\sigma}^{(n - 1)}
      \oplus (\text{null cascade at level }n)\bigr)\,\bigr]\,q^{n},$$

where $V_{\sigma}^{(n)}$ is the level-$n$ subspace of the Verma
module. The Kac null-vector propagation (applying $L_{-1}^{k}$ to
(6) for $k = 0, 1, 2, \ldots$) removes not only level 2 itself but
a whole family of higher levels, creating an intricate pattern of
"surviving" levels.

The full computation, carried out in Zamolodchikov 1989 §4 (see
also Delfino 2004 Proposition 1 and Di Francesco–Mathieu–Sénéchal
1997 §12.4 for streamlined re-derivations), yields the explicit
sequence of non-trivial surviving levels

$$s - 1 \;\in\; \{0,\,6,\,10,\,12,\,16,\,18,\,22,\,28\},$$

giving the conserved-charge spins

$$\boxed{\;
s \;\in\; \{1,\,7,\,11,\,13,\,17,\,19,\,23,\,29\}.
\;}$$

These are precisely the **exponents of the $E_{8}$ Lie algebra**
(Coxeter number $h = 30$; the exponents are the positive integers
less than $h$ coprime to $h$).  The full generating-function form,
corresponding to the $E_{8}$ "$W_{1+\infty}$ character" on the
perturbed theory, is Zamolodchikov 1989 eq. (4.5); we do not
reproduce its closed form here since the spin list above is all
that enters the bootstrap argument.

The Remark that ties this back to the $E_{8}$ structure:

> The spin spectrum $\{1, 7, 11, 13, 17, 19, 23, 29\}$ of the
> surviving conserved charges is precisely the multiset of
> exponents of the $E_{8}$ root system. By Dorey's theorem
> ([`e8-mass-spectrum-derivation`](e8-mass-spectrum-derivation.md)
> Step 4), any integrable field theory whose spin currents form
> exactly this multiset **is** the $E_{8}$ affine Toda theory, up
> to the overall mass scale.

So the counting identifies the perturbed Ising theory as the
$E_{8}$ massive integrable theory *without* constructing its
$S$-matrix explicitly; the $S$-matrix then follows from the
bootstrap program with $E_{8}$ structure imposed.

### Step 5 — Why $E_{8}$ exponents?

The exponents of $E_{8}$ are exactly $\{1, 7, 11, 13, 17, 19, 23,
29\}$, the numbers coprime to $30 = h_{E_{8}}$ in the range $[1,
h - 1]$. Coincidentally (though not coincidentally from the
representation-theory point of view), these are also the positive
integers $s < 30$ that:

- are not multiples of $2$, $3$, or $5$ — the primes dividing $30$,
- equivalently, are coprime to $h = 30$.

The appearance of these specific integers in Zamolodchikov's
null-vector counting for the $\sigma$-perturbation of $\mathcal{M}(3,
4)$ is a **structural coincidence** reflecting the hidden $E_{8}$
root lattice of the perturbed theory. Dorey 1991 turned this
coincidence into a theorem: any integrable perturbation by a
primary with the "right" null-vector structure gives an ADE-Toda
theory, and the specific perturbation $\sigma$ of $\mathcal{M}(3,
4)$ at level 2 is the one that produces $E_{8}$.

### Step 6 — Contrast: $\varepsilon$ perturbation → free massive Majorana

Perturbing by $\varepsilon$ instead gives:

$$\mathcal{A}_{\varepsilon} \;=\; \mathcal{A}_{\rm Ising\,CFT}
 \;+\; t\,\int d^{2}x\,\varepsilon(x),$$

with $t = \beta(T - T_{c})/J$ the thermal coupling (scaling
dimension $\Delta_{\varepsilon} = 1$). The $\varepsilon$ operator is
$\phi_{1, 2}$ in the Kac table with null vector at level
$r s = 1 \cdot 2 = 2$:

$$\bigl(L_{-2} \;-\; \tfrac{3}{4 h_{\varepsilon}}\,L_{-1}^{2}\bigr)
 \varepsilon \;=\; 0,
\qquad h_{\varepsilon} = \tfrac{1}{2},$$

i.e. $\bigl(L_{-2} - \tfrac{3}{2} L_{-1}^{2}\bigr)\varepsilon =
0$. The analogous counting for the $\varepsilon$-perturbation gives
surviving spins $s = 1, 3, 5, 7, \ldots$ — **all** odd spins
survive. This reflects the fact that the perturbed theory is a free
massive Majorana fermion, which has infinitely many conserved
currents (one for every odd spin, corresponding to products of the
fermion stress tensor) as any free field theory.

The difference between the two perturbations is geometric: the
$\varepsilon$-perturbation preserves the Bose-Fermi free structure
of the underlying Majorana, whereas the $\sigma$-perturbation
breaks the free-fermion symmetry (the $\sigma$ operator is
non-local in the fermion language — it is the **disorder
operator** dual to the free-fermion spin) and leaves behind only a
finite reminder of the infinite symmetry: the eight charges whose
spins match $E_{8}$ exponents.

### Step 7 — Mass gap scaling

The perturbation $\lambda\int\sigma$ has $[\lambda] =
[\text{mass}]^{2 - \Delta_{\sigma}} = [\text{mass}]^{15/8}$ by
dimensional analysis (since $\int d^{2}x$ has dimension $-2$ and
$\sigma(x)$ has dimension $\Delta_{\sigma} = 1/8$, and the action
must be dimensionless). The only mass scale in the perturbed
theory is $\lambda$, so

$$\boxed{\;
m_{1} \;\propto\; |\lambda|^{1/(2 - \Delta_{\sigma})}
 \;=\; |\lambda|^{1/(15/8)}
 \;=\; |\lambda|^{8/15}.
\;}
\tag{7}$$

The exact proportionality constant involves a beautiful integral
evaluation (Fateev 1994 eq. (1.1)):

$$m_{1} \;=\; C(\tfrac{8}{15})\,|\lambda|^{8/15},\qquad
 C(\tfrac{8}{15}) \;=\; \frac{2\,\Gamma(1/5)\,\Gamma(2/3)\,\Gamma(8/15)}
                              {\sqrt{\pi}\,\Gamma(3/5)\,\Gamma(5/6)\,\Gamma(8/15)}
 \;\times\;(\text{$E_{8}$-Toda overall normalisation}),$$

though the numerical value depends on the chosen normalisation of
$\sigma$ and $\lambda$.

### Step 8 — Experimental confirmation

The compound CoNb$_{2}$O$_{6}$ is a quasi-1D Ising ferromagnet
with residual inter-chain coupling acting as a small longitudinal
field. Near its quantum critical point ($B \approx B_{c} \approx
5.5\,\text{T}$ transverse field, $T \to 0$), Coldea et al. 2010
used inelastic neutron scattering to resolve two sharp
one-particle excitation modes. The observed mass ratio

$$\frac{m_{2}}{m_{1}} \;=\; 1.618 \pm 0.015$$

agrees with the $E_{8}$ prediction $\varphi = (1 + \sqrt{5})/2 =
1.6180\ldots$ within experimental error. A tentative identification
of the $m_{3}$ particle was reported (predicted $m_{3}/m_{1}
\approx 1.989$). This remains the most direct experimental
verification of a Lie algebra predicting particle masses in a
condensed-matter system.

---

## References

- A. B. Zamolodchikov, *Integrals of motion and S-matrix of the
  (scaled) $T = T_{c}$ Ising model with magnetic field*,
  Int. J. Mod. Phys. A **4**, 4235 (1989). The discovery paper;
  conserved-charge counting eq. (4.5).
- A. B. Zamolodchikov, *Infinite additional symmetries in
  two-dimensional conformal quantum field theory*, Theor. Math.
  Phys. **65**, 1205 (1985). Construction of higher-spin conserved
  currents in CFT (the $\mathcal{W}$-algebra).
- P. Dorey, *Root systems and purely elastic S-matrices*, Nucl.
  Phys. B **358**, 654 (1991). ADE/Toda integrable theories from
  null-vector perturbations.
- G. Delfino, *Integrable field theory and critical phenomena: the
  Ising model in a magnetic field*, J. Phys. A **37**, R45 (2004).
  Comprehensive review; Proposition 1 contains the surviving-spin
  counting with full null-vector accounting.
- V. A. Fateev, *The exact relations between the coupling constants
  and the masses of particles for the integrable perturbed
  conformal field theories*, Phys. Lett. B **324**, 45 (1994), eq.
  (1.1). Exact proportionality constant for $m_{1}(\lambda)$ in
  the $\sigma$-perturbed $\mathcal{M}(3, 4)$.
- J. L. Cardy, *Scaling and Renormalization in Statistical
  Physics*, Cambridge University Press (1996), §11.  OPE-based
  derivation of (1)–(3).
- P. Di Francesco, P. Mathieu, D. Sénéchal, *Conformal Field
  Theory*, Springer (1997), §6.6 and §8.1.  General OPE machinery
  (4)–(5), Ising null-vector equation (6).
- R. Coldea et al., *Quantum criticality in an Ising chain:
  experimental evidence for emergent E$_{8}$ symmetry*, Science
  **327**, 177 (2010). CoNb$_{2}$O$_{6}$ neutron-scattering
  experiment.

## Used by

- [`e8-mass-spectrum-derivation`](e8-mass-spectrum-derivation.md) —
  the $E_{8}$ mass formula takes Zamolodchikov's surviving-charge
  list as the input that identifies the perturbed theory as the
  $E_{8}$ Toda theory.
- [`ising-cft-primary-operators`](ising-cft-primary-operators.md) —
  the $\sigma$ primary and its null vector are derived / enumerated
  there.
- [$E_{8}$ universality class](../universalities/e8.md) — the
  mass-ratio universality class realised by this perturbation.
- [TFIM model page](../models/quantum/tfim.md) — lattice origin of
  the perturbation (small residual longitudinal field at $h = J$).

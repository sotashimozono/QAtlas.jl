# Kramers–Wannier Duality

## Main result

For the isotropic 2D classical Ising model on a square lattice with
reduced coupling $K = \beta J$, the high- and low-temperature
expansions of the partition function map into each other under the
**Kramers–Wannier duality**

$$\boxed{\;
\tanh K^{*} \;=\; e^{-2K}, \qquad
\sinh(2K)\,\sinh(2K^{*}) \;=\; 1.
\;}$$

The self-dual fixed point $K = K^{*}$ gives the exact critical
temperature of the 2D Ising model,

$$\boxed{\;
K_{c} \;=\; \tfrac{1}{2}\ln\!\bigl(1 + \sqrt{2}\bigr),
\qquad
T_{c} \;=\; \frac{2 J}{\ln(1 + \sqrt{2})}
\;\approx\; 2.269\,J.
\;}$$

The 1D quantum TFIM $H = -J\sum \sigma^z_i \sigma^z_{i+1} - h\sum
\sigma^x_i$ inherits an **operator-level duality** from the time-
continuum limit of the 2D Ising transfer matrix: the ordered
($h < J$) and disordered ($h > J$) phases are interchanged by a
unitary $U$ with

$$\boxed{\;
U\,H(J, h)\,U^{\dagger} \;=\; H(h, J),\qquad
U\,\sigma^z_i \sigma^z_{i+1}\,U^{\dagger} \;=\; \sigma^x_{i+1},
\qquad
U\,\sigma^x_i\,U^{\dagger} \;=\; \sigma^z_{i-1}\sigma^z_{i}.
\;}$$

The quantum critical point $h = J$ is the self-dual point of this
operator duality.

---

## Setup

### 2D classical Ising model

On an $L_x \times L_y$ square lattice with periodic boundary conditions,
classical spins $s_v \in \{\pm 1\}$ live on vertices $v$ and the
Hamiltonian is

$$H_{\rm cl} = -J\sum_{\langle u, v\rangle} s_u s_v,
\qquad J > 0,$$

with the sum over nearest-neighbour bonds. The partition function
is

$$Z = \sum_{\{s\}} e^{-\beta H_{\rm cl}}
    = \sum_{\{s\}} \prod_{\langle u, v\rangle} e^{K s_u s_v},
\qquad K \equiv \beta J.$$

### 1D quantum TFIM

The 1D transverse-field Ising model

$$H = -J\sum_{i=1}^{N}\sigma^z_i\sigma^z_{i+1}
      - h\sum_{i=1}^{N}\sigma^x_i$$

(PBC) is obtained from the 2D Ising model by taking the **continuous-
time limit** of the row-to-row transfer matrix: the 2D partition
function becomes $Z = \mathrm{Tr}\,e^{-\beta_{\rm eff} H}$
with $\beta_{\rm eff} \propto L_y$ and an anisotropic identification
between $J, h$ and the two 2D couplings (Mattis §§3.4–3.5; Sachdev
§5.2).  The self-dual structure of 2D Ising therefore descends to an
operator symmetry of the 1D TFIM.

### Goal

Derive $\tanh K^{*} = e^{-2K}$ from first principles by computing
both the high-$T$ and low-$T$ expansions of $Z$, identifying their
common form, and reading off the duality relation.  Then extract the
critical temperature and translate the duality to the 1D TFIM
operator algebra.

---

## Calculation

### Step 1 — High-temperature expansion

For a single bond $(u, v)$ with $s_u, s_v \in \{\pm 1\}$, the product
$s_u s_v \in \{\pm 1\}$ too. Expanding the exponential,

$$e^{K s_u s_v}
 = \cosh K + s_u s_v \sinh K
 = \cosh K\,\bigl(1 + s_u s_v\,\tanh K\bigr),$$

using $\cosh K\pm\sinh K = e^{\pm K}$ and splitting by parity in
$s_u s_v$. The full Boltzmann weight factorises:

$$e^{-\beta H_{\rm cl}}
 = (\cosh K)^{N_b}\,
   \prod_{\langle u, v\rangle}\bigl(1 + s_u s_v\,\tanh K\bigr),$$

where $N_b = L_x L_y \cdot 2$ is the total bond count (two bonds per
vertex on the square lattice, PBC). Summing over spins,

$$Z = (\cosh K)^{N_b}
      \sum_{\{s\}}\prod_{\langle u, v\rangle}
      \bigl(1 + s_u s_v \tanh K\bigr).$$

Expand the product into $2^{N_b}$ terms, one for each subset
$\Gamma \subseteq E$ of the bond set $E$. Each term looks like

$$(\tanh K)^{|\Gamma|}\cdot\prod_{(u,v)\in\Gamma} s_u s_v,$$

so

$$Z = (\cosh K)^{N_b}
      \sum_{\Gamma \subseteq E}(\tanh K)^{|\Gamma|}
      \sum_{\{s\}}\prod_{(u,v)\in\Gamma} s_u s_v.$$

The spin sum factorises over vertices. Each vertex $v$ appears in
the product as $s_v^{n_v(\Gamma)}$ where $n_v(\Gamma)$ is the number
of edges of $\Gamma$ incident to $v$. Because $s_v \in \{\pm 1\}$,

$$\sum_{s_v = \pm 1} s_v^{n_v} =
\begin{cases} 2 & n_v\text{ even},\\ 0 & n_v\text{ odd}.\end{cases}$$

So $\sum_{\{s\}}\prod_{(u,v)\in\Gamma} s_u s_v$ vanishes unless
**every** vertex has even $\Gamma$-degree. Such subsets $\Gamma$ are
called *closed graphs* (each vertex is an even-degree endpoint,
equivalent to a disjoint union of closed cycles on the lattice). If
$\Gamma$ is closed, the sum equals $2^{N_v}$ where $N_v = L_x L_y$
is the number of vertices.  Therefore

$$\boxed{\;
Z = (\cosh K)^{N_b}\,2^{N_v}\,
     \sum_{\Gamma \in \mathcal{C}} (\tanh K)^{|\Gamma|},
\;}
\tag{HT}$$

where $\mathcal{C}$ is the set of all closed subgraphs of the square
lattice. Equation (HT) is the *high-temperature expansion*: the
factor $\tanh K$ is small for $K \ll 1$, so the leading contributions
come from $\Gamma = \emptyset$ (weight 1), then short closed loops
(weight $\tanh^4 K$ for an elementary plaquette, etc.).

### Step 2 — Low-temperature expansion

Now expand $Z$ around the ordered ground state $s_v = +1$ for all
$v$, of energy $E_0 = -J N_b$. A configuration $\{s\}$ is specified
by the set of *domain walls*: bonds $(u, v)$ across which $s_u
s_v = -1$. Place a dual-lattice edge perpendicular to each domain
wall; this dual edge lives on the **dual square lattice** $\Lambda^{*}$
(another square lattice, shifted by $(1/2, 1/2)$).

A closed domain-wall configuration is a subset $\Gamma^{*} \subseteq
E^{*}$ of dual edges such that every dual vertex has even
$\Gamma^{*}$-degree (domain walls are *oriented* only through the
sign flip, but the boundary of an Ising-spin region is a closed
cycle of dual edges). There is a 2-to-1 correspondence between
$\{s\}$ configurations and $\Gamma^{*}$ configurations (flipping all
spins gives the same $\Gamma^{*}$), so

$$\sum_{\{s\}}\dotsb = 2\sum_{\Gamma^{*} \in \mathcal{C}^{*}}\dotsb,$$

with $\mathcal{C}^{*}$ the set of closed subgraphs of the dual
lattice.

Each domain wall costs an energy $2 J$ relative to the ordered
reference (the bond $(u, v)$ had $s_u s_v = -1$ instead of $+1$), so
for a configuration with $|\Gamma^{*}|$ broken bonds,

$$E(\{s\}) = E_0 + 2 J\,|\Gamma^{*}|,\qquad
e^{-\beta E(\{s\})} = e^{-\beta E_0}\,e^{-2 K\,|\Gamma^{*}|}.$$

Hence

$$\boxed{\;
Z = 2\,e^{-\beta E_0}
     \sum_{\Gamma^{*} \in \mathcal{C}^{*}}
     (e^{-2K})^{|\Gamma^{*}|}
 = 2\,e^{K N_b}
     \sum_{\Gamma^{*} \in \mathcal{C}^{*}}
     (e^{-2K})^{|\Gamma^{*}|}.
\;}
\tag{LT}$$

Equation (LT) is the *low-temperature expansion*: the factor
$e^{-2K}$ is small for $K \gg 1$, so the leading contributions come
from $\Gamma^{*} = \emptyset$ (all spins aligned), then single
flipped spins (a rectangle of 4 dual edges), etc.

### Step 3 — Matching (HT) and (LT)

The square lattice is **self-dual**: the dual of the square lattice
is another square lattice (up to translation). So $\mathcal{C} =
\mathcal{C}^{*}$ as abstract graph-combinatorial objects — both are
"closed subgraphs of a square lattice with the same number of
vertices and edges".

The HT expansion has weight $\tanh K$ per edge and prefactor
$(\cosh K)^{N_b} 2^{N_v}$. The LT expansion has weight $e^{-2K^{*}}$
per edge (reading $K^{*}$ instead of $K$ in (LT), which corresponds
to the *dual* coupling) and prefactor $2 e^{K^{*} N_b}$. Equating
the two sums over $\Gamma$ term-by-term forces

$$\boxed{\;\tanh K^{*} \;=\; e^{-2K}.\;}
\tag{1}$$

This is the Kramers–Wannier duality. Its geometric content is that
the high-$T$ closed-loop expansion of the original lattice becomes
the low-$T$ domain-wall expansion of the dual lattice, and on a
self-dual lattice the two are structurally identical with coupling
renamed.

**Equivalent symmetric form.** Taking $\sinh(2K) = 2\sinh K \cosh
K$ and $\cosh(2K) = \cosh^2 K + \sinh^2 K$, one finds

$$\sinh(2K)\,\sinh(2K^{*})
 = 2\,\frac{\sinh K}{\cosh K}\cdot\cosh^2 K
   \cdot 2\,\frac{\sinh K^{*}}{\cosh K^{*}}\cdot\cosh^2 K^{*}
 = 4\,\tanh K\,\tanh K^{*}\,\cosh^2 K\,\cosh^2 K^{*}.$$

Using (1) to substitute $\tanh K^{*} = e^{-2K}$ and $\tanh K =
e^{-2K^{*}}$ (by symmetry of (1)):

$$4\tanh K\,\tanh K^{*}\,\cosh^2 K\,\cosh^2 K^{*}
 = 4\,e^{-2K^{*}}\,e^{-2K}\,\cosh^2 K\,\cosh^2 K^{*}.$$

Now $\cosh^2 K = \tfrac{1}{4}(e^K + e^{-K})^2 = \tfrac{1}{4}
(e^{2K} + 2 + e^{-2K})$ and likewise for $K^{*}$, so
$4\,e^{-2K^{*}}\cosh^2 K^{*} = e^{-2K^{*}}(e^{2K^{*}} + 2 + e^{-2K^{*}})
= 1 + 2 e^{-2K^{*}} + e^{-4K^{*}}$. Using $e^{-2K^{*}} = \tanh K$
this simplifies to $(1 + \tanh K)^2 = 4\cosh^2 K/(e^K + e^{-K})^2
\cdot (e^K)^2 / 1$ — which gets cumbersome. A cleaner route is

$$\sinh(2K^{*})
 = \frac{2\tanh K^{*}}{1 - \tanh^2 K^{*}}
 = \frac{2\,e^{-2K}}{1 - e^{-4K}}
 = \frac{2\,e^{-2K}\cdot e^{2K}}{e^{2K} - e^{-2K}}
 = \frac{2}{2\sinh(2K)}
 = \frac{1}{\sinh(2K)}.$$

Hence

$$\sinh(2K)\,\sinh(2K^{*}) = 1,
\tag{2}$$

the symmetric form of the duality.

### Step 4 — Self-dual fixed point

Imposing $K = K^{*}$ in (2) gives $\sinh(2K_c)^2 = 1$, i.e.
$\sinh(2K_c) = 1$ (positive root for physical $K > 0$). Inverting,

$$2 K_c = \ln\!\bigl(1 + \sqrt{2}\bigr)
\quad\Longleftrightarrow\quad
K_c = \tfrac{1}{2}\ln\!\bigl(1 + \sqrt{2}\bigr),$$

using $\sinh x = 1 \Leftrightarrow e^x - e^{-x} = 2 \Leftrightarrow
e^{2x} - 2 e^x - 1 = 0 \Leftrightarrow e^x = 1 + \sqrt{2}$ (positive
root).  Converting to temperature $K_c = J/(k_B T_c)$ with $k_B = 1$,

$$T_c = \frac{J}{K_c} = \frac{2 J}{\ln(1 + \sqrt{2})}
      \;\approx\; 2.269\,J.$$

Onsager 1944 confirmed independently that this self-dual point is
indeed the phase-transition temperature (not just a self-dual
fixed point — the duality argument alone shows that *if* there is a
unique critical point, it must sit at $K_c$, but not that $K_c$ is
critical). The independent proof is in the
[`yang-magnetization-toeplitz`](yang-magnetization-toeplitz.md)
derivation, which exhibits the non-analyticity explicitly.

### Step 5 — Operator duality for the 1D TFIM

The 1D TFIM inherits the 2D duality via the transfer-matrix
identification. We derive the operator form directly by constructing
dual Pauli operators on the bonds of the original chain.

**Dual operators.** Labelling the bonds by $b = 1, \dots, N$ (bond $b$
connects sites $b$ and $b+1$ with PBC identification $N+1 \equiv 1$),
define

$$\tau^x_b \;\equiv\; \sigma^z_b\,\sigma^z_{b+1},
\qquad
\tau^z_b \;\equiv\; \prod_{k \le b}\sigma^x_k.
\tag{3}$$

**Pauli algebra on each bond.** Check $(\tau^x_b)^2 = (\sigma^z_b)^2
(\sigma^z_{b+1})^2 = 1$ and $(\tau^z_b)^2 = \prod_k (\sigma^x_k)^2 =
1$. Anticommutation on the same bond:

$$\tau^x_b\,\tau^z_b
 = \sigma^z_b\sigma^z_{b+1}\prod_{k\le b}\sigma^x_k
 = -\sigma^z_b\sigma^x_b\sigma^z_{b+1}\prod_{k < b}\sigma^x_k
 = -\prod_{k\le b}\sigma^x_k \,\sigma^z_b\,\sigma^z_{b+1}
 = -\tau^z_b\,\tau^x_b,$$

using $\sigma^z\sigma^x = -\sigma^x\sigma^z$ on site $b$ once and
commuting the remaining factors (all acting on different sites). The
same computation as in the
[`jw-tfim-bdg`](jw-tfim-bdg.md) derivation shows that operators on
different bonds commute. Therefore $\{\tau^x_b, \tau^z_b\}$ satisfy
the spin-$\tfrac{1}{2}$ algebra.

**Rewriting $H$.** The Ising bond term is immediate from (3):

$$-J\sum_{i=1}^{N}\sigma^z_i\sigma^z_{i+1}
 = -J\sum_{b=1}^{N}\tau^x_b.$$

For the transverse field, use $\tau^z_{b-1}\tau^z_b =
\prod_{k \le b-1}\sigma^x_k \cdot \prod_{k \le b}\sigma^x_k =
\sigma^x_b$ (the overlap on sites $1, \dots, b-1$ squares to 1,
leaving just the $b$-th factor of the longer product):

$$-h\sum_{i=1}^{N}\sigma^x_i
 = -h\sum_{b=1}^{N}\tau^z_{b-1}\,\tau^z_b,$$

with $\tau^z_0 \equiv 1$ (empty product). Hence

$$H \;=\; -J\sum_{b=1}^{N}\tau^x_b
       \;-\; h\sum_{b=1}^{N}\tau^z_{b-1}\,\tau^z_b.
\tag{4}$$

Compare to the original Hamiltonian in $\sigma$ operators:

$$H(J, h)\bigl|_{\sigma} \;=\; -J\sum\sigma^z\sigma^z - h\sum\sigma^x.$$

Relabelling bonds as sites ($b \to i$), the transformation
$\sigma \to \tau$ interchanges the Ising bond term and the
transverse-field term with the couplings swapped:

$$\boxed{\;
H(J, h)\bigl|_{\sigma} \;=\; H(h, J)\bigl|_{\tau}.
\;}
\tag{5}$$

So the duality is a **unitary operator $U$** that implements $\sigma
\to \tau$ in (3): $U^{\dagger}\sigma U = \tau$. Under $U$,

$$U\,\sigma^z_i\sigma^z_{i+1}\,U^{\dagger}
 = \tau^x_i = \sigma^x_{i+1}\bigl|_{\rm after\ relabelling},
\qquad
U\,\sigma^x_i\,U^{\dagger}
 = \tau^z_{i-1}\tau^z_i
 = \sigma^z_{i-1}\sigma^z_i\bigl|_{\rm after\ relabelling}.$$

### Step 6 — Self-dual TFIM critical point

The self-dual point of the TFIM operator duality (5) is $J = h$.
The 2D-duality ancestor of this point is the 2D Ising $K = K^{*}$
fixed point, $\sinh(2K_c) = 1$. The descent goes via the anisotropic
2D Ising transfer matrix, which has two couplings $(K_x, K_y)$ and a
self-dual condition $\sinh(2K_x)\sinh(2K_y) = 1$ — identifying $K_x
\leftrightarrow J$ and $\tau_y^{-1}\ln\tanh K_y \leftrightarrow h$
in the continuous-time limit $\tau_y \to 0$ with $K_y, \tau_y$
tuned such that $\tanh K_y \sim \tau_y\cdot h$, the self-dual
condition reduces to $J = h$ (Sachdev 2011 §5.2).

### Step 7 — Limiting-case checks

**Duality maps the phases.** From (1), $K \to 0$ sends $K^{*} \to
\infty$ and vice versa. So high-$T$ (disordered, $\langle s\rangle =
0$) on the original lattice corresponds to low-$T$ (ordered,
$\langle s\rangle \ne 0$) on the dual, and the two lattices share
critical properties at the self-dual point.

**Quantum check: $h = 0$.** At $h = 0$ the TFIM is the classical
Ising chain with ground state $|{\uparrow\cdots\uparrow}\rangle$ or
$|{\downarrow\cdots\downarrow}\rangle$. Under (5) this maps to the
$J = 0$ paramagnetic chain with ground state $|{\leftarrow\cdots
\leftarrow}\rangle$, the trivial paramagnet. The duality interchanges
ordered / disordered as expected.

**Quantum check: $h = J$.** The critical TFIM; the duality fixes
this point. The spectrum is gapless and conformally invariant with
central charge $c = 1/2$ — the same Ising CFT that describes the
2D classical Ising point at $K_c$.

**Graphene / square-lattice 3D check.** The square lattice is
self-dual; the honeycomb and triangular lattices are *dual* to each
other (not self-dual). The Kramers–Wannier relation for those pairs
has $\sinh(2 K_{\rm hc})\sinh(2 K_{\rm tri}) = 1$ with different
critical couplings on the two lattices, recovering the famous Wannier
1945 result for the triangular-lattice Ising critical coupling.

---

## References

- H. A. Kramers and G. H. Wannier, *Statistics of the two-
  dimensional ferromagnet. Part I*, Phys. Rev. **60**, 252 (1941).
  Original high-temperature / low-temperature duality argument.
- L. Onsager, *Crystal statistics. I. A two-dimensional model with
  an order–disorder transition*, Phys. Rev. **65**, 117 (1944).
  Independent exact solution confirming that the self-dual point is
  the phase-transition temperature.
- F. Wegner, *Duality in generalized Ising models and phase
  transitions without local order parameters*, J. Math. Phys. **12**,
  2259 (1971). Generalisation to gauge Ising models.
- D. C. Mattis, *Statistical Mechanics Made Simple*, 2nd ed. (World
  Scientific, 2008), §§3.4–3.5. Transfer-matrix derivation of the
  1D-quantum / 2D-classical correspondence.
- S. Sachdev, *Quantum Phase Transitions*, 2nd ed. (Cambridge
  University Press, 2011), §5.2. Continuous-time-limit derivation
  of the 1D TFIM from the 2D Ising transfer matrix and the
  TFIM self-duality $J \leftrightarrow h$.
- J. B. Kogut, *An introduction to lattice gauge theory and spin
  systems*, Rev. Mod. Phys. **51**, 659 (1979), §III–§IV. Pedagogical
  exposition of 1D TFIM self-duality from the 2D Ising perspective.

## Used by

- [`jw-tfim-bdg.md`](jw-tfim-bdg.md) — the KW step is the operator
  prerequisite for Jordan–Wigner on the TFIM.
- [IsingSquare model page](../models/classical/ising-square.md) —
  self-dual critical temperature $T_c = 2J/\ln(1 + \sqrt{2})$.
- [TFIM model page](../models/quantum/tfim.md) — ordered /
  disordered phase interchange at $h = J$.

# Bloch Hamiltonian Builder

## What It Is

The generic Bloch Hamiltonian builder `bloch_tb_spectrum` constructs
the momentum-space tight-binding Hamiltonian $H(\mathbf{k})$ for
**any** 2D lattice topology defined by a `Lattice2D` object. It
reads the unit cell structure (sublattice positions, hopping
vectors) from the lattice and builds an $n_{\mathrm{sub}} \times
n_{\mathrm{sub}}$ Bloch matrix at each allowed wavevector
$\mathbf{k}$, where $n_{\mathrm{sub}}$ is the number of sublattices.

---

## How It Works

### 1. Unit Cell from Lattice2D

The builder extracts the unit cell information from a `Lattice2D`
topology:

- **Sublattice count** $n_{\mathrm{sub}}$: number of sites per unit
  cell (e.g., 2 for honeycomb/kagome structures, 3 for kagome, etc.)
- **Hopping vectors** $\{\boldsymbol{\delta}_\alpha\}$: real-space
  displacements connecting sublattice sites, with sublattice indices

### 2. Bloch Matrix Construction

For each wavevector $\mathbf{k} = (k_x, k_y)$ on the grid of
allowed momenta ($L_x \times L_y$ points in the Brillouin zone), the
$n_{\mathrm{sub}} \times n_{\mathrm{sub}}$ Bloch Hamiltonian is

$$H_{ab}(\mathbf{k}) = t \sum_{\boldsymbol{\delta}: a \to b} e^{i\mathbf{k}\cdot\boldsymbol{\delta}}$$

where the sum runs over all hopping vectors connecting sublattice
$a$ to sublattice $b$.

### 3. Diagonalization

The eigenvalues of $H(\mathbf{k})$ at each $\mathbf{k}$ give the
band energies $E_n(\mathbf{k})$. The full spectrum is the collection
of all eigenvalues over all allowed $\mathbf{k}$ points, sorted.

---

## Key Design Principle: Topology-Agnostic

The builder does **not** contain any lattice-specific formulas. It
works identically for honeycomb (graphene), kagome, Lieb, triangular,
and any future `Lattice2D` topology. This generality is what makes
it valuable as a verification tool: it provides an independent
computation path that is completely generic.

---

## Verification: Three Independent Paths

For each tight-binding model in QAtlas, three independent
computations of the spectrum exist:

| Path | Method | Lattice-specific? |
|------|--------|-------------------|
| **A**: Hardcoded Bloch formula | Closed-form $E_n(\mathbf{k})$ in `src/` | Yes |
| **B**: Generic Bloch builder | `bloch_tb_spectrum(Topology, Lx, Ly, t)` | No |
| **C**: Real-space ED | `build_tight_binding(bonds, N; t)` + `eigvals` | No |

All three paths agree to machine precision for every topology.
This three-way agreement provides strong evidence that both the
analytical formulas and the generic builder are correct.

### Verified Topologies

| Topology | $n_{\mathrm{sub}}$ | Key spectral feature | Agreement |
|----------|---------------------|----------------------|-----------|
| Honeycomb (graphene) | 2 | Dirac cones at $K, K'$ | $< 10^{-12}$ |
| Kagome | 3 | Flat band at $E = +2t$ | $< 10^{-12}$ |
| Lieb | 3 | Flat band at $E = 0$ | $< 10^{-12}$ |
| Triangular | 1 | Van Hove singularities, range $[-6t, +3t]$ | $< 10^{-12}$ |

---

## API

```julia
using QAtlas

# Generic builder: works for ANY Lattice2D topology
spectrum = bloch_tb_spectrum(Honeycomb, 6, 6, 1.0)

# Compare against hardcoded formula
spectrum_exact = QAtlas.fetch(Graphene(), BlochSpectrum(); Lx=6, Ly=6, t=1.0)

# They agree to machine precision
@assert spectrum ≈ spectrum_exact
```

---

## References

- N. W. Ashcroft, N. D. Mermin, *Solid State Physics*
  (Holt, Rinehart and Winston, 1976), Ch. 8 --- Bloch's theorem.
- A. H. Castro Neto, F. Guinea, N. M. R. Peres, K. S. Novoselov,
  A. K. Geim, "The electronic properties of graphene",
  Rev. Mod. Phys. **81**, 109 (2009) --- honeycomb tight-binding.

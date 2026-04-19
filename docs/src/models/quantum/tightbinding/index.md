# Tight-Binding Models

## Overview

The nearest-neighbor tight-binding (TB) Hamiltonian on a periodic
lattice is one of the simplest quantum models, yet it already captures
essential physics: Dirac cones (graphene), flat bands (kagome, Lieb),
frustration-induced band asymmetry (triangular), and Van Hove
singularities.

$$H = -t \sum_{\langle i,j \rangle} \bigl(c^\dagger_i c_j + c^\dagger_j c_i\bigr)$$

where $t > 0$ is the nearest-neighbor hopping amplitude and the sum
runs over all nearest-neighbor pairs $\langle i,j \rangle$ determined
by the lattice topology.

**Parameters**: Hopping amplitude $t$ (default 1.0).

**Boundary conditions**: All TB spectra in QAtlas use periodic boundary
conditions (PBC) in both directions on an $L_x \times L_y$ supercell
of unit cells.

---

## Two Computation Paths

QAtlas implements two independent routes to the TB spectrum, enabling
cross-validation:

1. **Hardcoded Bloch formulas** (`src/models/quantum/tightbinding/`):
   For each lattice with a dedicated source file, the closed-form
   dispersion relation is evaluated analytically at each allowed
   $\mathbf{k}$-point. This is maximally efficient and transparent.

2. **Generic Bloch builder** (`test/util/bloch.jl`):
   Given *any* `Lattice2D` topology, the function `bloch_tb_spectrum`
   reads the unit cell definition (sublattice positions and
   connections) via `get_unit_cell`, constructs the $n_{\text{sub}}
   \times n_{\text{sub}}$ Bloch Hamiltonian

   $$H_{\alpha\beta}(\mathbf{k}) = -t \sum_{\text{connections}\; \beta \to \alpha} e^{i(d_x \theta_1 + d_y \theta_2)} + \text{h.c.}$$

   at each of the $L_x \cdot L_y$ allowed momenta, and diagonalizes it
   numerically. This requires **no model-specific formulas** and
   serves as a fully generic cross-check.

When both paths agree (verified in `test_bloch_generic.jl`), the
hardcoded formulas and the generic builder are mutually validated.

For lattices without dedicated source files (Dice, Shastry--Sutherland,
Union Jack, Square), only the generic builder and real-space ED are
used; both agree in every tested case.

---

## Lattice Catalogue

| Lattice | $n_{\text{sub}}$ | Bipartite? | Flat band? | `src/` formula? | Page |
|---------|:-----------------:|:----------:|:----------:|:---------------:|:----:|
| Square | 1 | Yes | No | No | -- |
| Honeycomb (Graphene) | 2 | Yes | No | Yes | [honeycomb.md](honeycomb.md) |
| Kagome | 3 | No | Yes ($+2t$) | Yes | [kagome.md](kagome.md) |
| Lieb | 3 | Yes | Yes ($0$) | Yes | [lieb.md](lieb.md) |
| Triangular | 1 | No | No | Yes | [triangular.md](triangular.md) |
| Dice ($T_3$) | 3 | Yes | Yes ($0$) | No | -- |
| Shastry--Sutherland | 4 | No | No | No | -- |
| Union Jack | 2 | No | No | No | -- |

**Bipartite**: The lattice graph admits a two-colouring (A/B sublattice
decomposition). Bipartite lattices have chiral (sublattice) symmetry,
so the spectrum is symmetric about $E = 0$.

**Flat band**: A completely dispersionless band (zero bandwidth).
Flat bands arise from destructive interference in the hopping network,
typically associated with a sublattice imbalance in bipartite lattices
(Lieb, Dice) or geometric frustration (Kagome).

---

## Generic Bloch Builder

The generic builder in `test/util/bloch.jl` works for any 2D periodic
topology that `Lattice2D` supports:

```julia
include("test/util/bloch.jl")

# Spectrum of the 3×3 honeycomb TB model
λ = bloch_tb_spectrum(Honeycomb, 3, 3, 1.0)

# Works identically for any topology
λ_dice = bloch_tb_spectrum(Dice, 4, 4, 1.0)
λ_ss   = bloch_tb_spectrum(ShastrySutherland, 3, 3, 1.0)
```

The builder reads `Lattice2D.get_unit_cell(Topology)`, which returns
the sublattice positions and a list of `Connection` objects
$(s_{\text{src}}, s_{\text{dst}}, d_x, d_y)$. No hand-derived
formulas are needed.

---

## Verification

| Test file | Scope | What is checked |
|-----------|-------|-----------------|
| `test_bloch_generic.jl` | All 4 lattices with `src/` code | Generic Bloch builder $=$ hardcoded formula |
| `test_graphene_tight_binding.jl` | Honeycomb | Closed form $=$ real-space ED |
| `test_kagome_tight_binding.jl` | Kagome | Closed form $=$ real-space ED, flat band count |
| `test_lieb_tight_binding.jl` | Lieb | Closed form $=$ real-space ED, flat band + M-point |
| `test_triangular_tight_binding.jl` | Triangular | Closed form $=$ real-space ED, frustration |
| `test_dice_tight_binding.jl` | Dice | Generic Bloch $=$ real-space ED |
| `test_shastry_sutherland_tight_binding.jl` | Shastry--Sutherland | Generic Bloch $=$ real-space ED |
| `test_unionjack_tight_binding.jl` | Union Jack | Generic Bloch $=$ real-space ED |

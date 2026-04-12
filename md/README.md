# QAtlas.jl — Internal Design Notes

Internal notes for QAtlas development. NOT user-facing documentation
(that lives in `docs/src/`).

## Current State (v0.12.6)

### What's implemented

**Models** (src/):
- IsingSquare: partition function (transfer matrix), T_c (Onsager), M(T) (Yang)
- TFIM: BdG spectrum (OBC, Infinite), thermal observables, dynamics, correlators
- Heisenberg1D: dimer spectrum, 4-site PBC, Bethe ansatz e₀
- Tight-binding: Graphene, Kagome, Lieb, Triangular (hardcoded Bloch)

**Universality classes** (src/universalities/):
- Universality{C} parametric type with dimension d keyword
- Ising (d=2 exact, d=3 bootstrap, MF), Percolation, Potts 3/4, KPZ, XY, Heisenberg, MF
- E8 mass ratios

**Test infrastructure** (test/util/):
- classical_partition.jl: brute-force 2^N enumeration
- tight_binding.jl: real-space TB Hamiltonian from bonds
- spinhalf_ed.jl: Kronecker embed, build_heisenberg, build_tfim, entanglement entropy
- bloch.jl: generic Bloch builder from unit cell

**Test suites** (850+ tests):
- standalone/: special values, scaling relations, Bethe ansatz
- verification/: Lattice2D cross-check (8 topologies), AD thermodynamics, gap closure, entanglement, disordered, universality cross-checks

**Documentation** (docs/src/, 46 pages):
- Zettelkasten calc/ notes (15 atomic derivation notes)
- Model pages, universality pages, verification philosophy, methods (Physical/Computational)

### Two API styles coexist

1. **Old style**: `Model{:TFIM}` + `Quantity{:energy}` + `OBC()` — used by TFIM, E8
2. **New style**: simple struct tags (`IsingSquare()`, `Universality(:Ising)`) — used by new models

No immediate plan to unify; both work. New additions should use the simple struct style.

## Pages in this directory

- [roadmap.md](roadmap.md) — future models, open issues, development priorities
- [documentation-design.md](documentation-design.md) — docs architecture blueprint
- [api.md](api.md) — API type hierarchy (mostly stable, see note about two styles)
- [testing.md](testing.md) — test philosophy (superseded by docs/src/verification/)
- [models/](models/) — model-specific notes (partially superseded by docs/src/models/)
- [universalities/](universalities/) — universality notes (superseded by docs/src/universalities/)

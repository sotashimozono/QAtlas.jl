# ─────────────────────────────────────────────────────────────────────────────
# test/util/classical_partition.jl
#
# Brute-force exact partition function for classical Ising models.
#
# Dependencies (expected to be `using`'d by the including test file):
#   Lattice2D  — provides `bonds(lat)`, `num_sites(lat)`
# ─────────────────────────────────────────────────────────────────────────────

"""
    exact_partition(lat, J::Float64, β::Float64) -> Float64

Exactly compute Z = Σ_σ exp(-β E(σ)) by enumerating all 2^N Ising spin
configurations on `lat`.

The energy is E(σ) = -J Σ_{b ∈ bonds(lat)} σ_{b.i} σ_{b.j}, using the
pre-computed bond list from `Lattice2D`. This adopts the same bond-counting
convention as `IsingSquare`'s transfer matrix: for PBC lattices where
Lx = 2 or Ly = 2, each unique physical bond appears twice in the bond list,
which doubles the effective coupling relative to an Lx,Ly ≥ 3 lattice.

Both `exact_partition` and `fetch(IsingSquare(), PartitionFunction(); ...)`
are internally consistent under this convention.

# Arguments
- `lat`: a finite `AbstractLattice` from `Lattice2D` (e.g., `build_lattice(Square, Lx, Ly)`)
- `J::Float64`: Ising coupling constant (J > 0 ferromagnetic)
- `β::Float64`: inverse temperature

# Complexity
O(2^N) — exact for any finite lattice; practical for N ≤ 20.
"""
function exact_partition(lat, J::Float64, β::Float64)
    N = num_sites(lat)
    bond_pairs = [(b.i, b.j) for b in bonds(lat)]
    Z = 0.0
    for σ_idx in 0:(2 ^ N - 1)
        σ = Int[((σ_idx >> j) & 1) == 1 ? 1 : -1 for j in 0:(N - 1)]
        E = -J * sum(σ[i] * σ[j] for (i, j) in bond_pairs)
        Z += exp(-β * E)
    end
    return Z
end

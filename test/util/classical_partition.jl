# ─────────────────────────────────────────────────────────────────────────────
# test/util/classical_partition.jl
#
# Brute-force exact partition function for classical Ising models on a square
# lattice with PBC in both directions.
#
# The bond list is constructed here (not taken from `Lattice2D.bonds`) to keep
# the brute-force enumeration in sync with the transfer-matrix PBC convention:
# for each site (i, j), the two bonds (i, j)–(i, (j mod Ly)+1) and
# (i, j)–((i mod Lx)+1, j) are enumerated. When a lattice dimension has length
# 2, each unique physical bond is enumerated twice, matching the
# `σ_j σ_{(j mod L)+1}` sum used inside `_ising_sq_transfer_matrix` and the
# two-edge cyclic trace `tr(T^Lx)`. This internal consistency is what makes
# `Z_bruteforce == Z_transfer-matrix` hold for small PBC lattices.
# ─────────────────────────────────────────────────────────────────────────────

"""
    square_pbc_bond_pairs(Lx, Ly) -> Vector{Tuple{Int,Int}}

Enumerate site-pair endpoints for the `Lx × Ly` square lattice with PBC in both
directions, in the transfer-matrix convention used by
`_ising_sq_transfer_matrix`: each site contributes one horizontal bond
`(i, j)–(i, (j mod Ly)+1)` and one vertical bond
`(i, j)–((i mod Lx)+1, j)`. When `Lx == 2` or `Ly == 2`, this enumerates each
unique physical bond twice — intentional, so that the brute-force sum matches
`tr(T^Lx)` element-for-element.

Sites are row-major indexed: `(i, j) ↦ (i - 1) * Ly + j`, `i ∈ 1:Lx`,
`j ∈ 1:Ly`.
"""
function square_pbc_bond_pairs(Lx::Int, Ly::Int)
    pairs = Vector{Tuple{Int,Int}}()
    sizehint!(pairs, 2 * Lx * Ly)
    idx(i, j) = (i - 1) * Ly + j
    for i in 1:Lx
        for j in 1:Ly
            a = idx(i, j)
            b_h = idx(i, (j % Ly) + 1)
            b_v = idx((i % Lx) + 1, j)
            push!(pairs, (a, b_h))
            push!(pairs, (a, b_v))
        end
    end
    return pairs
end

"""
    exact_partition(Lx, Ly, J, β) -> Float64

Exactly compute `Z = Σ_σ exp(-β E(σ))` by enumerating all `2^(Lx*Ly)` Ising
configurations. The energy is `E(σ) = -J Σ_{(a,b) ∈ pairs} σ_a σ_b`, using the
bond list from [`square_pbc_bond_pairs`](@ref).

Complexity: `O(2^(Lx*Ly))` — practical up to `N = Lx*Ly ≤ 20`.
"""
function exact_partition(Lx::Int, Ly::Int, J::Float64, β::Float64)
    N = Lx * Ly
    pairs = square_pbc_bond_pairs(Lx, Ly)
    Z = 0.0
    for σ_idx in 0:(2 ^ N - 1)
        σ = Int[((σ_idx >> k) & 1) == 1 ? 1 : -1 for k in 0:(N - 1)]
        E = -J * sum(σ[a] * σ[b] for (a, b) in pairs)
        Z += exp(-β * E)
    end
    return Z
end

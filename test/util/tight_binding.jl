# ─────────────────────────────────────────────────────────────────────────────
# test/util/tight_binding.jl
#
# Real-space single-particle tight-binding Hamiltonian constructor.
#
# Dependencies (expected to be `using`'d by the including test file):
#   Lattice2D  — provides `bonds(lat)`, `num_sites(lat)`
# ─────────────────────────────────────────────────────────────────────────────

"""
    build_tight_binding(lat, t::Real) -> Matrix{Float64}

Build the real-space single-particle tight-binding Hamiltonian

    H = -t Σ_{⟨i,j⟩} (c†_i c_j + c†_j c_i)

from the pre-computed bond list of `lat`. Each bond `b ∈ bonds(lat)`
contributes `H[b.i, b.j] -= t` and `H[b.j, b.i] -= t`.

This follows the same bond-counting convention as `Lattice2D.bonds`:
for small periodic systems where an edge is traversed in more than one
direction by the connection list (e.g. Lx = 2 square), the corresponding
matrix entries are naturally accumulated, consistent with the physical
"double-bond" of a 2-site periodic ring.

Returns a dense `Matrix{Float64}` of size `num_sites(lat) × num_sites(lat)`.
"""
function build_tight_binding(lat, t::Real)
    N = num_sites(lat)
    H = zeros(Float64, N, N)
    for b in bonds(lat)
        H[b.i, b.j] -= t
        H[b.j, b.i] -= t
    end
    return H
end

# ─────────────────────────────────────────────────────────────────────────────
# Graphene — nearest-neighbor tight-binding on the honeycomb lattice
#
# Hamiltonian (single-particle, spinless):
#   H = -t Σ_{⟨i,j⟩ ∈ A-B} (c†_i c_j + c†_j c_i)
#
# The Bloch Hamiltonian on a bipartite (A, B) basis is:
#
#   H(k) = [  0      f(k)  ]
#          [ f(k)*   0     ]
#
# with f(k) = -t [exp(i k·δ₁) + exp(i k·δ₂) + exp(i k·δ₃)], where δᵢ are
# the three nearest-neighbor displacement vectors from an A site to its
# three B neighbors. The two bands are ±|f(k)|.
#
# Using the unit-cell basis (a₁, a₂) adopted by `Lattice2D`'s honeycomb
# topology (a₁ = (√3, 0), a₂ = (√3/2, 3/2), δ's derived accordingly):
#
#   |f(k)|² = t² [3 + 2cos(k·a₁) + 2cos(k·a₂) + 2cos(k·(a₂−a₁))]
#
# References:
#   P. R. Wallace, "The Band Theory of Graphite", Phys. Rev. 71, 622 (1947).
#   A. H. Castro Neto et al., "The electronic properties of graphene",
#     Rev. Mod. Phys. 81, 109 (2009).
# ─────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# Dispatch tags
# ═══════════════════════════════════════════════════════════════════════════════

"""
    Graphene

Dispatch tag for the nearest-neighbor tight-binding model on the
honeycomb lattice (spinless, single-orbital).

Hamiltonian: H = -t Σ_{⟨i,j⟩} (c†_i c_j + h.c.)
"""
struct Graphene end

"""
    TightBindingSpectrum

Dispatch tag for the single-particle spectrum of a tight-binding model.
"""
struct TightBindingSpectrum end

# ═══════════════════════════════════════════════════════════════════════════════
# fetch: closed-form Bloch spectrum for Lx × Ly honeycomb PBC
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::Graphene, ::TightBindingSpectrum; Lx, Ly, t=1.0) -> Vector{Float64}

Sorted single-particle spectrum of the nearest-neighbor tight-binding
Hamiltonian on an Lx × Ly honeycomb lattice with periodic boundary
conditions in both directions.

The spectrum is obtained in closed form by diagonalizing the 2×2 Bloch
Hamiltonian at the `Lx·Ly` allowed momenta

    k_{mn} = (m/Lx) b₁ + (n/Ly) b₂,   m ∈ {0,…,Lx−1}, n ∈ {0,…,Ly−1}

where (b₁, b₂) is the reciprocal basis dual to the real-space
basis (a₁, a₂) used by `Lattice2D`'s `Honeycomb` topology.

Explicitly, using k·a₁ = 2π m/Lx and k·a₂ = 2π n/Ly,

    E_{mn,±} = ± t · √(3 + 2cos(2π m/Lx) + 2cos(2π n/Ly) + 2cos(2π(n/Ly − m/Lx)))

The chiral (sublattice) symmetry of the bipartite honeycomb lattice
guarantees that the spectrum is symmetric about zero. A small negative
floor is applied inside the square root to absorb rounding error at
Dirac points where the argument vanishes exactly.

# Arguments
- `Lx::Int`: number of unit cells in the first lattice direction
- `Ly::Int`: number of unit cells in the second lattice direction
- `t::Real`: nearest-neighbor hopping amplitude (default 1.0)

# Return
A sorted `Vector{Float64}` of length `2·Lx·Ly`.

# References
    P. R. Wallace, Phys. Rev. 71, 622 (1947).
    A. H. Castro Neto et al., Rev. Mod. Phys. 81, 109 (2009).
"""
function fetch(::Graphene, ::TightBindingSpectrum; Lx::Int, Ly::Int, t::Real=1.0)
    eigs = Float64[]
    sizehint!(eigs, 2 * Lx * Ly)
    for m in 0:(Lx - 1), n in 0:(Ly - 1)
        θ1 = 2π * m / Lx
        θ2 = 2π * n / Ly
        # |f(k)|² / t² = 3 + 2cos(θ1) + 2cos(θ2) + 2cos(θ2 - θ1)
        val = 3 + 2cos(θ1) + 2cos(θ2) + 2cos(θ2 - θ1)
        energy = abs(t) * sqrt(max(val, 0.0))
        push!(eigs, -energy)
        push!(eigs, energy)
    end
    return sort!(eigs)
end

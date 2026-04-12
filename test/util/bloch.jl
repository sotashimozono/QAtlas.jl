# ─────────────────────────────────────────────────────────────────────────────
# test/util/bloch.jl
#
# Generic Bloch Hamiltonian builder for nearest-neighbor tight-binding
# on any Lattice2D periodic topology.  Given a topology type (Square,
# Honeycomb, Kagome, Lieb, Triangular, Dice, …) and a system size
# (Lx, Ly), this util constructs the n_sub × n_sub Bloch Hamiltonian
#
#   H_{αβ}(k) = -t Σ_{connections c: β→α at offset (dx,dy)}  exp(i(dx θ₁ + dy θ₂))
#             + h.c.
#
# at each of the Lx·Ly allowed k-points (θ₁ = 2π m/Lx, θ₂ = 2π n/Ly),
# diagonalizes it, and returns the full sorted spectrum.
#
# This provides an INDEPENDENT cross-check of the hand-derived Bloch
# formulas in src/models/quantum/tightbinding/*.jl: if the generic
# builder agrees with the hardcoded formula, both are correct.
#
# Dependencies (expected to be `using`'d by the including test file):
#   Lattice2D  — `get_unit_cell`, `Connection`, topology types
#   LinearAlgebra — `eigvals`, `Hermitian`
# ─────────────────────────────────────────────────────────────────────────────

"""
    bloch_tb_spectrum(Topology, Lx, Ly, t) -> Vector{Float64}

Compute the sorted single-particle spectrum of the nearest-neighbor
tight-binding Hamiltonian on an `Lx × Ly` periodic lattice of topology
`Topology` by diagonalizing the `n_sub × n_sub` Bloch Hamiltonian at
each allowed momentum.

This function reads the unit-cell definition (sublattice positions,
connection list) directly from `Lattice2D.get_unit_cell(Topology)` and
requires no hand-derived formulas — it is fully generic over any 2D
periodic topology that `Lattice2D` supports.

# Arguments
- `Topology`: a Lattice2D topology type (e.g., `Honeycomb`, `Kagome`)
- `Lx::Int`, `Ly::Int`: number of unit cells in each direction
- `t::Real`: nearest-neighbor hopping amplitude

# Return
A sorted `Vector{Float64}` of length `n_sub · Lx · Ly`.
"""
function bloch_tb_spectrum(
    Topology::Type{<:AbstractTopology{2}}, Lx::Int, Ly::Int, t::Real
)
    uc = get_unit_cell(Topology)
    nsub = length(uc.sublattice_positions)
    eigs = Float64[]
    sizehint!(eigs, nsub * Lx * Ly)

    for m in 0:(Lx - 1), n in 0:(Ly - 1)
        θ1 = 2π * m / Lx
        θ2 = 2π * n / Ly

        H = zeros(ComplexF64, nsub, nsub)
        for c in uc.connections
            phase = exp(im * (c.dx * θ1 + c.dy * θ2))
            H[c.dst_sub, c.src_sub] -= t * phase
            H[c.src_sub, c.dst_sub] -= t * conj(phase)
        end

        append!(eigs, eigvals(Hermitian(H)))
    end

    return sort!(eigs)
end

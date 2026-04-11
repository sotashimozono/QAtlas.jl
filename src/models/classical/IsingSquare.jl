# ─────────────────────────────────────────────────────────────────────────────
# Classical 2D Ising model on the square lattice — exact partition function
#
# Hamiltonian:
#   H = -J Σ_{⟨i,j⟩} σᵢ σⱼ,   σᵢ ∈ {−1, +1}
#
# Solved exactly via the transfer-matrix method. For an Lx × Ly torus
# the partition function is Z = Tr(T^Lx) where T is the 2^Ly × 2^Ly
# transfer matrix.
#
# References:
#   L. Onsager, "Crystal Statistics. I. A Two-Dimensional Model with an
#   Order-Disorder Transition", Phys. Rev. 65, 117 (1944).
#   B. M. McCoy and T. T. Wu, "The Two-Dimensional Ising Model",
#   Harvard University Press (1973).
# ─────────────────────────────────────────────────────────────────────────────

using LinearAlgebra: Symmetric, eigvals

# ═══════════════════════════════════════════════════════════════════════════════
# Dispatch tags
# ═══════════════════════════════════════════════════════════════════════════════

"""
    IsingSquare

Dispatch tag for the classical 2D Ising model on a square lattice with
periodic boundary conditions (PBC) in both directions.

Hamiltonian: H = -J Σ_{⟨i,j⟩} σᵢ σⱼ, σᵢ ∈ {-1, +1}.

See also: [`PartitionFunction`](@ref), [`fetch(::IsingSquare, ::PartitionFunction)`](@ref).
"""
struct IsingSquare end

"""
    PartitionFunction

Dispatch tag for the canonical partition function Z = Σ_σ exp(-β H(σ)).

See also: [`IsingSquare`](@ref).
"""
struct PartitionFunction end

# ═══════════════════════════════════════════════════════════════════════════════
# Internal: transfer matrix for a row of Ly spins (PBC in y)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _ising_sq_transfer_matrix(Ly, β, J) -> Matrix{Float64}

Return the 2^Ly × 2^Ly symmetric transfer matrix for a row of `Ly` Ising
spins with PBC along the row direction.

The symmetric (Hermitian) form ensures real eigenvalues and is defined by:

    T_{σ,σ'} = exp(βJ/2 · Eₕ(σ)) · exp(βJ · Eᵥ(σ,σ')) · exp(βJ/2 · Eₕ(σ'))

where
    Eₕ(σ) = Σⱼ₌₁^Ly σⱼ σ_{(j mod Ly)+1}   (horizontal bonds within row, PBC)
    Eᵥ(σ,σ') = Σⱼ₌₁^Ly σⱼ σ'ⱼ             (vertical bonds between rows)

Note: for Ly = 2, each horizontal bond is counted twice by the PBC sum
(σ₁σ₂ + σ₂σ₁ = 2σ₁σ₂). This is consistent with the convention used by
`Lattice2D`'s bond list, which also lists each bond of a 2-site periodic
row twice.
"""
function _ising_sq_transfer_matrix(Ly::Int, β::Float64, J::Float64)
    dim = 2^Ly

    # Precompute spin vectors for each row state index (0-based binary encoding)
    spins_of = Vector{Vector{Int}}(undef, dim)
    for σ_idx in 0:(dim - 1)
        spins_of[σ_idx + 1] = Int[((σ_idx >> j) & 1) == 1 ? 1 : -1 for j in 0:(Ly - 1)]
    end

    # Diagonal weight: exp(βJ/2 · Eₕ(σ))
    h_weight = Vector{Float64}(undef, dim)
    for σ_idx in 0:(dim - 1)
        σ = spins_of[σ_idx + 1]
        e_h = sum(σ[j] * σ[(j % Ly) + 1] for j in 1:Ly)
        h_weight[σ_idx + 1] = exp(β * J / 2 * e_h)
    end

    # Build transfer matrix
    T = Matrix{Float64}(undef, dim, dim)
    for σ_idx in 0:(dim - 1)
        σ = spins_of[σ_idx + 1]
        wσ = h_weight[σ_idx + 1]
        for σp_idx in 0:(dim - 1)
            σp = spins_of[σp_idx + 1]
            e_v = sum(σ[j] * σp[j] for j in 1:Ly)
            T[σ_idx + 1, σp_idx + 1] = wσ * exp(β * J * e_v) * h_weight[σp_idx + 1]
        end
    end
    return T
end

# ═══════════════════════════════════════════════════════════════════════════════
# fetch: partition function via transfer matrix
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::IsingSquare, ::PartitionFunction; Lx, Ly, β, J=1.0) -> Float64

Exact partition function Z = Tr(T^Lx) for the classical 2D Ising model on an
Lx × Ly square lattice with periodic boundary conditions in both directions.

The transfer matrix T acts along the Lx direction (row-to-row transfer), with
each row containing Ly spins and PBC along the Ly direction. Explicitly:

    Z = Σᵢ λᵢ^Lx

where λᵢ are the eigenvalues of the 2^Ly × 2^Ly symmetric transfer matrix.

# Bond-counting convention

This function adopts the same bond-counting convention as `Lattice2D`'s
pre-computed bond list: for small PBC systems where Lx = 2 or Ly = 2, each
unique physical bond appears twice in the enumeration. Both the transfer-matrix
result and the brute-force enumeration via `Lattice2D.bonds` are internally
consistent under this convention.

# Special values

- β = 0 (any Lx, Ly, J): Z = 2^(Lx·Ly)  — all configurations equally weighted
- J = 0 (any β, Lx, Ly): Z = 2^(Lx·Ly)  — no interactions, same as β = 0

# Arguments
- `Lx::Int`: number of rows (transfer direction)
- `Ly::Int`: number of columns (row length, PBC)
- `β::Float64`: inverse temperature (β = 1/(k_B T))
- `J::Float64`: Ising coupling constant (default 1.0; J > 0 ferromagnetic)

# References
    L. Onsager, Phys. Rev. 65, 117 (1944).
    B. M. McCoy and T. T. Wu, "The Two-Dimensional Ising Model" (1973).
"""
function fetch(
    ::IsingSquare, ::PartitionFunction; Lx::Int, Ly::Int, β::Float64, J::Float64=1.0
)
    T = _ising_sq_transfer_matrix(Ly, β, J)
    λs = eigvals(Symmetric(T))
    return sum(λ^Lx for λ in λs)
end

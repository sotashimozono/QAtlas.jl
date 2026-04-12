# ─────────────────────────────────────────────────────────────────────────────
# test/util/spinhalf_ed.jl
#
# Dense exact-diagonalization helpers for spin-1/2 many-body Hamiltonians.
# Each site carries a 2-dimensional local Hilbert space; the total space is
# 2^N-dimensional. Basis ordering: Julia `kron` convention, i.e. site 1 is
# the outermost (most significant) factor in the tensor product.
#
# Dependencies (expected to be `using`'d by the including test file):
#   LinearAlgebra  — `I`, `kron`
#   Lattice2D      — `num_sites(lat)`, `bonds(lat)`
# ─────────────────────────────────────────────────────────────────────────────

"""
    embed_two_site(A, B, i, j, N) -> Matrix{Float64}

Embed a tensor-product two-site operator `A_i ⊗ B_j` into the full 2^N
Hilbert space of `N` spin-1/2 sites. Both `A` and `B` are 2×2 operators
acting on a single site; all other sites carry identity.

If `i > j`, the arguments are swapped (including the operators) so that
the smaller index is handled first. This is correct for Hermitian
bilinears such as `S_i · S_j` which are symmetric under `i ↔ j`:
`(Sp, Sm, i, j)` with `i > j` is silently reinterpreted as
`(Sm, Sp, j, i)`, which is the same bond contribution.
"""
function embed_two_site(A::AbstractMatrix, B::AbstractMatrix, i::Int, j::Int, N::Int)
    @assert 1 <= i <= N && 1 <= j <= N "site indices out of range"
    @assert i != j "two-site embed requires distinct sites"
    if i > j
        i, j = j, i
        A, B = B, A
    end
    left = Matrix{Float64}(I, 2^(i - 1), 2^(i - 1))
    middle = Matrix{Float64}(I, 2^(j - i - 1), 2^(j - i - 1))
    right = Matrix{Float64}(I, 2^(N - j), 2^(N - j))
    return kron(kron(kron(kron(left, A), middle), B), right)
end

"""
    embed_single_site(A, i, N) -> Matrix{Float64}

Embed a single-site operator `A` (2×2) on site `i` into the full 2^N
Hilbert space of `N` spin-1/2 sites. All other sites carry identity.
"""
function embed_single_site(A::AbstractMatrix, i::Int, N::Int)
    @assert 1 <= i <= N "site index out of range"
    left = Matrix{Float64}(I, 2^(i - 1), 2^(i - 1))
    right = Matrix{Float64}(I, 2^(N - i), 2^(N - i))
    return kron(kron(left, A), right)
end

"""
    build_tfim(lat, J, h) -> Matrix{Float64}

Construct the full 2^N × 2^N dense Hamiltonian of the 1D transverse-field
Ising model

    H = -J Σ_{⟨i,j⟩} σᶻ_i σᶻ_j  -  h Σ_i σˣ_i

on the lattice `lat`. The Ising interaction iterates over `bonds(lat)`;
the transverse field acts on all `num_sites(lat)` sites.

Uses Pauli matrices (σᶻ, σˣ) directly — the resulting matrix is real
symmetric.
"""
function build_tfim(lat, J::Real, h::Real)
    N = num_sites(lat)
    dim = 2^N
    H = zeros(Float64, dim, dim)

    σz = [1.0 0.0; 0.0 -1.0]
    σx = [0.0 1.0; 1.0 0.0]

    for b in bonds(lat)
        H .-= J * embed_two_site(σz, σz, b.i, b.j, N)
    end
    for i in 1:N
        H .-= h * embed_single_site(σx, i, N)
    end
    return H
end

"""
    build_spinhalf_heisenberg(lat, J) -> Matrix{Float64}

Construct the full 2^N × 2^N dense Hamiltonian of the spin-1/2
antiferromagnetic Heisenberg model

    H = J Σ_{⟨i,j⟩} S_i · S_j
      = J Σ_{⟨i,j⟩} [ Sᶻ_i Sᶻ_j + (1/2)(S⁺_i S⁻_j + S⁻_i S⁺_j) ]

on the lattice `lat`, iterating over `bonds(lat)`. Uses ladder operators
to avoid complex arithmetic — the resulting matrix is real symmetric.

Each bond contributes once; no de-duplication is performed. Small PBC
systems where `Lattice2D` lists an edge multiply (e.g. Lx = 2 ring)
accumulate the multiplicity, matching `Lattice2DMonteCarlo`'s classical
convention.
"""
function build_spinhalf_heisenberg(lat, J::Real)
    N = num_sites(lat)
    dim = 2^N
    H = zeros(Float64, dim, dim)

    # Spin-1/2 ladder operators and Sz in the {|↑⟩, |↓⟩} basis.
    Sp = [0.0 1.0; 0.0 0.0]
    Sm = [0.0 0.0; 1.0 0.0]
    Sz = [0.5 0.0; 0.0 -0.5]

    for b in bonds(lat)
        i, j = b.i, b.j
        H .+= J * embed_two_site(Sz, Sz, i, j, N)
        H .+= (J / 2) * embed_two_site(Sp, Sm, i, j, N)
        H .+= (J / 2) * embed_two_site(Sm, Sp, i, j, N)
    end
    return H
end

# ═══════════════════════════════════════════════════════════════════════════════
# Entanglement entropy
# ═══════════════════════════════════════════════════════════════════════════════

"""
    reduced_density_matrix(ψ, l, N) -> Matrix{Float64}

Compute the reduced density matrix ρ_A = tr_B |ψ⟩⟨ψ| for the first
`l` sites (subsystem A) of an `N`-site spin-1/2 state `ψ`.

The state vector `ψ` has length 2^N (full Hilbert space). The result
is a 2^l × 2^l density matrix.

The Hilbert space is decomposed as H = H_A ⊗ H_B with dim_A = 2^l,
dim_B = 2^(N−l). The state is reshaped into a (dim_A, dim_B) matrix
and ρ_A = M M†.
"""
function reduced_density_matrix(ψ::AbstractVector, l::Int, N::Int)
    @assert length(ψ) == 2^N
    @assert 1 <= l < N
    dim_A = 2^l
    dim_B = 2^(N - l)
    M = reshape(ψ, dim_A, dim_B)
    return M * M'
end

"""
    entanglement_entropy(ψ, l, N) -> Float64

Von Neumann entanglement entropy S(l) = −tr(ρ_A ln ρ_A) for the
bipartition of an N-site state ψ into subsystem A (sites 1:l) and
subsystem B (sites l+1:N).

For a critical 1D system of length N with OBC, the Calabrese-Cardy
formula predicts:

    S(l) = (c/3) ln[(2N/π) sin(πl/N)] + const

where c is the central charge. This formula enables extraction of c
from a single ground-state wavefunction by fitting S(l) vs l.

# References
    P. Calabrese, J. Cardy, "Entanglement entropy and quantum field
    theory", J. Stat. Mech. 0406, P06002 (2004).
"""
function entanglement_entropy(ψ::AbstractVector, l::Int, N::Int)
    ρ = reduced_density_matrix(ψ, l, N)
    λs = eigvals(Hermitian(ρ))
    S = 0.0
    for λ in λs
        if λ > 1e-15
            S -= λ * log(λ)
        end
    end
    return S
end

# ─────────────────────────────────────────────────────────────────────────────
# test/util/sparse_ed.jl
#
# Sparse-matrix exact-diagonalization helpers for spin-1/2 many-body
# Hamiltonians.  Same basis convention as spinhalf_ed.jl (Julia `kron`,
# site 1 outermost) but the embed functions return `SparseMatrixCSC`,
# and `ground_state_krylov` extracts only the lowest eigenpair via
# Lanczos, avoiding the O(D³) cost of a full dense diagonalisation.
#
# Intended for tests that build many disorder realisations of a large
# Hilbert space and need only the ground state — dense ED quickly
# becomes the test-suite bottleneck.
#
# Dependencies (expected to be `using`'d by the including test file):
#   LinearAlgebra  — `I`
#   SparseArrays   — `SparseMatrixCSC`, `sparse`, `spzeros`
#   KrylovKit      — `eigsolve`
#   Random         — `AbstractRNG` (optional, for deterministic start vectors)
# ─────────────────────────────────────────────────────────────────────────────

"""
    embed_two_site_sparse(A, B, i, j, N) -> SparseMatrixCSC{Float64}

Sparse analogue of `embed_two_site` (see `spinhalf_ed.jl`): embeds
`A_i ⊗ B_j` into the full 2^N spin-1/2 Hilbert space.  For typical
Pauli-like operators the result has O(2^N) non-zeros rather than the
dense 2^N × 2^N.
"""
function embed_two_site_sparse(A::AbstractMatrix, B::AbstractMatrix, i::Int, j::Int, N::Int)
    @assert 1 <= i <= N && 1 <= j <= N "site indices out of range"
    @assert i != j "two-site embed requires distinct sites"
    if i > j
        i, j = j, i
        A, B = B, A
    end
    As = sparse(A)
    Bs = sparse(B)
    left = sparse(I, 2^(i - 1), 2^(i - 1))
    middle = sparse(I, 2^(j - i - 1), 2^(j - i - 1))
    right = sparse(I, 2^(N - j), 2^(N - j))
    return kron(kron(kron(kron(left, As), middle), Bs), right)
end

"""
    embed_single_site_sparse(A, i, N) -> SparseMatrixCSC{Float64}

Sparse analogue of `embed_single_site` (see `spinhalf_ed.jl`).
"""
function embed_single_site_sparse(A::AbstractMatrix, i::Int, N::Int)
    @assert 1 <= i <= N "site index out of range"
    As = sparse(A)
    left = sparse(I, 2^(i - 1), 2^(i - 1))
    right = sparse(I, 2^(N - i), 2^(N - i))
    return kron(kron(left, As), right)
end

"""
    build_tfim_sparse(lat, J, h) -> SparseMatrixCSC{Float64}

Sparse analogue of `build_tfim` from `spinhalf_ed.jl`.  Memory footprint
is O(N · 2^N) instead of O(4^N), enabling N ≳ 14 entanglement tests
without blowing out dense eigen() allocations.
"""
function build_tfim_sparse(lat, J::Real, h::Real)
    N = num_sites(lat)
    dim = 2^N
    H = spzeros(Float64, dim, dim)
    σz = [1.0 0.0; 0.0 -1.0]
    σx = [0.0 1.0; 1.0 0.0]
    for b in bonds(lat)
        H -= J * embed_two_site_sparse(σz, σz, b.i, b.j, N)
    end
    for i in 1:N
        H -= h * embed_single_site_sparse(σx, i, N)
    end
    return H
end

"""
    build_spinhalf_heisenberg_sparse(lat, J) -> SparseMatrixCSC{Float64}

Sparse analogue of `build_spinhalf_heisenberg` from `spinhalf_ed.jl`.
Same Hamiltonian, O(N · 2^N) nonzeros.
"""
function build_spinhalf_heisenberg_sparse(lat, J::Real)
    N = num_sites(lat)
    dim = 2^N
    H = spzeros(Float64, dim, dim)
    Sp = [0.0 1.0; 0.0 0.0]
    Sm = [0.0 0.0; 1.0 0.0]
    Sz = [0.5 0.0; 0.0 -0.5]
    for b in bonds(lat)
        H += J * embed_two_site_sparse(Sz, Sz, b.i, b.j, N)
        H += (J / 2) * embed_two_site_sparse(Sp, Sm, b.i, b.j, N)
        H += (J / 2) * embed_two_site_sparse(Sm, Sp, b.i, b.j, N)
    end
    return H
end

"""
    ground_state_krylov(H; rng = nothing, tol = 1e-10, krylovdim = 30) -> Vector{Float64}

Return the ground-state eigenvector of a symmetric / Hermitian linear
operator `H` via `KrylovKit.eigsolve` (Lanczos for real symmetric input).
Only the lowest eigenpair is computed, so for sparse `H` the cost is
O(nnz(H) · k) per Lanczos iteration rather than O(D³).

Passing an `AbstractRNG` makes the random start vector deterministic,
which keeps disorder-averaged tests reproducible across runs.
"""
function ground_state_krylov(
    H; rng::Union{AbstractRNG,Nothing}=nothing, tol::Float64=1e-10, krylovdim::Int=30
)
    D = size(H, 1)
    x₀ = rng === nothing ? randn(D) : randn(rng, D)
    vals, vecs, info = eigsolve(
        H, x₀, 1, :SR; issymmetric=true, tol=tol, krylovdim=krylovdim
    )
    info.converged < 1 && @warn "KrylovKit ground state failed to converge" info
    return vecs[1]
end

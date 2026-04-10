# ─────────────────────────────────────────────────────────────────────────────
# Transverse Field Ising Model — exact real-time correlation functions
#
# Hamiltonian (OBC, N sites):
#
#   H = -J Σ_i σᶻ_i σᶻ_{i+1}  -  h Σ_i σˣ_i
#
# Jordan–Wigner mapping (standard): σˣ diagonal, σᶻ nonlocal.  Introducing
# Majorana operators
#
#   γ_{2i-1} = c_i + c_i†
#   γ_{2i}   = i (c_i† - c_i)            {γ_a, γ_b} = 2 δ_{ab}
#
# one finds
#
#   σˣ_i                = -i γ_{2i-1} γ_{2i}
#   σᶻ_i σᶻ_{i+1}       = -i γ_{2i}   γ_{2i+1}
#   σᶻ_i                = (-i)^{i-1} γ_1 γ_2 … γ_{2i-2} γ_{2i-1}  (JW string,
#                                                                  consecutive
#                                                                  Majoranas)
#
# (the single-operator expression contains an odd number of Majoranas — its
# expectation value in a Gaussian state vanishes, so ⟨σᶻ⟩_GS = 0 on OBC.)
#
# The Hamiltonian takes the quadratic Majorana form
#
#   H = (i/4) Σ_{ab} h_{ab} γ_a γ_b
#
# with h a 2N × 2N real antisymmetric matrix whose only non-zero entries
# (upper triangle) are
#
#   h_{2i-1, 2i} = 2h    (i = 1 … N)
#   h_{2i  , 2i+1} = 2J  (i = 1 … N-1)
#
# The Heisenberg-picture evolution reads
#
#   γ(t) = exp(h t) γ(0),          exp(h t) ∈ SO(2N)
#
# and the ground-state Majorana 2-point function is
#
#   ⟨γ_a γ_b⟩_GS = δ_{ab} + i Σ_{ab},       Σ = -i · sign(i h)
#
# (derivation: in the canonical basis where h is block-diagonal with 2×2
# blocks [0 Λ; -Λ 0], the ground-state block is `Σ̃_block = [0 1; -1 0]`,
# while a direct computation gives `i sign(i h)|block = [0 -1; 1 0]`, hence
# the overall minus sign).
# Higher-point correlators reduce to a Pfaffian via Wick's theorem.
#
# The quantities exposed to `fetch` are
#
#   :sz_sz_correlation  — ⟨σᶻ_i(t) σᶻ_j(0)⟩_GS          (Pfaffian, size 2(i+j)-2)
#   :sx_sx_correlation  — ⟨σˣ_i(t) σˣ_j(0)⟩_GS          (4×4 Pfaffian)
#   :sz_sz_spreading    — matrix C[it, ix] = ⟨σᶻ_{ix}(t) σᶻ_center(0)⟩
#                          over all sites × supplied `times`
# ─────────────────────────────────────────────────────────────────────────────

using LinearAlgebra: Hermitian, eigen, Diagonal, I

# ═══════════════════════════════════════════════════════════════════════════════
# Building blocks
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _majorana_ham(N, J, h) -> Matrix{Float64}

Return the 2N × 2N real antisymmetric Majorana Hamiltonian matrix `h`
of the OBC TFIM such that `H = (i/4) Σ h_{ab} γ_a γ_b`.
"""
function _majorana_ham(N::Int, J::Float64, h::Float64)
    M = zeros(Float64, 2N, 2N)
    @inbounds for i in 1:N
        M[2i - 1, 2i] = 2h
        M[2i, 2i - 1] = -2h
    end
    @inbounds for i in 1:(N - 1)
        M[2i, 2i + 1] = 2J
        M[2i + 1, 2i] = -2J
    end
    return M
end

"""
    _majorana_covariance_gs(h::AbstractMatrix) -> Matrix{Float64}

Ground-state Majorana covariance Σ of a quadratic Hamiltonian whose
Majorana matrix is `h`.  Σ is real antisymmetric, with
`⟨γ_a γ_b⟩_GS = δ_{ab} + i Σ_{ab}`.  The formula `Σ = -i · sign(i h)` is
used, evaluated through an eigendecomposition of the Hermitian matrix
`i h`.
"""
function _majorana_covariance_gs(h::AbstractMatrix{<:Real})
    M = im .* h                                # Hermitian 2N × 2N
    F = eigen(Hermitian((M + M') / 2))
    s = [λ > 0 ? 1.0 : (λ < 0 ? -1.0 : 0.0) for λ in F.values]
    sM = F.vectors * Diagonal(s) * F.vectors'
    Σ = real(-im .* sM)
    # Symmetrise to strict antisymmetry to cancel round-off drift.
    return (Σ - Σ') / 2
end

"""
    _majorana_evolution(h::AbstractMatrix, t::Real) -> Matrix{Float64}

`exp(h * t)` for a real antisymmetric `h`, returned as a real
orthogonal matrix.  Tiny imaginary noise from the matrix exponential
is projected away.
"""
function _majorana_evolution(h::AbstractMatrix{<:Real}, t::Real)
    if t == 0
        return Matrix{Float64}(I, size(h)...)
    end
    return real(exp(h .* t))
end

# Majorana index list for σᶻ_k = (-i)^{k-1} γ_1 γ_2 … γ_{2k-2} γ_{2k-1}.
# (Equivalent to σˣ_k in the standard Kitaev convention; QAtlas's H uses
# Z-Z coupling, X-field, which is unitarily related by a Hadamard rotation.)
_sz_majorana_indices(k::Int) = collect(1:(2k - 1))

# For σˣ_k = -i γ_{2k-1} γ_{2k}, the two Majoranas are (2k-1, 2k).
_sx_majorana_indices(k::Int) = Int[2k - 1, 2k]

"""
    _build_wick_matrix(idx_t, idx_0, Σ, R, RΣ) -> Matrix{ComplexF64}

Assemble the Wick contraction matrix `F` for a product

    γ_{idx_t[1]}(t) γ_{idx_t[2]}(t) … γ_{idx_0[1]}(0) γ_{idx_0[2]}(0) …

from the ground-state Majorana covariance `Σ`, the time evolution
matrix `R = exp(h t)`, and the precomputed product `RΣ = R * Σ`.
`F` is complex antisymmetric; its Pfaffian gives `⟨…⟩_GS`.
"""
function _build_wick_matrix(
    idx_t::Vector{Int},
    idx_0::Vector{Int},
    Σ::AbstractMatrix,
    R::AbstractMatrix,
    RΣ::AbstractMatrix,
)
    na = length(idx_t)
    nb = length(idx_0)
    M = na + nb
    F = zeros(ComplexF64, M, M)

    # (t, t) block
    @inbounds for k in 1:na, l in (k + 1):na
        a = idx_t[k]
        b = idx_t[l]
        val = (a == b ? 1.0 : 0.0) + im * Σ[a, b]
        F[k, l] = val
        F[l, k] = -val
    end

    # (t, 0) cross block
    @inbounds for k in 1:na, l in 1:nb
        a = idx_t[k]
        b = idx_0[l]
        val = R[a, b] + im * RΣ[a, b]
        F[k, na + l] = val
        F[na + l, k] = -val
    end

    # (0, 0) block
    @inbounds for k in 1:nb, l in (k + 1):nb
        a = idx_0[k]
        b = idx_0[l]
        val = (a == b ? 1.0 : 0.0) + im * Σ[a, b]
        F[na + k, na + l] = val
        F[na + l, na + k] = -val
    end

    return F
end

# ═══════════════════════════════════════════════════════════════════════════════
# σᶻ σᶻ and σˣ σˣ real-time correlation functions
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _sz_sz_corr(N, J, h, i, j, t) -> ComplexF64

`⟨GS| σᶻ_i(t) σᶻ_j(0) |GS⟩` for the OBC TFIM.

Implementation: Wick / Pfaffian over the `(2i-1) + (2j-1) = 2(i+j)-2`
Majorana operators constituting the product.  The overall phase is
`i^{i+j-2}` (from the two `-i^{k-1}` factors in each `σᶻ_k`).
"""
function _sz_sz_corr(N::Int, J::Float64, h::Float64, i::Int, j::Int, t::Real)
    (1 ≤ i ≤ N && 1 ≤ j ≤ N) || throw(ArgumentError("site indices out of range"))
    hmat = _majorana_ham(N, J, h)
    Σ = _majorana_covariance_gs(hmat)
    R = _majorana_evolution(hmat, t)
    return _sz_sz_corr_from_cached(Σ, R, i, j)
end

function _sz_sz_corr_from_cached(Σ::AbstractMatrix, R::AbstractMatrix, i::Int, j::Int)
    RΣ = R * Σ
    idx_t = _sz_majorana_indices(i)
    idx_0 = _sz_majorana_indices(j)
    F = _build_wick_matrix(idx_t, idx_0, Σ, R, RΣ)
    pf = pfaffian(F)
    # σᶻ_k = (-i)^{k-1} γ_1 … γ_{2k-1}, so the prefactor of the product is
    # (-i)^{i+j-2}.
    return ((-im)^(i + j - 2)) * pf
end

"""
    _sx_sx_corr(N, J, h, i, j, t) -> ComplexF64

`⟨GS| σˣ_i(t) σˣ_j(0) |GS⟩` for the OBC TFIM.  Reduces to a 4×4
Pfaffian since `σˣ_k = -i γ_{2k-1} γ_{2k}`.
"""
function _sx_sx_corr(N::Int, J::Float64, h::Float64, i::Int, j::Int, t::Real)
    (1 ≤ i ≤ N && 1 ≤ j ≤ N) || throw(ArgumentError("site indices out of range"))
    hmat = _majorana_ham(N, J, h)
    Σ = _majorana_covariance_gs(hmat)
    R = _majorana_evolution(hmat, t)
    return _sx_sx_corr_from_cached(Σ, R, i, j)
end

function _sx_sx_corr_from_cached(Σ::AbstractMatrix, R::AbstractMatrix, i::Int, j::Int)
    RΣ = R * Σ
    idx_t = _sx_majorana_indices(i)
    idx_0 = _sx_majorana_indices(j)
    F = _build_wick_matrix(idx_t, idx_0, Σ, R, RΣ)
    pf = pfaffian(F)
    # σˣ carries factor (-i)·(-i) = -1 from the two -i γγ prefactors.
    return -pf
end

# ═══════════════════════════════════════════════════════════════════════════════
# Single-operator expectations
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _sx_expect(Σ, i) -> Float64

`⟨σˣ_i⟩_GS = Σ[2i-1, 2i]`, since `σˣ_i = -i γ_{2i-1} γ_{2i}`.
"""
_sx_expect(Σ::AbstractMatrix, i::Int) = Σ[2i - 1, 2i]

# ⟨σᶻ_i⟩_GS = 0 for OBC (odd-product expectation in Gaussian state).
_sz_expect(::AbstractMatrix, ::Int) = 0.0

# ═══════════════════════════════════════════════════════════════════════════════
# Spreading-correlation convenience
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _sz_sz_spreading(N, J, h, center, times) -> Matrix{ComplexF64}

Return `C[it, ix] = ⟨σᶻ_ix(t_it) σᶻ_center(0)⟩_GS` for `ix ∈ 1:N`,
`t_it ∈ 1:length(times)`.  Reuses the Majorana Hamiltonian, covariance,
and (per time-step) evolution matrix, so the cost is
`O(length(times) · N · M³)` with `M = 2(center + N) - 2`.
"""
function _sz_sz_spreading(
    N::Int, J::Float64, h::Float64, center::Int, times::AbstractVector{<:Real}
)
    (1 ≤ center ≤ N) || throw(ArgumentError("center site out of range"))
    hmat = _majorana_ham(N, J, h)
    Σ = _majorana_covariance_gs(hmat)
    nt = length(times)
    C = zeros(ComplexF64, nt, N)
    for (it, t) in enumerate(times)
        R = _majorana_evolution(hmat, t)
        for ix in 1:N
            C[it, ix] = _sz_sz_corr_from_cached(Σ, R, ix, center)
        end
    end
    return C
end

# ═══════════════════════════════════════════════════════════════════════════════
# fetch dispatch
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::Model{:TFIM}, ::Quantity{:sz_sz_correlation}, ::OBC;
          i::Int, j::Int, t::Float64) -> ComplexF64

Exact ground-state `⟨σᶻ_i(t) σᶻ_j(0)⟩` for the OBC TFIM.
"""
function fetch(
    model::Model{:TFIM},
    ::Quantity{:sz_sz_correlation},
    ::OBC;
    i::Int,
    j::Int,
    t::Float64,
    kwargs...,
)
    N = Int(model.params[:N])
    J = Float64(model.params[:J])
    h = Float64(model.params[:h])
    return _sz_sz_corr(N, J, h, i, j, t)
end

"""
    fetch(model::Model{:TFIM}, ::Quantity{:sx_sx_correlation}, ::OBC;
          i::Int, j::Int, t::Float64) -> ComplexF64

Exact ground-state `⟨σˣ_i(t) σˣ_j(0)⟩` for the OBC TFIM.
"""
function fetch(
    model::Model{:TFIM},
    ::Quantity{:sx_sx_correlation},
    ::OBC;
    i::Int,
    j::Int,
    t::Float64,
    kwargs...,
)
    N = Int(model.params[:N])
    J = Float64(model.params[:J])
    h = Float64(model.params[:h])
    return _sx_sx_corr(N, J, h, i, j, t)
end

"""
    fetch(model::Model{:TFIM}, ::Quantity{:sz_sz_spreading}, ::OBC;
          center::Int, times::AbstractVector{<:Real}) -> Matrix{ComplexF64}

Exact ground-state spreading correlation `C[it, ix] = ⟨σᶻ_ix(t_it) σᶻ_center(0)⟩`
for all sites `ix ∈ 1:N` and all `t_it ∈ times`.
"""
function fetch(
    model::Model{:TFIM},
    ::Quantity{:sz_sz_spreading},
    ::OBC;
    center::Int,
    times::AbstractVector{<:Real},
    kwargs...,
)
    N = Int(model.params[:N])
    J = Float64(model.params[:J])
    h = Float64(model.params[:h])
    return _sz_sz_spreading(N, J, h, center, times)
end
